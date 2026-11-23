!*****************************************************************************************
!>
!  Compare the two Jacobian methods of the halo problem at the initial guess:
!  finite differences (`jacobian_mode=1`) and the analytic state transition
!  matrix partials (`jacobian_mode=2`).
!
!  Usage (from the halo root directory):
!
!    fpm test jacobian_test --profile release --flag "-fopenmp" -- [config.json] [n_revs] [jac_file] [ref_config.json]
!
!  With `jac_file`, the two Jacobians (scaled, in the sparsity pattern order) are also
!  written to that file (stream access: `nnz` (int32), then the `nnz` finite-difference values,
!  then the `nnz` STM values (float64)), e.g., to compare runs with other tolerances.
!
!  With `ref_config.json` (e.g., the same problem with tighter integration tolerances), a
!  reference analytic Jacobian is computed with that config *at the same variables `x`*,
!  and appended to `jac_file` (another `nnz` float64 values). The initial guess depends
!  slightly on the integration tolerances (the guess generation integrates), so Jacobians
!  from separate runs are at slightly different points, and their differences then
!  include that, not only the integration error. Comparing at the same `x` avoids it.
!
!  The default config file is `examples/example_sparse.json`. The problem is set
!  up with `jacobian_mode=2` (the finite differences are also available then).

    program jacobian_test

    use halo_module
    use halo_kinds_module, only: wp
    use iso_fortran_env, only: output_unit, int32
!$  use omp_lib

    implicit none

    type(my_solver_type) :: solver, ref_solver
    real(wp),dimension(:),allocatable :: x, x_ref, jac_fd, jac_stm, jac_ref, colnorm
    integer,dimension(:),allocatable :: irow, icol, keep
    integer,dimension(:),allocatable :: removed
    character(len=1024) :: arg
    character(len=:),allocatable :: config
    integer :: i, k, n, n_full, n_revs, fc, nk
    real(wp) :: t0, t1, t_fd, t_stm, d, err_state, err_epoch, num_state, num_epoch, &
                worst_state, worst_epoch, rel
    logical :: is_epoch
    integer :: iunit

    config = 'examples/example_sparse.json'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        config = trim(arg)
    end if

    ! read the config, and set up the problem with the analytic partials:
    call solver%read_config_file(config)
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg,*) n_revs
        solver%mission%n_revs = n_revs
    end if
    solver%mission%jacobian_mode = 2
    solver%mission%use_splined_ephemeris = .true.
    call solver%mission%init(x)

    ! finite differences:
    do i = 1, size(solver%mission%segs)
        call solver%mission%segs(i)%cache()
    end do
    t0 = wall_time()
    call solver%mission%compute_jacobian(x, jac_fd)
    t1 = wall_time()
    t_fd = t1 - t0
    do i = 1, size(solver%mission%segs)
        call solver%mission%segs(i)%uncache()
    end do

    ! analytic:
    t0 = wall_time()
    call solver%mission%compute_jacobian_stm(x, jac_stm)
    t1 = wall_time()
    t_stm = t1 - t0

    ! compare, column by column (scaled Jacobian): the error relative to the column's norm,
    ! separately for the epoch columns and the state columns
    call solver%mission%get_sparsity_pattern(irow, icol)
    call solver%mission%define_problem_size(n=n)
    call solver%mission%define_problem_size(n=n_full, full_problem=.true.)
    removed = solver%mission%removed_columns()
    allocate(keep(n))
    nk = 0
    do i = 1, n_full
        if (any(removed==i)) cycle
        nk = nk + 1
        keep(nk) = i
    end do
    allocate(colnorm(n))
    colnorm = 0.0_wp
    do k = 1, size(irow)
        colnorm(icol(k)) = colnorm(icol(k)) + jac_fd(k)**2
    end do
    colnorm = sqrt(colnorm)

    err_state = 0.0_wp; num_state = 0.0_wp; worst_state = 0.0_wp
    err_epoch = 0.0_wp; num_epoch = 0.0_wp; worst_epoch = 0.0_wp
    do k = 1, size(irow)
        fc = keep(icol(k))
        is_epoch = mod(fc-1,7) == 0
        d = jac_stm(k) - jac_fd(k)
        rel = abs(d) / max(colnorm(icol(k)), tiny(1.0_wp))
        if (is_epoch) then
            err_epoch = err_epoch + d**2
            num_epoch = num_epoch + jac_fd(k)**2
            worst_epoch = max(worst_epoch, rel)
        else
            err_state = err_state + d**2
            num_state = num_state + jac_fd(k)**2
            worst_state = max(worst_state, rel)
        end if
    end do

    write(output_unit,'(A)') ''
    write(output_unit,'(A,I0,A,I0,A,I0)') ' revs: ', solver%mission%n_revs, '   n: ', n, '   nonzeros: ', size(irow)
    write(output_unit,'(A,F10.4,A)') ' finite-difference Jacobian time: ', t_fd,  ' s'
    write(output_unit,'(A,F10.4,A)') ' STM Jacobian time:               ', t_stm, ' s'
    write(output_unit,'(A,F10.2)')   ' speedup:                         ', t_fd/t_stm
    write(output_unit,'(A)') ' difference (STM - finite differences), scaled Jacobian:'
    write(output_unit,'(A,ES10.2,A,ES10.2)') '   state columns: ||dJ||/||J|| = ', sqrt(err_state/num_state), &
                                             '   max |dJ|/||column|| = ', worst_state
    write(output_unit,'(A,ES10.2,A,ES10.2)') '   epoch columns: ||dJ||/||J|| = ', sqrt(err_epoch/num_epoch), &
                                             '   max |dJ|/||column|| = ', worst_epoch

    ! the reference Jacobian (another config, at the same x):
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg)
        call ref_solver%read_config_file(trim(arg))
        ref_solver%mission%n_revs = solver%mission%n_revs
        ref_solver%mission%jacobian_mode = 2
        ref_solver%mission%use_splined_ephemeris = .true.
        call ref_solver%mission%init(x_ref)
        if (size(x_ref) /= size(x)) error stop 'the reference config has a different problem size'
        call ref_solver%mission%compute_jacobian_stm(x, jac_ref)  ! (at this run's x)
        write(output_unit,'(A)') ' relative to the reference (analytic, '//trim(arg)//', same x):'
        write(output_unit,'(A,ES10.2)') '   finite differences: ||dJ||/||J_ref|| = ', norm2(jac_fd-jac_ref)/norm2(jac_ref)
        write(output_unit,'(A,ES10.2)') '   analytic:           ||dJ||/||J_ref|| = ', norm2(jac_stm-jac_ref)/norm2(jac_ref)
        write(output_unit,'(A,ES10.2,A)') '   (the reference config''s own guess differs by max |x - x_ref| = ', &
                                          maxval(abs(x-x_ref)), ' (scaled), not used)'
    end if

    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg)
        open(newunit=iunit, file=trim(arg), access='stream', form='unformatted', status='replace')
        write(iunit) int(size(irow),int32)
        write(iunit) jac_fd
        write(iunit) jac_stm
        if (allocated(jac_ref)) write(iunit) jac_ref
        close(iunit)
    end if

    contains

        function wall_time() result(t)
        real(wp) :: t
        integer(8) :: count, rate
        t = 0.0_wp
!$      t = omp_get_wtime()
!$      return
        call system_clock(count, rate)
        t = real(count,wp) / real(rate,wp)
        end function wall_time

    end program jacobian_test
!*****************************************************************************************
