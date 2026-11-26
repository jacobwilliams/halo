!*****************************************************************************************
!>
!  Write the scaled analytic Jacobian of the continuity constraints at a point, for
!  analysis outside the program (e.g., its singular values).
!
!  Usage (from the halo root directory):
!
!    fpm test jacobian_dump_test --profile release --flag "-fopenmp" -- config.json out_file
!
!  The point is the config's initial guess (use `initial_guess_from_file` for a converged
!  solution). `out_file` gets (stream access): `m`, `n`, `nnz` (int32), then `irow` and `icol`
!  (`nnz` int32 each), the values (`nnz` float64), and the full-problem column of each column
!  of `x` (`n` int32: node `(c-1)/7+1`, component `mod(c-1,7)`: 0 the epoch, 1-6 the state).

    program jacobian_dump_test

    use halo_module
    use halo_kinds_module, only: wp
    use iso_fortran_env, only: output_unit, int32, real64

    implicit none

    type(my_solver_type) :: solver
    real(wp),dimension(:),allocatable :: x, jac
    integer,dimension(:),allocatable :: irow, icol, removed, keep
    character(len=1024) :: config, out_file
    integer :: n, m, n_full, i, k, iunit

    if (command_argument_count() < 2) then
        ! (a tool, run by hand: under a plain `fpm test` there are no arguments, and it is skipped)
        write(*,'(A)') 'usage: jacobian_dump_test config.json out_file (no arguments: skipped)'
        stop
    end if
    call get_command_argument(1, config)
    call get_command_argument(2, out_file)

    call solver%read_config_file(trim(config))
    solver%mission%jacobian_mode = 2
    solver%mission%use_splined_ephemeris = .true.
    call solver%mission%init(x)
    call solver%mission%define_problem_size(n=n, m=m)
    call solver%mission%define_problem_size(n=n_full, full_problem=.true.)
    call solver%mission%compute_jacobian_stm(x, jac)
    call solver%mission%get_sparsity_pattern(irow, icol)

    removed = solver%mission%removed_columns()
    allocate(keep(n))
    k = 0
    do i = 1, n_full
        if (any(removed==i)) cycle
        k = k + 1
        keep(k) = i
    end do

    open(newunit=iunit, file=trim(out_file), access='stream', form='unformatted', status='replace')
    write(iunit) int(m,int32), int(n,int32), int(size(irow),int32)
    write(iunit) int(irow,int32), int(icol,int32)
    write(iunit) real(jac,real64)
    write(iunit) int(keep,int32)
    close(iunit)
    write(output_unit,'(A,I0,A,I0,A,I0,A)') ' wrote the Jacobian: ', m, ' x ', n, ', ', size(irow), ' nonzeros'

    end program jacobian_dump_test
!*****************************************************************************************
