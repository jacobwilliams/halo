!*****************************************************************************************
!>
!  A reference analytic Jacobian at a given point, for the Jacobian accuracy tests: meant to
!  be built in quad precision (`-DREAL128`), with tight integration tolerances.
!
!  Usage (from the halo root directory):
!
!    fpm test jacobian_ref_test --profile release --flag "-fopenmp -DREAL128" -- config.json n_revs x_file jac_file
!
!  `x_file` is written by `jacobian_test.f90` (its fifth argument): the point `x` with the
!  `xoffset`, `xscale` and `fscale` of that (double-precision) run. They replace this run's own,
!  so the Jacobian is at exactly the same point and scaled the same way: the epoch variables
!  are offsets from the guess's epochs, which differ at round-off level between builds.
!
!  `jac_file` gets `nnz` (int32) and the `nnz` values of the scaled analytic Jacobian (float64),
!  in the sparsity pattern order.

    program jacobian_ref_test

    use halo_module
    use halo_kinds_module, only: wp
    use iso_fortran_env, only: output_unit, int32, real64
!$  use omp_lib

    implicit none

    type(my_solver_type) :: solver
    real(wp),dimension(:),allocatable :: x, jac
    real(real64),dimension(:),allocatable :: x_in, xoffset_in, xscale_in, fscale_in
    integer,dimension(:),allocatable :: irow, icol
    character(len=1024) :: arg
    integer :: n_revs, iunit
    integer(int32) :: n, m
    real(wp) :: t0

    if (command_argument_count() < 4) error stop 'usage: jacobian_ref_test config.json n_revs x_file jac_file'

    call get_command_argument(1, arg)
    call solver%read_config_file(trim(arg))
    call get_command_argument(2, arg)
    read(arg,*) n_revs
    solver%mission%n_revs = n_revs
    solver%mission%jacobian_mode = 2
    solver%mission%use_splined_ephemeris = .true.
    call solver%mission%init(x)

    ! the point and the scales of the double-precision run:
    call get_command_argument(3, arg)
    open(newunit=iunit, file=trim(arg), access='stream', form='unformatted', status='old')
    read(iunit) n
    if (n /= size(x)) error stop 'the x file has a different problem size'
    allocate(x_in(n), xoffset_in(n), xscale_in(n))
    read(iunit) x_in, xoffset_in, xscale_in
    read(iunit) m
    if (m /= size(solver%mission%fscale)) error stop 'the x file has a different number of constraints'
    allocate(fscale_in(m))
    read(iunit) fscale_in
    close(iunit)

    write(output_unit,'(A,I0,A,I0)') ' n: ', n, '   m: ', m
    write(output_unit,'(A,ES10.2)') ' this run''s own guess vs the given x:  max |dx| (scaled)       = ', &
                                    maxval(abs(x - real(x_in,wp)))
    write(output_unit,'(A,ES10.2)') ' this run''s own xoffset vs the given:   max |d xoffset| (days)  = ', &
                                    maxval(abs(solver%mission%xoffset - real(xoffset_in,wp)))
    write(output_unit,'(A,ES10.2)') ' this run''s own xscale vs the given:    max rel. difference     = ', &
                                    maxval(abs(solver%mission%xscale/real(xscale_in,wp) - 1.0_wp))
    write(output_unit,'(A,ES10.2)') ' this run''s own fscale vs the given:    max rel. difference     = ', &
                                    maxval(abs(solver%mission%fscale/real(fscale_in,wp) - 1.0_wp))
    x                       = real(x_in, wp)
    solver%mission%xoffset  = real(xoffset_in, wp)
    solver%mission%xscale   = real(xscale_in, wp)
    solver%mission%fscale   = real(fscale_in, wp)

    t0 = wall_time()
    call solver%mission%compute_jacobian_stm(x, jac)
    write(output_unit,'(A,F10.2,A)') ' analytic Jacobian time: ', wall_time() - t0, ' s'

    call solver%mission%get_sparsity_pattern(irow, icol)
    if (size(jac) /= size(irow)) error stop 'unexpected Jacobian size'
    call get_command_argument(4, arg)
    open(newunit=iunit, file=trim(arg), access='stream', form='unformatted', status='replace')
    write(iunit) int(size(jac),int32)
    write(iunit) real(jac, real64)
    close(iunit)

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

    end program jacobian_ref_test
!*****************************************************************************************
