!*****************************************************************************************
!>
!  Tests for `astro_partials_module`:
!
!  * Pines acceleration vs. the toolkit's `geopotential_model_pines`,
!  * Pines gravity gradient vs. central differences of the acceleration,
!  * point-mass and third-body gravity gradients vs. central differences,
!  * the two-body rotating frame matrices vs. the toolkit's frame transformation,
!    and their time derivatives vs. central differences in time.
!
!  Run from the halo root directory (it reads `data/grav/gggrx_0020pm_sha.tab`
!  and `data/eph/JPLEPH.421`):  `fpm test partials_test --profile release`

    program partials_test

    use halo_kinds_module, only: wp
    use astro_partials_module
    use splined_ephemeris_module
    use fortran_astrodynamics_toolkit, only: geopotential_model_pines, third_body_gravity, &
                                             icrf_frame, two_body_rotating_frame, &
                                             center_at_secondary_body, body_earth, body_moon

    implicit none

    character(len=*),parameter :: gravfile = 'data/grav/gggrx_0020pm_sha.tab'
    character(len=*),parameter :: ephfile  = 'data/eph/JPLEPH.421'
    real(wp),parameter :: et_test = 820317600.0_wp  !! 2025-12-29 22:00:00 TDB

    type(geopotential_model_pines) :: grav
    type(jpl_ephemeris_splined) :: eph
    real(wp),dimension(3,5) :: rtest
    real(wp),dimension(3) :: a_fat, a_pines, ap, am, r, rb, a1, a2
    real(wp),dimension(3,3) :: grad, grad_fd
    real(wp),dimension(6,6) :: m, mdot, n, ndot, mp, mm, np, nm, m_fat
    real(wp),dimension(6) :: rv, rvdot, e, out
    real(wp) :: h, err, scale, worst_acc, worst_grad, worst_3b, worst_m, worst_mdot, worst_n, worst_ndot, et
    integer :: i, j, k, nmax
    logical :: status_ok, all_ok
    type(icrf_frame) :: inertial
    type(two_body_rotating_frame) :: rotating

    all_ok = .true.

    !-----------------------------------------------------------
    ! gravity: Pines acceleration and gradient (8x8, as in the baseline)
    !-----------------------------------------------------------
    nmax = 8
    call grav%initialize(gravfile, nmax, nmax, status_ok)
    if (.not. status_ok) error stop 'error reading gravity file'

    ! test points [km, Moon body-fixed]: NRHO periapsis/apoapsis-like distances, and at the poles
    rtest(:,1) = [ 1000.0_wp,  -500.0_wp, -3100.0_wp]
    rtest(:,2) = [ 0.0_wp,      0.0_wp,   -3300.0_wp]   ! over the south pole (a singular point for latitude/longitude methods)
    rtest(:,3) = [ 0.0_wp,      0.0_wp,    1900.0_wp]   ! over the north pole, low altitude
    rtest(:,4) = [ 1740.0_wp,  900.0_wp,   200.0_wp]    ! low altitude
    rtest(:,5) = [-20000.0_wp, 35000.0_wp, 60000.0_wp]  ! apoapsis distance

    worst_acc = 0.0_wp
    worst_grad = 0.0_wp
    do k = 1, size(rtest,2)
        r = rtest(:,k)
        call grav%get_acc(r, nmax, nmax, a_fat)
        call pines_gravity(r, nmax, nmax, grav%re, grav%mu, grav%cnm, grav%snm, a_pines, grad)
        worst_acc = max(worst_acc, norm2(a_pines - a_fat) / norm2(a_fat))
        ! 4th-order central differences of the acceleration:
        h = 1.0e-3_wp * norm2(r)
        do j = 1, 3
            r = rtest(:,k); r(j) = r(j) + h
            call pines_gravity(r, nmax, nmax, grav%re, grav%mu, grav%cnm, grav%snm, ap)
            r = rtest(:,k); r(j) = r(j) - h
            call pines_gravity(r, nmax, nmax, grav%re, grav%mu, grav%cnm, grav%snm, am)
            r = rtest(:,k); r(j) = r(j) + 2.0_wp*h
            call pines_gravity(r, nmax, nmax, grav%re, grav%mu, grav%cnm, grav%snm, a1)
            r = rtest(:,k); r(j) = r(j) - 2.0_wp*h
            call pines_gravity(r, nmax, nmax, grav%re, grav%mu, grav%cnm, grav%snm, a2)
            grad_fd(:,j) = (-a1 + 8.0_wp*ap - 8.0_wp*am + a2) / (12.0_wp*h)
        end do
        worst_grad = max(worst_grad, maxval(abs(grad - grad_fd)) / maxval(abs(grad)))
    end do
    write(*,'(A,ES10.2)') 'Pines acceleration vs toolkit gravpot (max rel. diff):      ', worst_acc
    write(*,'(A,ES10.2)') 'Pines gradient vs 4th-order diff. of acceleration (max rel.):', worst_grad
    call check(worst_acc  < 1.0e-13_wp, 'Pines acceleration')
    call check(worst_grad < 1.0e-9_wp,  'Pines gradient')  ! (limited by the difference's truncation error)

    !-----------------------------------------------------------
    ! third-body gravity gradient (Earth, from the Moon):
    !-----------------------------------------------------------
    rb = [300000.0_wp, -200000.0_wp, 50000.0_wp]
    worst_3b = 0.0_wp
    do k = 1, size(rtest,2)
        grad = third_body_gravity_gradient(rtest(:,k), rb, 398600.0_wp)
        h = 1.0_wp
        do j = 1, 3
            r = rtest(:,k); r(j) = r(j) + h
            call third_body_gravity(r, rb, 398600.0_wp, a1)
            r = rtest(:,k); r(j) = r(j) - h
            call third_body_gravity(r, rb, 398600.0_wp, a2)
            grad_fd(:,j) = (a1 - a2) / (2.0_wp*h)
        end do
        worst_3b = max(worst_3b, maxval(abs(grad - grad_fd)) / maxval(abs(grad)))
    end do
    write(*,'(A,ES10.2)') 'third-body gradient vs central diff. (max rel.):            ', worst_3b
    call check(worst_3b < 1.0e-6_wp, 'third-body gradient')

    !-----------------------------------------------------------
    ! rotating frame matrices:
    !-----------------------------------------------------------
    call eph%initialize_splinded_ephemeris(filename=ephfile, status_ok=status_ok, &
                                           et0=et_test-86400.0_wp*10, dt=3600.0_wp, etf=et_test+86400.0_wp*10)
    if (.not. status_ok) error stop 'error initializing the splined ephemeris'

    worst_m = 0.0_wp; worst_n = 0.0_wp; worst_mdot = 0.0_wp; worst_ndot = 0.0_wp
    do k = 0, 4
        et = et_test + k*86400.0_wp*1.37_wp
        call frame(et, m, mdot, n, ndot)
        ! M vs. the toolkit's transformation (applied to each unit vector):
        inertial = icrf_frame(b=body_moon)
        rotating = two_body_rotating_frame(primary_body=body_earth, secondary_body=body_moon, &
                                           center=center_at_secondary_body, et=et)
        do j = 1, 6
            e = 0.0_wp; e(j) = 1.0_wp
            call inertial%transform(e, rotating, et, eph, out, status_ok)
            m_fat(:,j) = out
        end do
        worst_m = max(worst_m, maxval(abs(m - m_fat)) / maxval(abs(m_fat)))
        worst_n = max(worst_n, maxval(abs(matmul(m,n) - identity6())))
        ! time derivatives vs. central differences:
        h = 60.0_wp
        call frame(et+h, mp, n=np)
        call frame(et-h, mm, n=nm)
        scale = maxval(abs(mdot))
        worst_mdot = max(worst_mdot, maxval(abs(mdot - (mp-mm)/(2.0_wp*h))) / scale)
        scale = maxval(abs(ndot))
        worst_ndot = max(worst_ndot, maxval(abs(ndot - (np-nm)/(2.0_wp*h))) / scale)
    end do
    write(*,'(A,ES10.2)') 'frame M vs toolkit transform (max rel. diff):               ', worst_m
    write(*,'(A,ES10.2)') 'frame M*N - I (max abs.):                                   ', worst_n
    write(*,'(A,ES10.2)') 'frame dM/dt vs central diff. in time (max rel.):            ', worst_mdot
    write(*,'(A,ES10.2)') 'frame dN/dt vs central diff. in time (max rel.):            ', worst_ndot
    call check(worst_m    < 1.0e-12_wp, 'frame M')
    call check(worst_n    < 1.0e-12_wp, 'frame N = inverse of M')
    call check(worst_mdot < 1.0e-6_wp,  'frame dM/dt')
    call check(worst_ndot < 1.0e-6_wp,  'frame dN/dt')

    if (all_ok) then
        write(*,'(A)') 'all partials tests passed'
    else
        error stop 'some partials tests FAILED'
    end if

    contains

        subroutine frame(et, m, mdot, n, ndot)
        real(wp),intent(in) :: et
        real(wp),dimension(6,6),intent(out),optional :: m, mdot, n, ndot
        call eph%get_rv(et, body_moon, body_earth, rv, status_ok)
        call eph%get_rv_derivative(et, body_moon, body_earth, rvdot, status_ok)
        call two_body_rotating_frame_matrices(rv, rvdot(1:3), rvdot(4:6), m, mdot, n, ndot)
        end subroutine frame

        pure function identity6() result(x)
        real(wp),dimension(6,6) :: x
        integer :: i
        x = 0.0_wp
        do i = 1, 6
            x(i,i) = 1.0_wp
        end do
        end function identity6

        subroutine check(ok, name)
        logical,intent(in) :: ok
        character(len=*),intent(in) :: name
        if (.not. ok) then
            write(*,'(A)') 'FAILED: '//name
            all_ok = .false.
        end if
        end subroutine check

    end program partials_test
!*****************************************************************************************
