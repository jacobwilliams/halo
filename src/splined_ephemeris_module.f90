!*****************************************************************************************
!>
!  Simple splined ephemeris model for Earth/Sun wrt Moon.

    module splined_ephemeris_module

    use fortran_astrodynamics_toolkit
    use bspline_module, only: db1ink, db1val, bspline_order_cubic
    use iso_fortran_env
    use halo_kinds_module

    implicit none

    private

        ! spline parameters:
        integer,parameter :: kx = bspline_order_cubic !! spline interpolation order
        integer,parameter :: iknot = 0 !! use default knots

        !TODO: this should be moved into the fortran astrodynamics toolkit
        ! [also update the MU .... it isn't used here]
        type(celestial_body),parameter,public :: body_ssb = celestial_body(0, 'SSB', 0.0_wp ) !! solar-system barycenter

        ! the coefficients are stored as global variables,
        ! so all the segments can access them without having
        ! to have multiple copies of the ephemeris in each segment
        type :: body_eph
            integer :: id = 0
            real(wp),dimension(:,:),allocatable :: f !! function evals. not needed once the spline has been computed
            real(wp),dimension(:,:),allocatable :: tx
            real(wp),dimension(:,:),allocatable :: bcoef
            integer :: nx = 0
            contains
            procedure,public :: destroy => destroy_body_eph
        end type body_eph
        type(body_eph),target :: earth_eph
        type(body_eph),target :: sun_eph
        type(body_eph),target :: ssb_eph

        type :: body_eph_interface
            !! the interface to a splined ephemeris for a body
            !! used for calls to `db1val`.
            private
            integer :: inbvx = 0
            real(wp),dimension(:),allocatable :: w0  !! work array - dimension(3_ip*kx)
            type(body_eph),pointer :: eph => null() !! points to the ephemeris
            contains
            procedure,public :: get_r
            procedure,public :: get_rv
            procedure,public :: destroy => destroy_body_eph_interface
        end type body_eph_interface

        type,extends(jpl_ephemeris),public :: jpl_ephemeris_splined
            !! make this an extension of the [[jpl_ephemeris]],
            !! since it is also needed in tranformations.
            !! also, has the extra feature of a `get_r` method,
            !! since we don't need velocity for the gravity calculation.

            type(body_eph_interface) :: earth_eph_interface !! splined version of earth ephemeris
            type(body_eph_interface) :: sun_eph_interface !! splined version of sun ephemeris
            type(body_eph_interface) :: ssb_eph_interface !! splined version of ssb ephemeris
        contains
            procedure,public :: initialize_splinded_ephemeris
            procedure :: initialize_globals !! this is done once to initialize the global ephemeris vars
            procedure,public :: get_r  => get_r_splined
            procedure,public :: get_rv => get_rv_splined
        end type jpl_ephemeris_splined

    contains

    !TODO: really we also need a routine to destroy earth_eph and sun_eph.

    subroutine destroy_body_eph(me)
        class(body_eph),intent(out) :: me
    end subroutine destroy_body_eph

    subroutine destroy_body_eph_interface(me)
        class(body_eph_interface),intent(out) :: me
    end subroutine destroy_body_eph_interface

    subroutine initialize_globals(me, et0, dt, etf)

        !! initialize the earth and sun ephemeris.

        class(jpl_ephemeris_splined),intent(inout) :: me
        real(wp),intent(in) :: et0 !! initial ephemeris time [sec]
        real(wp),intent(in) :: dt  !! ephemeris time step [sec]
        real(wp),intent(in) :: etf !! final ephemeris time [sec]

        integer :: i !! counter
        real(wp) :: et !! ephemeris time
        integer :: nx !! number of points
        real(wp),dimension(6) :: rv_earth, rv_sun
        logical :: status_ok
        real(wp),dimension(:),allocatable :: et_vec
        integer :: iflag

        ! write(*,*) 'initialize_globals...'
        ! write(*,*) '  et0 = ', et0
        ! write(*,*) '  dt  = ', dt
        ! write(*,*) '  etf = ', etf

        call me%earth_eph_interface%destroy()
        call me%sun_eph_interface%destroy()
        call me%ssb_eph_interface%destroy()

        if (allocated(earth_eph%tx)) deallocate(earth_eph%tx)
        if (allocated(sun_eph%tx)) deallocate(sun_eph%tx)
        if (allocated(ssb_eph%tx)) deallocate(ssb_eph%tx)
        if (allocated(earth_eph%bcoef)) deallocate(earth_eph%bcoef)
        if (allocated(sun_eph%bcoef)) deallocate(sun_eph%bcoef)
        if (allocated(ssb_eph%bcoef)) deallocate(ssb_eph%bcoef)
        if (allocated(earth_eph%f)) deallocate(earth_eph%f)
        if (allocated(sun_eph%f)) deallocate(sun_eph%f)
        if (allocated(ssb_eph%f)) deallocate(ssb_eph%f)
        if (allocated(et_vec)) deallocate(et_vec)

        ! first, count the number of points and allocate the arrays:
        nx = 0
        et = et0
        do
            if (et>etf) exit
            nx = nx + 1
            et = et + dt
        end do
        allocate(earth_eph%tx(   nx+kx,6))  ! columns are for each state element (rx,ry,rz,vx,vy,vz)
        allocate(sun_eph%tx(     nx+kx,6))
        allocate(ssb_eph%tx(     nx+kx,6))
        allocate(earth_eph%bcoef(nx,6))
        allocate(sun_eph%bcoef(  nx,6))
        allocate(ssb_eph%bcoef(  nx,6))
        allocate(earth_eph%f(    nx,6))
        allocate(sun_eph%f(      nx,6))
        allocate(ssb_eph%f(      nx,6))
        allocate(et_vec(         nx))

        ! function calls from et0 to etf:
        i = 0 ! index in the arrays
        et = et0
        do
            if (et>etf) exit
            i = i + 1
            et = et + dt
            et_vec(i) = et
            ! get_rv_from_jpl_ephemeris(me,et,targ,obs,rv,status_ok)
            call me%jpl_ephemeris%get_rv(et,body_earth,body_moon,earth_eph%f(i,:),status_ok)
            if (.not. status_ok) error stop 'earth eph error'
            call me%jpl_ephemeris%get_rv(et,body_sun,body_moon,sun_eph%f(i,:),status_ok)
            if (.not. status_ok) error stop 'sun eph error'
            call me%jpl_ephemeris%get_rv(et,body_sun,body_moon,ssb_eph%f(i,:),status_ok)
            if (.not. status_ok) error stop 'ssb eph error'
        end do

        ! create the splines (one for each coordinate):
        do i = 1, 6
            call db1ink(et_vec, nx, earth_eph%f(:,i), kx, iknot, earth_eph%tx(:,i), earth_eph%bcoef(:,i), iflag)
            if (iflag/=0) then
                write(*,*) 'db1ink iflag = ', iflag
                error stop 'db1ink error'
            end if
            call db1ink(et_vec, nx, sun_eph%f(:,i), kx, iknot, sun_eph%tx(:,i), sun_eph%bcoef(:,i), iflag)
            if (iflag/=0) then
                write(*,*) 'db1ink iflag = ', iflag
                error stop 'db1ink error'
            end if
            call db1ink(et_vec, nx, ssb_eph%f(:,i), kx, iknot, ssb_eph%tx(:,i), ssb_eph%bcoef(:,i), iflag)
            if (iflag/=0) then
                write(*,*) 'db1ink iflag = ', iflag
                error stop 'db1ink error'
            end if
        end do
        deallocate(earth_eph%f) ! don't need these anymore
        deallocate(sun_eph%f)
        deallocate(ssb_eph%f)

        earth_eph%nx = nx
        sun_eph%nx = nx
        ssb_eph%nx = nx

    end subroutine initialize_globals

    subroutine initialize_splinded_ephemeris(me,filename,ksize,km,bary,status_ok,et0,dt,etf)

        !! the ephemeris must be initialized before this is called

        class(jpl_ephemeris_splined),intent(inout) :: me
        character(len=*),intent(in)        :: filename   !! ephemeris file name
        integer,intent(in),optional        :: ksize      !! corresponding `ksize`
        logical,intent(in),optional        :: km         !! defining physical units of the output states.
                                                         !! `km = .true.`  : km and km/sec [default],
                                                         !! `km = .false.` : au and au/day.
        logical,intent(in),optional        :: bary       !! logical flag defining output center.
                                                         !! only the 9 planets are affected.
                                                         !! `bary = .true.`  : center is solar-system barycenter,
                                                         !! `bary = .false.` : center is sun [default].
        logical,intent(out) :: status_ok                 !! true if there were no problems.
        real(wp),intent(in) :: et0 !! initial ephemeris time [sec]
        real(wp),intent(in) :: dt  !! ephemeris time step [sec]
        real(wp),intent(in) :: etf !! final ephemeris time [sec]

        write(*,'(A)') ' * Using splined ephemeris'

        ! have to first initialize the base one:
        call me%jpl_ephemeris%initialize(filename,ksize,km,bary,status_ok)

        status_ok = (etf>et0 .and. dt>0.0_wp)
        if (.not. status_ok) return

        ! initialize the global spline coefficients
        call me%initialize_globals(et0, abs(dt), etf)

        ! now, the local variables in this class
        allocate(me%earth_eph_interface%w0(3*kx))
        allocate(me%sun_eph_interface%w0(3*kx))
        allocate(me%ssb_eph_interface%w0(3*kx))
        me%earth_eph_interface%inbvx = 0
        me%sun_eph_interface%inbvx = 0
        me%ssb_eph_interface%inbvx = 0

        ! point to the global ephemeris:
        me%earth_eph_interface%eph => earth_eph
        me%sun_eph_interface%eph => sun_eph
        me%ssb_eph_interface%eph => ssb_eph

    end subroutine initialize_splinded_ephemeris

    subroutine get_rv_splined(me,et,targ,obs,rv,status_ok)

        class(jpl_ephemeris_splined),intent(inout) :: me
        real(wp),intent(in)                :: et         !! ephemeris time [sec]
        type(celestial_body),intent(in)    :: targ       !! target body
        type(celestial_body),intent(in)    :: obs        !! observer body
        real(wp),dimension(6),intent(out)  :: rv         !! state of targ w.r.t. obs [km,km/s] in ICRF frame
        logical,intent(out)                :: status_ok  !! true if there were no problems

        status_ok = .true.
        if (targ==body_earth .and. obs==body_moon) then
            rv = me%earth_eph_interface%get_rv(et)
        elseif (targ==body_sun .and. obs==body_moon) then
            rv = me%sun_eph_interface%get_rv(et)
        elseif (targ==body_ssb .and. obs==body_moon) then
            rv = me%ssb_eph_interface%get_rv(et)

        elseif (targ==body_moon .and. obs==body_earth) then  ! inverse are negative
            rv = -me%earth_eph_interface%get_rv(et)
        elseif (targ==body_moon .and. obs==body_sun) then
            rv = -me%sun_eph_interface%get_rv(et)
        elseif (targ==body_moon .and. obs==body_ssb) then
            rv = -me%ssb_eph_interface%get_rv(et)

        else
            write(*,*) 'targ = ', targ
            write(*,*) 'obs  = ', obs
            error stop 'error in get_rv_splined: this combo has not been splined'
            ! or could call me%jpl_ephemeris%get_rv(et,targ,obs,rv,status_ok)
        end if

    end subroutine get_rv_splined

    subroutine get_r_splined(me,et,targ,obs,r,status_ok)

        class(jpl_ephemeris_splined),intent(inout) :: me
        real(wp),intent(in)                :: et         !! ephemeris time [sec]
        type(celestial_body),intent(in)    :: targ       !! target body
        type(celestial_body),intent(in)    :: obs        !! observer body
        real(wp),dimension(3),intent(out)  :: r          !! r of targ w.r.t. obs [km] in ICRF frame
        logical,intent(out)                :: status_ok  !! true if there were no problems

        status_ok = .true.
        if (targ==body_earth .and. obs==body_moon) then
            r = me%earth_eph_interface%get_r(et)
        elseif (targ==body_sun .and. obs==body_moon) then
            r = me%sun_eph_interface%get_r(et)
        elseif (targ==body_ssb .and. obs==body_moon) then
            r = me%ssb_eph_interface%get_r(et)

        elseif (targ==body_moon .and. obs==body_earth) then  ! inverse are negative
            r = -me%earth_eph_interface%get_r(et)
        elseif (targ==body_moon .and. obs==body_sun) then
            r = -me%sun_eph_interface%get_r(et)
        elseif (targ==body_moon .and. obs==body_ssb) then
            r = -me%ssb_eph_interface%get_r(et)

        elseif (targ==body_sun .and. obs==body_ssb) then
            ! for this one we subtract these
            ! ssb -> sun = ssb -> moon + moon -> sun
            r = me%sun_eph_interface%get_r(et) - me%ssb_eph_interface%get_r(et)

        elseif (targ==body_earth .and. obs==body_ssb) then
            ! for this one we subtract these
            ! ssb -> earth = ssb -> moon + moon -> earth
            r = me%earth_eph_interface%get_r(et) - me%ssb_eph_interface%get_r(et)

        else
            write(*,*) 'targ = ', targ
            write(*,*) 'obs  = ', obs
            error stop 'error in get_r_splined: this combo has not been splined'
            ! or could call me%jpl_ephemeris%get_rv(et,targ,obs,rv,status_ok); r = rv(1:3)
        end if

    end subroutine get_r_splined

    function get_rv(me,et) result(rv)

        class(body_eph_interface),intent(inout) :: me
        real(wp),intent(in) :: et !! ephemeris time (sec)
        real(wp),dimension(6) :: rv !! position/velocity vector

        integer :: iflag, i

        do i = 1, 6
            call db1val(et, 0, me%eph%tx(:,i), me%eph%nx, kx, me%eph%bcoef(:,i), rv(i), iflag, &
                        me%inbvx, me%w0, extrap=.false.)
            if (iflag /= 0) then
                write(*,*) 'et    = ', et
                write(*,*) 'iflag = ', iflag
                error stop 'error calling get_rv'
            end if
        end do

    end function get_rv

    function get_r(me,et) result(r)

        class(body_eph_interface),intent(inout) :: me
        real(wp),intent(in) :: et !! ephemeris time (sec)
        real(wp),dimension(3) :: r !! position vector

        integer :: iflag, i

        do i = 1, 3
            call db1val(et, 0, me%eph%tx(:,i), me%eph%nx, kx, me%eph%bcoef(:,i), r(i), iflag, &
                        me%inbvx, me%w0, extrap=.false.)
            if (iflag /= 0) then
                write(*,*) 'iflag = ', iflag
                error stop 'error calling get_r'
            end if
        end do

    end function get_r

    end module splined_ephemeris_module
!*****************************************************************************************
