!*****************************************************************************************
!>
!  Splined lunar principal-axis (MOON_PA) frame: the rotation from J2000 to the
!  frame, interpolated from a table of the frame's roll, pitch and yaw angles vs.
!  ephemeris time (e.g., `data/moon_pa_2000_2100.csv`, generated with SPICE).
!
!  This is the same interpolation as the `moonframe` library's
!  `moon_frame_interpolater` (quartic B-splines of the roll, the pitch, and the
!  cosine and sine of the yaw angle, with the same knots), and gives identical
!  results, but with the same pattern as [[splined_ephemeris_module]]: the spline
!  coefficients are global variables, computed once and then only read, so all
!  the segments (and threads) share one copy. Each [[moon_pa_frame]] object only
!  holds its own small evaluation workspace (the interval indices and a work
!  array, which `db1val` modifies), so each segment (thread) has its own.
!  (With a copy of the coefficients in every segment, the memory grew with the
!  square of the number of revolutions.)

    module splined_moon_frame_module

    use bspline_module, only: db1ink, db1val, bspline_order_quartic
    use csv_module, only: csv_file
    use halo_kinds_module

    implicit none

    private

    integer,parameter :: kx = bspline_order_quartic  !! spline order (as in `moonframe`)
    integer,parameter :: iknot = 0  !! knots chosen by `db1ink` (as in `moonframe`)
    integer,parameter :: n_splines = 4  !! roll, pitch, cos(yaw), sin(yaw)

    ! the spline data: global variables, so there is only one copy
    ! (set by [[initialize_moon_pa_frame]], then only read):
    integer :: nx = 0  !! number of points
    real(wp),dimension(:,:),allocatable :: tx     !! `(nx+kx, n_splines)` knots
    real(wp),dimension(:,:),allocatable :: bcoef  !! `(nx, n_splines)` coefficients

    type,public :: moon_pa_frame
        !! The interface to the splined MOON_PA frame: only the evaluation workspace
        !! (the splines themselves are global). Make a copy for each thread.
        private
        integer,dimension(n_splines) :: inbvx = 1  !! `db1val`'s interval index for each spline
        real(wp),dimension(3*kx) :: w0 = 0.0_wp    !! `db1val` work array
    contains
        private
        procedure,public :: initialize => initialize_moon_pa_frame
        procedure,public :: j2000_to_frame
        procedure,public :: frame_to_j2000
        procedure,public :: destroy => destroy_moon_pa_frame
    end type moon_pa_frame

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the table and compute the (global) splines, for the time span `[et0, etf]`:
!  the rows from the last one at or before `et0` to the first one at or after
!  `etf` (as in `moonframe`). Then the object's workspace is reset; copies of
!  the object (`allocate(..., source=...)`) can be used by other threads.

    subroutine initialize_moon_pa_frame(me, filename, et0, etf)

    class(moon_pa_frame),intent(inout) :: me
    character(len=*),intent(in) :: filename !! csv file with the ephemeris time, roll, pitch and yaw (rad)
    real(wp),intent(in) :: et0 !! start ephemeris time of the span to spline [sec]
    real(wp),intent(in) :: etf !! end ephemeris time of the span to spline [sec]

    type(csv_file) :: f
    real(wp),dimension(:),allocatable :: et, roll, pitch, yaw, fcn
    integer :: n, istart, iend, i, j, iflag
    logical :: status_ok

    ! read the data from the csv file:
    call f%read(filename, header_row=1, status_ok=status_ok)
    if (.not. status_ok) error stop 'error reading file: '//trim(filename)
    call f%get(1, et,    status_ok); if (.not. status_ok) error stop 'error getting et from file: '//trim(filename)
    call f%get(2, roll,  status_ok); if (.not. status_ok) error stop 'error getting roll from file: '//trim(filename)
    call f%get(3, pitch, status_ok); if (.not. status_ok) error stop 'error getting pitch from file: '//trim(filename)
    call f%get(4, yaw,   status_ok); if (.not. status_ok) error stop 'error getting yaw from file: '//trim(filename)
    call f%destroy()

    ! the rows to use (the span, with the rows just outside it):
    n = size(et)
    istart = minloc(abs(et-et0), 1)
    if (et(istart) > et0) istart = max(1, istart - 1)
    iend = minloc(abs(et-etf), 1)
    if (et(iend) < etf) iend = min(n, iend + 1)
    if (et(iend) < et(istart)) error stop 'Error: end time is before start time'

    ! the global splines. roll and pitch have no discontinuities; the yaw angle wraps
    ! around once per lunar rotation, so its cosine and sine are splined instead.
    if (allocated(tx)) deallocate(tx)
    if (allocated(bcoef)) deallocate(bcoef)
    nx = iend - istart + 1
    allocate(tx(nx+kx, n_splines), bcoef(nx, n_splines), fcn(nx))
    do j = 1, n_splines
        do i = istart, iend  ! (a loop: no large temporaries on the stack)
            select case (j)
            case (1); fcn(i-istart+1) = roll(i)
            case (2); fcn(i-istart+1) = pitch(i)
            case (3); fcn(i-istart+1) = cos(yaw(i))
            case (4); fcn(i-istart+1) = sin(yaw(i))
            end select
        end do
        call db1ink(et(istart:iend), nx, fcn, kx, iknot, tx(:,j), bcoef(:,j), iflag)
        if (iflag/=0) error stop 'error initializing the MOON_PA spline: '//trim(filename)
    end do

    call me%destroy()  ! (reset the workspace)

    end subroutine initialize_moon_pa_frame
!*****************************************************************************************

!*****************************************************************************************
!>
!  The rotation matrix from J2000 to the MOON_PA frame at `et`.

    function j2000_to_frame(me, et) result(rot)

    class(moon_pa_frame),intent(inout) :: me
    real(wp),intent(in) :: et  !! ephemeris time [sec]
    real(wp),dimension(3,3) :: rot

    real(wp),dimension(n_splines) :: v  !! roll, pitch, cos(yaw), sin(yaw)
    integer :: j, iflag

    if (nx == 0) error stop 'error: the MOON_PA frame has not been initialized'
    do j = 1, n_splines
        call db1val(et, 0, tx(:,j), nx, kx, bcoef(:,j), v(j), iflag, me%inbvx(j), me%w0, extrap=.false.)
        if (iflag/=0) error stop 'error evaluating the MOON_PA spline.'
    end do
    call rpy_to_rot(v(1), v(2), atan2(v(4), v(3)), rot)

    end function j2000_to_frame
!*****************************************************************************************

!*****************************************************************************************
!>
!  The rotation matrix from the MOON_PA frame to J2000 at `et`.

    function frame_to_j2000(me, et) result(rot)
    class(moon_pa_frame),intent(inout) :: me
    real(wp),intent(in) :: et  !! ephemeris time [sec]
    real(wp),dimension(3,3) :: rot
    rot = transpose(me%j2000_to_frame(et))
    end function frame_to_j2000
!*****************************************************************************************

!*****************************************************************************************
!>
!  Roll, pitch and yaw (rad) to the rotation matrix (as in `moonframe`).

    pure subroutine rpy_to_rot(roll, pitch, yaw, r)

    real(wp),intent(in)  :: roll, pitch, yaw
    real(wp),dimension(3,3),intent(out) :: r

    real(wp) :: cr, sr, cp, sp, cy, sy

    cr = cos(roll)
    sr = sin(roll)
    cp = cos(pitch)
    sp = sin(pitch)
    cy = cos(yaw)
    sy = sin(yaw)

    r(1,1) = cy*cp
    r(1,2) = cy*sp*sr - sy*cr
    r(1,3) = cy*sp*cr + sy*sr
    r(2,1) = sy*cp
    r(2,2) = sy*sp*sr + cy*cr
    r(2,3) = sy*sp*cr - cy*sr
    r(3,1) = -sp
    r(3,2) = cp*sr
    r(3,3) = cp*cr

    end subroutine rpy_to_rot
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reset the object's workspace (the global splines are kept).

    subroutine destroy_moon_pa_frame(me)
    class(moon_pa_frame),intent(inout) :: me
    me%inbvx = 1
    me%w0 = 0.0_wp
    end subroutine destroy_moon_pa_frame
!*****************************************************************************************

    end module splined_moon_frame_module
!*****************************************************************************************
