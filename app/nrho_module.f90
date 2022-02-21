!*****************************************************************************************
!>
!  Main module with all the stuff to solve the NRHO problem.

    module nrho_module

    use parameters_module
    use fortran_astrodynamics_toolkit, fat_wp =>wp ! have to rename to avoid conflict [fix this]
    use ddeabm_module
    use nlesolver_module
    use numerical_differentiation_module
    use json_module
    use pyplot_module
    use csv_module

    implicit none

    public  ! just make it all public for now

    type segment_data

        !! data for a segment

        real(wp)              :: t0 = zero          !! initial time [days]
        real(wp)              :: t0_scale = one     !! opt var scale for `t0`

        real(wp),dimension(6) :: x0_rotating = zero !! initial state [km, km/s], rotating frame
        real(wp),dimension(6) :: x0_rotating_scale = one !! opt var scale for `x0_rotating`

        real(wp),dimension(6) :: x0 = zero          !! initial state [km, km/s], inertial frame
        real(wp)              :: tf = zero          !! final time [days]
        real(wp),dimension(6) :: xf = zero          !! final state [km, km/s], inertial frame

        real(wp),dimension(6) :: xf_rotating = zero !! final state [km, km/s], rotating frame
        real(wp),dimension(6) :: xf_rotating_scale = one !! func scale for `xf_rotating`

    end type segment_data

    type :: trajectory
        !! a segment trajectory (for plotting or output)
        real(wp),dimension(:),allocatable :: et
        real(wp),dimension(:),allocatable :: x
        real(wp),dimension(:),allocatable :: y
        real(wp),dimension(:),allocatable :: z
        real(wp),dimension(:),allocatable :: vx
        real(wp),dimension(:),allocatable :: vy
        real(wp),dimension(:),allocatable :: vz
        contains
        procedure :: destroy => destroy_trajectory
        procedure :: add => add_point_to_trajectory
    end type trajectory

    type,extends(ddeabm_class)  :: segment

        !! a ballistic segment in the mission
        !!
        !! * central body is the Moon (8x8 gravity),
        !!   with Earth and Sun as third-bodies.
        !! * Inputs are: `t0`, `tf', and `x0_rotating`.
        !! * Output are: `xf`, `xf_rotating`

        character(len=:),allocatable :: name  !! the segment name

        type(segment_data) :: data
        type(segment_data) :: cached_data  !! used when computing gradients

        ! note: originally I had grav and eph in the segments...
        !       ... but not sure about how that would work with coarrays...
        !
        !  These are pointers that are pointing to the ones in the mission ......

        type(geopotential_model_pines),pointer :: grav => null() !! central body geopotential model [global]
        type(jpl_ephemeris),pointer :: eph=> null()  !! the ephemeris [global]

        ! for saving the trajectory for plotting:
        type(trajectory) :: traj_inertial  !! in the inertial frame
        type(trajectory) :: traj_rotating  !! in the rotating frame

        contains

        procedure,public :: set_input   => set_segment_inputs
        procedure,public :: get_outputs => get_segment_outputs
        procedure,public :: set_outputs => set_segment_outputs
        procedure,public :: propagate   => propagate_segment

        procedure :: cache
        procedure :: uncache
        procedure :: put_data_in_segment

    end type segment

    type,extends(numdiff_type) :: mission_type

        !! the mission [this is a [[numdiff_type]] for
        !! organizational purposes.... rethink this maybe ...

        type(geopotential_model_pines),pointer :: grav => null()  !! central body geopotential model [global]
        type(jpl_ephemeris),pointer :: eph => null()                !! the ephemeris [global]
        type(segment),dimension(:),allocatable :: segs              !! the list of segments

        real(wp),dimension(:),allocatable :: xscale  !! opt var scale factors
        real(wp),dimension(:),allocatable :: fscale  !! function scale factors
        character(len=20),dimension(:),allocatable :: xname !! opt var names

    contains

        procedure,public :: init => initialize_the_mission
        procedure,public :: plot => plot_trajectory

        procedure,public :: write_optvars_to_file
        procedure,public :: constraint_violations
        procedure :: get_problem_arrays
        procedure :: put_x_in_segments
        procedure :: get_scales_from_segs
        procedure :: segs_to_propagate
        procedure :: get_sparsity_pattern

    end type mission_type

    type,extends(nlesolver_type) :: my_solver_type

        !! the solver class, which contains
        !! an instance of the mission

        type(mission_type) :: mission  !! the NRHO mission

    contains

        procedure,public :: init => initialize_the_solver

    end type my_solver_type

    ! test
    public :: define_problem_size

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destroy a trajectory (deallocate the data arrays).

    pure subroutine destroy_trajectory(me)

    implicit none

    class(trajectory),intent(inout) :: me

    if (allocated(me%et)) deallocate(me%et)
    if (allocated(me%x )) deallocate(me%x )
    if (allocated(me%y )) deallocate(me%y )
    if (allocated(me%z )) deallocate(me%z )
    if (allocated(me%vx)) deallocate(me%vx)
    if (allocated(me%vy)) deallocate(me%vy)
    if (allocated(me%vz)) deallocate(me%vz)

    end subroutine destroy_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the size variables for the "forward-backward" formulation of the problem.
!
!@note All the outputs are functions of `n_revs`,
!      the number of orbits we want to solve.

    pure subroutine define_problem_size(n,m,n_segs,n_nonzero)

    implicit none

    integer,intent(out),optional :: n          !! number of opt vars for the solver
    integer,intent(out),optional :: m          !! number of equality constraints for the solver
    integer,intent(out),optional :: n_segs     !! number of segments
    integer,intent(out),optional :: n_nonzero  !! number of nonzero elements in the Jacobian

    if (present(n))      n      = 28 * n_revs + 7
    if (present(m))      m      = 24 * n_revs
    if (present(n_segs)) n_segs = 8 * n_revs

    ! there are 4 blocks of nonzeros per rev
    ! (each block contains 84 elements)
    if (present(n_nonzero)) n_nonzero = n_revs * (84*4)

    end subroutine define_problem_size
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set *all* the data in a segment.

    subroutine put_data_in_segment(me,d)

    implicit none

    class(segment),intent(inout) :: me
    type(segment_data),intent(in) :: d

    me%data        = d
    me%cached_data = d ! also the cache

    end subroutine put_data_in_segment
!*****************************************************************************************

!*****************************************************************************************
!>
!  Cache all the data in a segment.

    subroutine cache(me)

    implicit none

    class(segment),intent(inout) :: me

    me%cached_data = me%data

    end subroutine cache
!*****************************************************************************************

!*****************************************************************************************
!>
!  Restore all the segment data form the cache.

    subroutine uncache(me)

    implicit none

    class(segment),intent(inout) :: me

    me%data = me%cached_data

    end subroutine uncache
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the mission.
!
!  This is the "forward-backward" formulation.
!  It is assumed that the mission has already been initialized.

    subroutine initialize_the_solver(me,config_file_name,x)

    implicit none

    class(my_solver_type),intent(inout) :: me
    character(len=*),intent(in) :: config_file_name !! the config file to read
    real(wp),dimension(:),allocatable,intent(out) :: x !! initial guess

    integer :: i          !! counter for index in `x`
    integer :: irev       !! counter for number of revs
    integer :: n          !! number of opt vars for the solver
    integer :: m          !! number of equality constraints for the solver
    logical :: status_ok  !! status flag for solver initialization
    integer :: iseg       !! segment number counter
    integer :: istat      !! status code from solver

    ! note: we don't know the problem size
    ! until we read the config file.

    ! first we have to read the config file:
    ! this will populate the global variables
    call read_config_file(config_file_name)

    call define_problem_size(n,m)

    ! initialize the solver:
    call me%initialize(     n                = n,            &
                            m                = m,            &
                            max_iter         = 100,          & ! maximum number of iteration
                            func             = nrho_func,    &
                            grad             = nrho_grad,    &
                            tol              = 1.0e-6_wp,    & ! tolerance
                            step_mode        = 4,            & ! 3-point "line search" (2 intervals)
                            n_intervals      = 2,            & ! number of intervals for step_mode=4
                            use_broyden      = .false.,      & ! broyden update
                            export_iteration = nrho_export   )

    call me%status(istat=istat)
    status_ok = istat == 0

    ! initialize the mission
    call me%mission%init(x)

    if (.not. status_ok) error stop 'error in initialize_the_solver'

    end subroutine initialize_the_solver
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a positive number the same magnitude as the input,
!  with only one significant digit.
!
!  If `mina` is present, then `max(mina,mag(a))` is returned
!
!  Examples:
!```
!     mag(1234.56)  -> 1000.0
!     mag(-999.99)  -> 900.0
!     mag(1.456e-4) -> 0.0001
!```

    pure elemental function mag(a,mina) result(m)

    implicit none

    real(wp),intent(in) :: a
    real(wp),intent(in),optional :: mina
    real(wp) :: m

    real(wp) :: x,tmp

    x = abs(a)

    if (x==0.0_wp) then
        if (.not. present(mina)) then
            m = 1.0_wp
        else
            m = mina
        end if
    else
        tmp = 10.0_wp ** floor(log10(x))
        m = tmp * floor(x/tmp)
        if (present(mina)) m = max(mina,m)
    end if

    end function mag
!*****************************************************************************************

!*****************************************************************************************
!>
!  Get the arrays for the problem

    subroutine get_problem_arrays(me,x,f)

    implicit none

    class(mission_type),intent(in) :: me
    real(wp),dimension(:),intent(out),optional :: x !! opt var vector [scaled]
    real(wp),dimension(:),intent(out),optional :: f !! constraint violations

    integer :: i          !! counter for index in `x`
    integer :: irev       !! counter for number of revs
    integer :: iseg       !! segment number counter
    integer :: j          !! counter

    if (present(x)) then

        ! x = [t01,x01,t02,x02,...,t0n,x0n] - scaled

        i = 0
        iseg = 0
        do irev = 1, n_revs

            if (irev==1) then

                ! the first one has an extra opt point at the initial periapsis passage

                x(i+1)       = me%segs(iseg+1)%data%t0
                x(i+2:i+7)   = me%segs(iseg+1)%data%x0_rotating

                x(i+8)       = me%segs(iseg+2)%data%t0
                x(i+9:i+14)  = me%segs(iseg+2)%data%x0_rotating

                x(i+15)      = me%segs(iseg+4)%data%t0
                x(i+16:i+21) = me%segs(iseg+4)%data%x0_rotating

                x(i+22)      = me%segs(iseg+6)%data%t0
                x(i+23:i+28) = me%segs(iseg+6)%data%x0_rotating

                x(i+29)      = me%segs(iseg+8)%data%t0
                x(i+30:i+35) = me%segs(iseg+8)%data%x0_rotating

                ! for next rev:
                i = i + 35
                iseg = iseg + 8

            else

                x(i+1)       = me%segs(iseg+2)%data%t0
                x(i+2:i+7)   = me%segs(iseg+2)%data%x0_rotating

                x(i+8)       = me%segs(iseg+4)%data%t0
                x(i+9:i+14)  = me%segs(iseg+4)%data%x0_rotating

                x(i+15)      = me%segs(iseg+6)%data%t0
                x(i+16:i+21) = me%segs(iseg+6)%data%x0_rotating

                x(i+22)      = me%segs(iseg+8)%data%t0
                x(i+23:i+28) = me%segs(iseg+8)%data%x0_rotating

                ! for next rev:
                i = i + 28
                iseg = iseg + 8

            end if

        end do

        ! scale the x vector:
        x = x / me%xscale

    end if

    if (present(f)) then

        ! f = [xf1-xf2, xf3-xf4, xf5-xf6, xf7-xf8, ... ]

        i = 0
        iseg = 0
        do irev = 1, n_revs

            f(i+1:i+6)   = -(me%segs(iseg+1)%data%xf_rotating - me%segs(iseg+2)%data%xf_rotating)
            f(i+7:i+12)  = -(me%segs(iseg+3)%data%xf_rotating - me%segs(iseg+4)%data%xf_rotating)
            f(i+13:i+18) = -(me%segs(iseg+5)%data%xf_rotating - me%segs(iseg+6)%data%xf_rotating)
            f(i+19:i+24) = -(me%segs(iseg+7)%data%xf_rotating - me%segs(iseg+8)%data%xf_rotating)

            i = i + 24
            iseg = iseg + 8

        end do

        !scale the f vector:
        f = f / me%fscale

    end if

    end subroutine get_problem_arrays
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take the `x` vector from the solver, and use it to populate all the segments

    subroutine put_x_in_segments(me,x_scaled)

    implicit none

    class(mission_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x_scaled !! opt var vector (scaled)

    integer :: i,j,irev,iseg  !! counter
    real(wp),dimension(8) :: t0 !! [days]
    real(wp),dimension(8) :: tf !! [days]
    real(wp),dimension(6,8) :: x0_rotating  !! rotating frame
    real(wp),dimension(size(x_scaled)) :: x !! opt var vector (unscaled)

    ! unscale the x vector:
    x = x_scaled * me%xscale

    ! first extract data from the opt var vector and put it into the segments:
    i = 0
    iseg = 0
    do irev = 1, n_revs

        ! extract the values:

        if (irev==1) then

            ! the first one has an extra opt point at the initial periapsis passage

            t0(1)               = x(i+1)
            x0_rotating(:,1)    = x(i+2:i+7)
            tf(1)               = t0(1) + period8

            t0(2)               = x(i+8)
            x0_rotating(:,2)    = x(i+9:i+14)
            tf(2)               = tf(1)

            t0(4)               = x(i+15)
            x0_rotating(:,4)    = x(i+16:i+21)
            tf(4)               = t0(4) - period8

            t0(3)               = t0(2)
            x0_rotating(:,3)    = x0_rotating(:,2)
            tf(3)               = tf(4)

            t0(5)               = t0(4)
            x0_rotating(:,5)    = x0_rotating(:,4)
            tf(5)               = t0(5) + period8

            t0(6)               = x(i+22)
            x0_rotating(:,6)    = x(i+23:i+28)
            tf(6)               = tf(5)

            t0(8)               = x(i+29)
            x0_rotating(:,8)    = x(i+30:i+35)
            tf(8)               = t0(8) - period8

            t0(7)               = t0(6)
            x0_rotating(:,7)    = x0_rotating(:,6)
            tf(7)               = tf(8)

            ! put the data for this rev into the segments:
            do j = 1, 8
                call me%segs(iseg+j)%set_input(t0(j),tf(j),x0_rotating(:,j))
            end do

            ! for next rev:
            i = i + 35
            iseg = iseg + 8

        else

            t0(1)               =  t0(8)                ! inherits from the previous rev
            x0_rotating(:,1)    =  x0_rotating(:,8)
            tf(1)               =  t0(1) + period8

            t0(2)               =  x(i+1)
            x0_rotating(:,2)    =  x(i+2:i+7)
            tf(2)               =  tf(1)

            t0(4)               =  x(i+8)
            x0_rotating(:,4)    =  x(i+9:i+14)
            tf(4)               =  t0(4) - period8

            t0(3)               =  t0(2)
            x0_rotating(:,3)    =  x0_rotating(:,2)
            tf(3)               =  tf(4)

            t0(5)               =  t0(4)
            x0_rotating(:,5)    =  x0_rotating(:,4)
            tf(5)               =  t0(5) + period8

            t0(6)               =  x(i+15)
            x0_rotating(:,6)    =  x(i+16:i+21)
            tf(6)               =  tf(5)

            t0(8)               =  x(i+22)
            x0_rotating(:,8)    =  x(i+23:i+28)
            tf(8)               =  t0(8) - period8

            t0(7)               =  t0(6)
            x0_rotating(:,7)    =  x0_rotating(:,6)
            tf(7)               =  tf(8)

            ! put the data for this rev into the segments:
            do j = 1, 8
                call me%segs(iseg+j)%set_input(t0(j),tf(j),x0_rotating(:,j))
            end do

            ! for next rev:
            i = i + 28
            iseg = iseg + 8

        end if

    end do

    end subroutine put_x_in_segments
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the sparsity pattern for the "forward-backward" NRHO problem.

    subroutine get_sparsity_pattern(me,irow,icol,&
                                    linear_irow,linear_icol,linear_vals,&
                                    maxgrp,ngrp)

    implicit none

    class(mission_type),intent(inout) :: me
    integer,dimension(:),allocatable,intent(out)  :: irow  !! sparsity pattern nonzero elements row indices
    integer,dimension(:),allocatable,intent(out)  :: icol  !! sparsity pattern nonzero elements column indices
    integer,dimension(:),allocatable,intent(out),optional  :: linear_irow !! linear sparsity pattern
                                                                          !! nonzero elements row indices
    integer,dimension(:),allocatable,intent(out),optional  :: linear_icol !! linear sparsity pattern nonzero
                                                                          !! elements column indices
    real(wp),dimension(:),allocatable,intent(out),optional :: linear_vals !! linear sparsity values (constant
                                                                          !! elements of the Jacobian)
    integer,intent(out),optional                           :: maxgrp      !! DSM sparsity partition
    integer,dimension(:),allocatable,intent(out),optional  :: ngrp        !! DSM sparsity partition

    ! class(mission_type),intent(inout) :: me
    ! integer,dimension(:),allocatable,intent(out) :: irow
    ! integer,dimension(:),allocatable,intent(out) :: icol

    integer :: i,j,k,ii,jj,icol_start,irow_start
    integer :: n_nonzero    !! number of nonzero elements in the jacobian
    integer :: iblock

    call define_problem_size(n_nonzero=n_nonzero)

    allocate(irow(n_nonzero))
    allocate(icol(n_nonzero))

    k = 0
    do iblock = 1, n_revs*4 ! block loop

        icol_start = (iblock-1)*7 + 1       ! [1-14, 8-21, 15-29, ...]
        irow_start = (iblock-1)*6 + 1       ! [1-6,  7-12, 13-18, ...]

        do ii = icol_start,icol_start+13
            do jj = irow_start,irow_start+5
                k = k + 1
                irow(k) = jj
                icol(k) = ii
            end do
        end do

    end do

    end subroutine get_sparsity_pattern
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the mission.

    subroutine initialize_the_mission(me,x)

    implicit none

    class(mission_type),intent(inout) :: me
    real(wp),dimension(:),allocatable,intent(out) :: x !! initial guess

    logical :: status_ok                        !! status flag
    integer :: i,j                              !! counter
    integer :: n_segs                           !! number of segments
    integer :: iseg                             !! segment counter
    real(wp) :: t_periapsis                     !! time of last periapsis crossing (in days from ref epoch)
    integer :: irev                             !! rev counter
    real(wp),dimension(8) :: t0                 !! [days]
    real(wp),dimension(8) :: tf                 !! [days]
    real(wp),dimension(6,8) :: x0_rotating      !! rotating frame
    integer :: n                                !! number of opt vars
    integer :: m                                !! number of constraints
    real(wp),dimension(:),allocatable :: dpert  !! perturbation step size array
    real(wp),dimension(:),allocatable :: xlow   !! lower bounds array
    real(wp),dimension(:),allocatable :: xhigh  !! upper bounds
    character(len=10) :: seg_name               !! segment name
    integer,dimension(:),allocatable :: irow    !! sparsity pattern nonzero elements row indices
    integer,dimension(:),allocatable :: icol    !! sparsity pattern nonzero elements column indices

    write(*,*) 'initialize_the_mission'

    ! problem size:
    call define_problem_size(n=n,m=m,n_segs=n_segs)

    write(*,*) ''
    write(*,*) 'n:                 ',n
    write(*,*) 'm:                 ',m
    write(*,*) 'number of segments:',n_segs
    write(*,*) 'number of revs:    ',n_revs
    write(*,*) ''

    ! arrays for the gradient computations:
    allocate(dpert(n))
    allocate(xlow(n))
    allocate(xhigh(n))
    dpert = 1.0e-06_wp

    ! need to construct some bounds so the sparsity can be computed.
    ! [note: these are not used since we are inputing the sparsity below]
    xlow  = -10.0_wp
    xhigh = 10.0_wp

    ! first set up the gradient computation in the base class:
    call me%initialize(n,m, &
                        xlow                       = xlow,&
                        xhigh                      = xhigh,&
                        dpert                      = dpert,&
                        problem_func               = my_func,&
                        info                       = nrho_grad_info,&
                        sparsity_mode              = 3,&        ! specified below
                        jacobian_method            = 3,&        ! standard central diff
                        perturb_mode               = 1,&        ! absolute mode
                        partition_sparsity_pattern = .true.,&   ! partition the pattern
                        cache_size                 = 1000 )

    ! generate and set the sparsity pattern for this problem:
    call me%get_sparsity_pattern(irow,icol)
    call me%set_sparsity_pattern(irow,icol)

    ! pointers:
    allocate(me%eph)
    allocate(me%grav)

    ! set up the ephemeris:
    write(*,*) 'loading ephemeris file: '//trim(ephemeris_file)
    call me%eph%initialize(filename=ephemeris_file,status_ok=status_ok)
    if (.not. status_ok) error stop 'error initializing ephemeris'

    ! set up the force model [main body is moon]:
    call me%grav%initialize(gravfile,grav_n,grav_m,status_ok)
    if (.not. status_ok) error stop 'error initializing gravity model'

    ! now, we set up the segment structure for the problem we are solving:
    ! This is for the "forward-backward" method from the paper (see Figure 2b):
    allocate(me%segs(n_segs))

    ! set up the integrators:
    do i = 1, n_segs

        call me%segs(i)%initialize(n_eoms,maxnum,ballistic_derivs,&
                                   [rtol],[atol],&
                                   report=trajectory_export_func)

        ! for now, each seg points to the global ones for the whole mission
        ! .... reconsider this ....
        me%segs(i)%grav => me%grav
        me%segs(i)%eph  => me%eph

        ! name all the segments:
        write(seg_name,'(I10)') i
        me%segs(i)%name = trim(adjustl(seg_name))

    end do

    ! load the initial guess from the patch point file:

    i = 0
    iseg = 0
    t_periapsis = 0.0_wp
    do irev = 1, n_revs

        ! these are all unscaled values:

        t0(1) = t_periapsis + periapsis%t
        x0_rotating(:,1) = periapsis%rv
        tf(1) = t0(1) + period8

        t0(2) = t_periapsis + quarter%t
        x0_rotating(:,2) = quarter%rv
        tf(2) = tf(1)

        t0(4) = t_periapsis + apoapsis%t
        x0_rotating(:,4) = apoapsis%rv
        tf(4) = t0(4) - period8

        t0(3) = t0(2)
        x0_rotating(:,3) = x0_rotating(:,2)
        tf(3) = tf(4)

        t0(5) = t0(4)
        x0_rotating(:,5) = x0_rotating(:,4)
        tf(5) = t0(5) + period8

        t0(6) = t_periapsis + apoapsis%t + quarter%t
        x0_rotating(:,6) = [quarter%rv(1),-quarter%rv(2),quarter%rv(3),-quarter%rv(4),quarter%rv(5),-quarter%rv(6)]
        tf(6) = tf(5)

        t0(8) = t_periapsis + period
        x0_rotating(:,8) = periapsis%rv
        tf(8) = t0(8) - period8

        t0(7) = t0(6)
        x0_rotating(:,7) = x0_rotating(:,6)
        tf(7) = tf(8)

        do j = 1, 8

            ! put the data for this rev into the segments:
            call me%segs(iseg+j)%set_input(t0(j),tf(j),x0_rotating(:,j))

            ! also set the scales:
            ! the t0 gradually grows larger so use the initial value for each segment.
            ! for the states, just use the constant values

            me%segs(iseg+j)%data%t0_scale = abs(mag(2.0_wp*me%segs(iseg+j)%data%t0))    ! magic number

            ! do the same for the states, but just in case, specify the min values:
            !me%segs(iseg+j)%data%x0_rotating_scale = xscale_x0   ! original
            me%segs(iseg+j)%data%x0_rotating_scale = abs(mag(2.0_wp*me%segs(iseg+j)%data%x0_rotating, xscale_x0)) ! new

            me%segs(iseg+j)%data%xf_rotating_scale = fscale_xf  ! these are all the same for the constraints

        end do

        ! for next rev:
        i = i + 35
        iseg = iseg + 8
        t_periapsis = t_periapsis + period ! add another period

    end do

    allocate(x(n))    ! size the opt var vector
    allocate(me%xscale(n))
    allocate(me%xname(n))
    allocate(me%fscale(m))

    ! also set the scale factors:
    call me%get_scales_from_segs()

    ! get the initial guess from the mission:
    call me%get_problem_arrays(x=x)

    write(*,*) 'Done with initialize_the_mission.'

    end subroutine initialize_the_mission
!*****************************************************************************************

!*****************************************************************************************
!>
!  Populate the `xscale` and `fscale` problem arrays from the segment data.

    subroutine get_scales_from_segs(me)

    implicit none

    class(mission_type),intent(inout) :: me

    integer :: i,j    !! counter
    integer :: iseg   !! segment number counter
    integer :: n_segs !! number of segments
    character(len=10) :: iseg_str  !! segment number string

    character(len=*),parameter :: t0_label = 'T0 (day)'
    character(len=9),dimension(6),parameter :: x0_label = ['Rx (km)  ',&
                                                           'Ry (km)  ',&
                                                           'Rz (km)  ',&
                                                           'Vx (km/s)',&
                                                           'Vy (km/s)',&
                                                           'Vz (km/s)']

    ! iseg loop: x: [1,2,4,...n_segs]
    !            f: [  2,4,...n_segs]

    call define_problem_size(n_segs=n_segs)

    ! x scales - segment 1:
    me%xscale(1)   = me%segs(1)%data%t0_scale
    me%xscale(2:7) = me%segs(1)%data%x0_rotating_scale
    me%xname(1) = 'SEG1 '//t0_label
    do j=2,7
        me%xname(j) = 'SEG1 '//x0_label(j-1)
    end do
    i = 8

    ! x scales - the rest:
    do iseg = 2, n_segs, 2
        me%xscale(i)       = me%segs(iseg)%data%t0_scale
        me%xscale(i+1:i+6) = me%segs(iseg)%data%x0_rotating_scale

        write(iseg_str,'(I10)') iseg
        me%xname(i) = 'SEG'//trim(adjustl(iseg_str))//' '//t0_label
        do j=1,6
            me%xname(i+j) = 'SEG'//trim(adjustl(iseg_str))//' '//x0_label(j)
        end do

        i = i + 7
    end do

    ! f scales:
    i = 1
    do iseg = 2, n_segs, 2
        me%fscale(i:i+5)  = me%segs(iseg)%data%xf_rotating_scale
        i = i + 6
    end do

    !write(*,*) ''
    !write(*,*) 'xscale: ', me%xscale
    !write(*,*) ''
    !write(*,*) 'fscale: ', me%fscale
    !write(*,*) ''

    end subroutine get_scales_from_segs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Equations of motion for a ballistic segment.
!
!@note In this example, all the segments have the same equations of motion!

    subroutine ballistic_derivs(me,t,x,xdot)

    implicit none

    class(ddeabm_class),intent(inout) :: me
    real(wp),intent(in)               :: t    !! time [sec from epoch]
    real(wp),dimension(:),intent(in)  :: x    !! state [r,v] in inertial frame (moon-centered)
    real(wp),dimension(:),intent(out) :: xdot !! derivative of state (\( dx/dt \))

    real(wp),dimension(3) :: r,rb,v
    reaL(wp),dimension(6) :: rv_earth_wrt_moon,rv_sun_wrt_moon
    real(wp),dimension(3,3) :: rotmat
    real(wp),dimension(3) :: a_geopot
    real(wp),dimension(3) :: a_earth
    real(wp),dimension(3) :: a_sun
    real(wp),dimension(3) :: a_third_body
    real(wp) :: et !! ephemeris time of `t`
    logical :: status_ok

    select type (me)

    class is (segment)

        ! get state:
        r = x(1:3)
        v = x(4:6)

        ! compute ephemeris time [sec]:
        et = et_ref + t

        ! geopotential gravity:
        rotmat = icrf_to_iau_moon(et)   ! rotation matrix from inertial to body-fixed Moon frame
        rb = matmul(rotmat,r)           ! r in body-fixed frame
        call me%grav%get_acc(rb,grav_n,grav_m,a_geopot)  ! get the acc due to the geopotential
        a_geopot = matmul(transpose(rotmat),a_geopot)    ! convert acc back to inertial frame

        ! third-body state vectors (wrt the central body, which is the moon in this case):
        ! [inertial frame]
        call me%eph%get_rv(et,body_earth,body_moon,rv_earth_wrt_moon,status_ok)
        call me%eph%get_rv(et,body_sun,body_moon,rv_sun_wrt_moon,status_ok)

        ! third-body perturbation (earth & sun):
        a_third_body = 0.0_wp
        call third_body_gravity(r,rv_earth_wrt_moon(1:3),mu_earth,a_earth)
        call third_body_gravity(r,rv_sun_wrt_moon(1:3),mu_sun,a_sun)
        a_third_body = a_earth + a_sun

        !total derivative vector:
        xdot(1:3) = v
        xdot(4:6) = a_geopot + a_third_body

    class default

        error stop 'invalid class in ballistic_derivs'

    end select

    end subroutine ballistic_derivs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sets all the info in a segment for it to be propagated.

    subroutine set_segment_inputs(me,t0,tf,x0_rotating)

    implicit none

    class(segment),intent(inout) :: me
    real(wp),intent(in) :: t0 !! [days]
    real(wp),intent(in) :: tf !! [days]
    real(wp),dimension(6),intent(in) :: x0_rotating  !! state in rotating frame

    real(wp),dimension(6) :: x0_inertial !! inertial frame
    type(icrf_frame) :: inertial
    type(two_body_rotating_frame) :: rotating
    real(wp) :: et0
    logical :: status_ok

    ! note that the position vectors are in the rotating frame,
    ! have to transform to inertial to integrate.
    ! create the frames for the transformation:
    et0 = et_ref + t0 * day2sec  ! convert to ephemeris time [sec]
    rotating = two_body_rotating_frame(primary_body=body_earth,&
                                        secondary_body=body_moon,&
                                        center=center_at_secondary_body,&
                                        et=et0)
    inertial = icrf_frame(b=body_moon)

    ! from rotating to inertial:
    call rotating%transform(x0_rotating,inertial,et0,me%eph,x0_inertial,status_ok)
    if (.not. status_ok) error stop 'transformation error in set_segment_inputs'

    ! write(*,*) ''
    ! write(*,*) 'rotating to inertial'
    ! write(*,'(A/,*(E30.16/))') 'et0:        ', et0
    ! write(*,'(A/,*(E30.16/))') 'x0_rotating:', x0_rotating
    ! write(*,'(A/,*(E30.16/))') 'x0_inertial:', x0_inertial
    ! write(*,*) ''
    !PAUSE

    ! set the inputs (needed by the propagator)
    me%data%t0          = t0
    me%data%x0_rotating = x0_rotating
    me%data%x0          = x0_inertial
    me%data%tf          = tf

    end subroutine set_segment_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  After propagating a segment, this gets the outputs.

    subroutine get_segment_outputs(me,xf,xf_rotating)

    implicit none

    class(segment),intent(in) :: me
    real(wp),dimension(6),intent(out) :: xf             !!  inertial frame
    real(wp),dimension(6),intent(out) :: xf_rotating    !!  rotating frame

    xf = me%data%xf
    xf_rotating = me%data%xf_rotating

    end subroutine get_segment_outputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the outputs of a segment, assuming it has been propagated elsewhere (e.g., in a coarray?)

    subroutine set_segment_outputs(me,xf,xf_rotating)

    implicit none

    class(segment),intent(inout) :: me
    real(wp),dimension(6),intent(in) :: xf  !! inertial frame
    real(wp),dimension(6),intent(in) :: xf_rotating  !! inertial frame

    me%data%xf = xf
    me%data%xf_rotating = xf_rotating

    end subroutine set_segment_outputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Propagate a segment (assumes the inputs have already been populated)

    subroutine propagate_segment(me,mode)

    implicit none

    class(segment),intent(inout) :: me
    integer,intent(in),optional :: mode  !! 1 - don't report steps, 2 - report steps (for plotting)

    integer  :: idid
    real(wp) :: t
    real(wp) :: tf
    real(wp),dimension(6) :: x
    type(icrf_frame) :: inertial
    type(two_body_rotating_frame) :: rotating
    real(wp) :: etf
    real(wp),dimension(6) :: xf !! inertial - from the propagation
    real(wp),dimension(6) :: xf_rotating
    logical :: status_ok
    integer :: integration_mode

    if (present(mode)) then
        integration_mode = mode
    else
        integration_mode = 1
    end if

    t  = me%data%t0 * day2sec  ! initial time in seconds from epoch
    tf = me%data%tf * day2sec  ! final time in seconds from epoch
    x  = me%data%x0   ! inertial state

    !write(*,*) ''
    !write(*,*) 'propagate segment '//me%name, t0, tf
    !write(*,*) 't0:',t0
    !write(*,'(A,*(F15.3,1X))') 'x0:',x
    !write(*,*) 'tf:',tf

    call me%first_call()  !restarting the integration
    call me%integrate(t,x,tf,idid=idid,integration_mode=integration_mode)
    if (idid<0) then
        write(*,'(A,*(I5/))')    'idid: ',idid
        error stop 'error in integrator'
    end if

    xf = x  ! final state [inertial frame]

    ! also save the rotating frame state at tf (to compute the constraint violations):
    etf = et_ref + tf  ! convert to ephemeris time [sec]
    inertial = icrf_frame(b=body_moon)
    rotating = two_body_rotating_frame(primary_body=body_earth,&
                                       secondary_body=body_moon,&
                                       center=center_at_secondary_body,&
                                       et=etf)

    ! from inertial to rotating:
    call inertial%transform(x,rotating,etf,me%eph,xf_rotating,status_ok)
    if (.not. status_ok) error stop 'transformation error in propagate_segment'

    ! put final states in the segment:
    call me%set_outputs(xf,xf_rotating)

    end subroutine propagate_segment
!*****************************************************************************************

!*****************************************************************************************
!>
!  Given the functions to be computed, returns the segments that need to be propagated.

    subroutine segs_to_propagate(me,funcs_to_compute,isegs_to_propagate)

    implicit none

    class(mission_type),intent(inout) :: me
    integer,dimension(:),intent(in),optional :: funcs_to_compute !! the indices of f to compute
    integer,dimension(:),allocatable,intent(out) :: isegs_to_propagate  !! segment numbers

    integer :: i,j  !! counters

    ! In this example, there are no dependencies, so we assume they
    ! can be propagated in any order.

    if (present(funcs_to_compute)) then

        ! this is for the "forward-backward" problem formulation:

        !       1-6      7-12    13-19     20-26
        ! f = [xf1-xf2, xf3-xf4, xf5-xf6, xf7-xf8, ... ]
        !
        !      0     1     2
        ! f = [123456123456123456]
        !
        ! 0: 1-2
        ! 1: 3-4
        ! 2: 5-6
        ! 3: 7-8

        do i = 1, size(funcs_to_compute)
            j = (funcs_to_compute(i)-1) / 6 ! 0,1,2,...
            call add_it(j*2+1,isegs_to_propagate) ! ... is this correct ?
            call add_it(j*2+2,isegs_to_propagate)
        end do

    else
        ! propagate all the segments:
        isegs_to_propagate = [(i, i=1,size(me%segs))]
    end if

    !write(*,*) ''
    !write(*,*) 'funcs_to_compute:',funcs_to_compute
    !write(*,*) ''
    !write(*,*) 'isegs_to_propagate:',isegs_to_propagate
    !write(*,*) ''

    contains

        !********************************
        subroutine add_it(segnum,array)

        implicit none

        integer,intent(in) :: segnum
        integer,dimension(:),allocatable,intent(inout) :: array

        if (allocated(array)) then
            if (.not. any(array==segnum)) then
                array = [array,segnum]
            end if
        else
            array = [segnum] ! auto LHS allocation
        end if

        end subroutine add_it
        !********************************

    end subroutine segs_to_propagate
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the constraint violation vector for the mission.

    subroutine constraint_violations(me,x,f,funcs_to_compute)

    implicit none

    class(mission_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x     !! opt var vector for the mission [n]
    real(wp),dimension(:),intent(out) :: f     !! constraint violation vector for the mission [m]
    integer,dimension(:),intent(in),optional :: funcs_to_compute !! the indices of f to compute

    integer :: i,j  !! counters
    integer,dimension(:),allocatable :: isegs

    if (.not. allocated(me%segs)) error stop 'error: segs is not allocated'

    ! first extract data from the opt var vector
    ! and put it into all the segments:
    call me%put_x_in_segments(x)

    ! get the list of segments that needs to be propagated:
    ! this depends on the functions that are being computed,
    ! and the sparsity pattern
    call me%segs_to_propagate(funcs_to_compute,isegs)

    !====================================
    ! now propagate the segments:
    do i = 1, size(isegs)

        call me%segs(isegs(i))%propagate()

        ! ... is this sort of how the coarray part would work ??
        !
        ! call co_segs[...]%data = me%segs(isegs(i))%data   ! copy data to coarray [or use put_data_in_segment]
        ! call co_segs[...]%propagate()                     ! propagate in the coarray
        ! call co_segs[...]%get_segment_outputs(xf)         ! get the results form coarray
        ! call me%segs(isegs(i))%set_segment_outputs(xf) ! copy results back to real segment

        ! ... but, this is deep inside the solver, we don't want all the images
        !     to be running the solver, we just want all the images to be
        !     available so we can propagate the segments ...

    end do
    !====================================

    ! now that all the segments are propagated, we
    ! compute the function (constraint violations for the
    ! final states of segment segment)
    call me%get_problem_arrays(f=f)

    end subroutine constraint_violations
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the solver function (all the constraint violations)

    subroutine nrho_func(me,x,f)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:),intent(out)    :: f

    integer :: i

    select type (me)
    class is (my_solver_type)
        call me%mission%constraint_violations(x,f)
    class default
        error stop 'invalid class in nrho_func'
    end select

    end subroutine nrho_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the gradient of the solver function (Jacobian matrix)

    subroutine nrho_grad(me,x,g)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:,:),intent(out)  :: g

    real(wp),dimension(:,:),allocatable :: jac !! the jacobian matrix returned by `numdiff`
        ! ...note: need to modify so it doesn't
        !          have to return an allocatable array

    integer :: i  !! seg number counter

    select type (me)
    class is (my_solver_type)

        ! first let's cache all the segment data:
        do i=1,size(me%mission%segs)
            call me%mission%segs(i)%cache()
        end do

        ! use numdiff to compute the jacobian matrix (dense version)
        call me%mission%compute_jacobian_dense(x,jac)
        g = jac

        ! restore data just in case:
        do i=1,size(me%mission%segs)
            call me%mission%segs(i)%uncache()
        end do

    class default
        error stop 'invalid class in nrho_grad'
    end select

    end subroutine nrho_grad
!*****************************************************************************************

!*****************************************************************************************
!>
!  Info function for the gradients.
!
!  Used to restore the segments back to their unperturbed
!  state when perturbing the opt vars.
!
!@note This may not be strictly necessary for this particular problem...
!      ... need to think about this some more ...

    subroutine nrho_grad_info(me,column,i,x)

    implicit none

    class(numdiff_type),intent(inout)  :: me
    integer,dimension(:),intent(in) :: column !! the columns being computed.
    integer,intent(in) :: i                   !! perturbing these columns for the `i`th time (1,2,...)
    real(wp),dimension(:),intent(in)  :: x    !! the nominal variable vector

    integer :: iseg !! segment number counter

    select type (me)
    class is (mission_type)

        ! restore data to all the segments
        ! from the unperturbed mission:
        ! [just in case]
        do iseg=1,size(me%segs)
            call me%segs(iseg)%uncache()
        end do

    class default
        error stop 'invalid class in nrho_grad_info'
    end select

    end subroutine nrho_grad_info
!*****************************************************************************************

!*****************************************************************************************
!>
!  Wrapper for the constraint violation function (use by NumDiff)

    subroutine my_func(me,x,f,funcs_to_compute)

    implicit none

    class(numdiff_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)  :: x
    real(wp),dimension(:),intent(out) :: f
    integer,dimension(:),intent(in)   :: funcs_to_compute !! the function indices to
                                                          !! compute in the full `f` vector

    select type (me)
    class is (mission_type)
        call me%constraint_violations(x,f,funcs_to_compute)
    class default
        error stop 'invalid class in nrho_grad'
    end select

    end subroutine my_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  export an iteration from the solver

    subroutine nrho_export(me,x,f,iter)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:),intent(in)     :: f
    integer,intent(in)                   :: iter !! iteration number

    if (iter==1) then
        write(*,'(A4,1X,*(A30,1X))') 'ITER', 'NORM(X)', 'NORM(F)'
    end if

    write(*,'(I4,1X,*(F30.16,1X))') iter, norm2(x), norm2(f)
    !write(*,'(I4,1X,*(F30.16,1X))') iter, x

    end subroutine nrho_export
!*****************************************************************************************

!*****************************************************************************************
!>
!  Plot the trajectory using matplotlib.
!
!@note It is assumed that all the data is present in the segments needed to propagate.

    subroutine plot_trajectory(me,filename,export_trajectory)

    implicit none

    class(mission_type),intent(inout) :: me
    character(len=*),intent(in) :: filename !! plot file name [without extension]
    logical,intent(in),optional :: export_trajectory    !! if true, a text trajectory
                                                        !! file is also produced
                                                        !! (that can be read by MKSPK)

    type(pyplot) :: plt  !! for generating the plots
    integer :: iseg  !! segment number counter
    integer :: istat !! pyplot status code
    character(len=10) :: iseg_str !! string version of segment number
    logical :: export  !! if the txt file is to be produced as well.
    integer :: iunit !! file unit for the txt file
    integer :: i !! counter for txt file write
    integer :: istart,iend,istep,iendprev  !! index counters

    logical :: plot_rotating = .true.    !! if true, the rotating state is plotted.
                                         !! if false, the inertial state is plotted.

    ! optional arguments:
    if (present(export_trajectory)) then
        export = export_trajectory
    else
        export = .false.
    end if

    ! initialize the plot:
    call plt%initialize(grid=.true.,&
                        xlabel='\n\n\nx (km)',&
                        ylabel='\n\n\ny (km)',&
                        zlabel='\n\n\nz (km)',&
                        figsize=[20,20],&
                        font_size       = 20,&
                        axes_labelsize  = 25,&
                        xtick_labelsize = 20,&
                        ytick_labelsize = 20,&
                        ztick_labelsize = 20,&
                        title='NRHO Trajectory: '//trim(filename),&
                        legend=.false.,&
                        axis_equal=.true.,&
                        mplot3d=.true.)

    if (export) then
        open(newunit=iunit,file=trim(filename)//'_'//get_case_name()//'.txt',status='REPLACE',iostat=istat)
        if (istat/=0) error stop 'error opening trajectory file.'
    end if

    do iseg = 1, size(me%segs)
        ! destroy all trajectories:
        call me%segs(iseg)%traj_inertial%destroy()
        call me%segs(iseg)%traj_rotating%destroy()
    end do

    do iseg = 1, size(me%segs)

        write(iseg_str,'(I10)') iseg

        ! generate the trajectory for this segment:
        call me%segs(iseg)%propagate(mode=2)    ! [export points]

        ! export to plot:
        if (plot_rotating) then
            call plt%add_3d_plot(x=me%segs(iseg)%traj_rotating%x,&
                                 y=me%segs(iseg)%traj_rotating%y,&
                                 z=me%segs(iseg)%traj_rotating%z,&
                                 label='seg'//trim(adjustl(iseg_str)),&
                                 linestyle='r-',linewidth=1,istat=istat)
        else
            call plt%add_3d_plot(x=me%segs(iseg)%traj_inertial%x,&
                                 y=me%segs(iseg)%traj_inertial%y,&
                                 z=me%segs(iseg)%traj_inertial%z,&
                                 label='seg'//trim(adjustl(iseg_str)),&
                                 linestyle='r-',linewidth=1,istat=istat)
        end if

        if (export) then

            ! write the trajectory data:     ! Note: THIS HAS TO BE THE INERTIAL FRAME

            if (me%segs(iseg)%traj_inertial%et(2) > me%segs(iseg)%traj_inertial%et(1)) then
                ! forward propagated
                istart = 1
                iend   = size(me%segs(iseg)%traj_inertial%et)
                istep  = 1
            else
                ! backward propagated
                istart = size(me%segs(iseg)%traj_inertial%et)
                iend   = 1
                istep  = -1
            end if

            if (iseg>1) then
                ! if the last one written is the same as the first one here, don't write it
                ! [MKSPK does not allow duplicate time tags]
                if (me%segs(iseg-1)%traj_inertial%et(iendprev) >= me%segs(iseg)%traj_inertial%et(istart)) then
                    !write(*,*) 'duplicate point: ', me%segs(iseg-1)%traj_inertial%et(iendprev)
                    istart = istart + istep
                end if
            end if

            do i=istart,iend,istep
                ! inertial trajectory for SPK:
                write(iunit,'(*(E30.16E3,A,1X))',iostat=istat) &
                                               me%segs(iseg)%traj_inertial%et(i),';',&
                                               me%segs(iseg)%traj_inertial%x(i) ,';',&
                                               me%segs(iseg)%traj_inertial%y(i) ,';',&
                                               me%segs(iseg)%traj_inertial%z(i) ,';',&
                                               me%segs(iseg)%traj_inertial%vx(i),';',&
                                               me%segs(iseg)%traj_inertial%vy(i),';',&
                                               me%segs(iseg)%traj_inertial%vz(i),''
            end do
            iendprev = iend
        end if

    end do

    ! add the moon as a sphere:
    call plt%add_sphere(r=r_moon,xc=zero,yc=zero,zc=zero,istat=istat,color='Grey')

    ! save the plot:
    ! call plt%showfig(pyfile=trim(filename)//'.py',istat=istat)
    call plt%savefig(trim(filename)//'_'//get_case_name()//'.png',&
                     pyfile=trim(filename)//'_'//get_case_name()//'.py',dpi='300',&
                     istat=istat)

    ! cleanup:
    call plt%destroy()
    if (export) then
        close(iunit,iostat=istat)
    end if

    end subroutine plot_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add a point to a trajectory variable.

    subroutine add_point_to_trajectory(me,et,x)

    implicit none

    class(trajectory),intent(inout) :: me !! the trajectory
    real(wp),intent(in) :: et
    real(wp),dimension(6),intent(in) :: x

    if (.not. allocated(me%et)) then
        allocate(me%et(1)); me%et(1) = et
        allocate(me%x(1));  me%x(1)  = x(1)
        allocate(me%y(1));  me%y(1)  = x(2)
        allocate(me%z(1));  me%z(1)  = x(3)
        allocate(me%vx(1)); me%vx(1) = x(4)
        allocate(me%vy(1)); me%vy(1) = x(5)
        allocate(me%vz(1)); me%vz(1) = x(6)
    else
        me%et = [me%et, et   ] ! auto lhs allocation
        me%x  = [me%x,  x(1) ]
        me%y  = [me%y,  x(2) ]
        me%z  = [me%z,  x(3) ]
        me%vx = [me%vx, x(4) ]
        me%vy = [me%vy, x(5) ]
        me%vz = [me%vz, x(6) ]
    end if

    end subroutine add_point_to_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Export a trajectory point during propagation.
!
!@note It saves both the inertial and rotating trajectory
!      for plotting and exporting later.

    subroutine trajectory_export_func(me,t,x)

    implicit none

    class(ddeabm_class),intent(inout) :: me
    real(wp),intent(in)               :: t    !! time from integrator (sec from sim epoch)
    real(wp),dimension(:),intent(in)  :: x    !! state [moon-centered inertial frame]

    real(wp) :: et !! ephemeris time [sec]
    real(wp),dimension(6) :: x_rotating !! the state to export
    type(icrf_frame) :: inertial
    type(two_body_rotating_frame) :: rotating
    logical :: status_ok !! transformation status flag

    select type (me)
    class is (segment)

        et = et_ref + t  ! convert to ephemeris time [sec]

        ! convert state to the moon-centered rotating frame
        inertial = icrf_frame(b=body_moon)
        rotating = two_body_rotating_frame(primary_body=body_earth,&
                                           secondary_body=body_moon,&
                                           center=center_at_secondary_body,&
                                           et=et)
        ! from inertial to rotating:
        call inertial%transform(x,rotating,et,me%eph,x_rotating,status_ok)
        if (.not. status_ok) error stop 'transformation error in propagate_segment'

        ! save both the inertial and rotating trajectories:
        call me%traj_inertial%add(et,x)
        call me%traj_rotating%add(et,x_rotating)

    class default
        error stop 'invalid class in trajectory_export_func'
    end select

    end subroutine trajectory_export_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the config file that defines the problem
!  to be solved and set all the global variables.
!
!@todo Add patch point interpolation based on \( r_p \)
!      if the specified value is not in the file.

    subroutine read_config_file(filename)

    implicit none

    character(len=*),intent(in) :: filename  !! the JSON config file to read

    type(json_file) :: json !! the config file structure
    type(csv_file) :: csv   !! the patch point file
    logical :: found
    logical :: status_ok
    real(wp),dimension(:),allocatable :: rp_vec
    real(wp),dimension(:),allocatable :: t_vec
    real(wp),dimension(:),allocatable :: x_vec
    real(wp),dimension(:),allocatable :: y_vec
    real(wp),dimension(:),allocatable :: z_vec
    real(wp),dimension(:),allocatable :: vx_vec
    real(wp),dimension(:),allocatable :: vy_vec
    real(wp),dimension(:),allocatable :: vz_vec
    integer :: i, irow
    logical :: tmp  !! for optional logical inputs

    call json%initialize()
    write(*,*) 'Loading config file: '//trim(filename)
    call json%load_file(filename=filename)
    if (json%failed()) error stop 'error loading config file'

    call json%get('rp',            rp,            found)
    if (.not. found) error stop 'rp variable not found in config file'
    call json%get('N_or_S',        N_or_S,        found)
    if (.not. found) error stop 'N_or_S variable not found in config file'
    call json%get('L1_or_L2',      L1_or_L2,      found)
    if (.not. found) error stop 'L1_or_L2 variable not found in config file'
    call json%get('year',          year,          found)
    if (.not. found) error stop 'year variable not found in config file'
    call json%get('month',         month,         found)
    if (.not. found) error stop 'month variable not found in config file'
    call json%get('day',           day,           found)
    if (.not. found) error stop 'day variable not found in config file'
    call json%get('hour',          hour,          found)
    if (.not. found) error stop 'hour variable not found in config file'
    call json%get('minute',        minute,        found)
    if (.not. found) error stop 'minute variable not found in config file'
    call json%get('sec',           sec,           found)
    if (.not. found) error stop 'sec variable not found in config file'
    call json%get('n_revs',        n_revs,        found)
    if (.not. found) error stop 'n_revs variable not found in config file'
    call json%get('ephemeris_file',ephemeris_file,found)
    if (.not. found) error stop 'ephemeris_file variable not found in config file'
    call json%get('gravfile',gravfile,found)
    if (.not. found) error stop 'gravfile variable not found in config file'
    call json%get('patch_point_file',patch_point_file,found)
    if (.not. found) error stop 'patch_point_file variable not found in config file'

    ! optional ones:
    ! [if not present, then the defaults are used]
    call json%get('generate_plots', tmp, found)
    if (found) generate_plots = tmp
    call json%get('generate_trajectory_files', tmp, found)
    if (found) generate_trajectory_files = tmp

    if (json%failed()) error stop 'error loading config file'
    call json%destroy() ! cleanup

    ! compute reference epoch from the date (TDB):
    et_ref = jd_to_et(julian_date(year,month,day,hour,minute,sec))

    ! read the patch point file and load it:
    ! [TODO we can also add the L1 file]
    if (L1_or_L2/='L2') error stop 'only L2 is currently supported'

    ! read the CSV file:
    call csv%read(patch_point_file,skip_rows=[1,2],status_ok=status_ok)
    if (.not. status_ok) error stop 'error reading csv file'

    ! get the data from the CSV file:
    call csv%get(1,rp_vec,status_ok)  ! rp values
    if (.not. status_ok) error stop 'error reading rp column from csv file'

    ! find the correct row with the specified rp value:
    irow = 0
    do i=1,size(rp_vec) ! skip the first two header rows
        if (rp_vec(i)==rp) then
            irow = i
            exit
        end if
    end do
    if (irow==0) error stop 'Rp value not found in file'

    ! populate all the patch point structures:
    call get_patch_point(periapsis,4)
    call get_patch_point(quarter, 20)
    call get_patch_point(apoapsis,36)

    ! compute some time variables:
    period  = apoapsis%t * 2.0_wp  ! NRHO period [days]
    period8 = period / 8.0_wp      ! 1/8 of NRHO period [days]

    ! cleanup:
    call csv%destroy()

    contains
!*****************************************************************************************

    !**********************************************
    !>
    !  Populate the patch point structure.

        subroutine get_patch_point(pp,istart)

        implicit none

        type(patch_point),intent(out) :: pp  !! the patch point data to populate
        integer,intent(in) :: istart !! column index where data starts

        call csv%get(istart  ,t_vec,status_ok)
        if (.not. status_ok) error stop 'error reading t_vec from csv file'
        call csv%get(istart+1,x_vec,status_ok)
        if (.not. status_ok) error stop 'error reading x_vec from csv file'
        call csv%get(istart+2,y_vec,status_ok)
        if (.not. status_ok) error stop 'error reading y_vec from csv file'
        call csv%get(istart+3,z_vec,status_ok)
        if (.not. status_ok) error stop 'error reading z_vec from csv file'
        call csv%get(istart+4,vx_vec,status_ok)
        if (.not. status_ok) error stop 'error reading vx_vec from csv file'
        call csv%get(istart+5,vy_vec,status_ok)
        if (.not. status_ok) error stop 'error reading vy_vec from csv file'
        call csv%get(istart+6,vz_vec,status_ok)
        if (.not. status_ok) error stop 'error reading vz_vec from csv file'

        if (N_or_S=='S') then
            ! the data in the file is for the South family
            pp = patch_point(t = t_vec(irow),&
                             rv = [ x_vec(irow),&
                                    y_vec(irow),&
                                    z_vec(irow),&
                                    vx_vec(irow),&
                                    vy_vec(irow),&
                                    vz_vec(irow)])
        elseif (N_or_S=='N') then
            pp = patch_point(t = t_vec(irow),&
                             rv = [ x_vec(irow),&
                                    y_vec(irow),&
                                    -z_vec(irow),&
                                    vx_vec(irow),&
                                    vy_vec(irow),&
                                    -vz_vec(irow)])
        else
            error stop 'invalid value for N_or_S'
        end if

        end subroutine get_patch_point
    !**********************************************

    end subroutine read_config_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a string that can be used for file names, etc.
!  (read the values of the global variables).
!
!@note We are only saving `rp` and `sec` as an integer.

    function get_case_name() result(case_name)

    implicit none

    character(len=:),allocatable :: case_name

    character(len=10) :: istr

    case_name = int_to_string(year,4)//int_to_string(month,2)//&
                int_to_string(day,2)//int_to_string(hour,2)//&
                int_to_string(minute,2)//int_to_string(int(sec),2)//&
                '_RP='//int_to_string(int(rp))//'_'//&
                L1_or_L2//'_'//N_or_S//'_NREVS='//int_to_string(n_revs)

    contains

    !*************************************************
        function int_to_string(i,ipad) result(str)
        !! integer to string with zero padding
        implicit none
        integer,intent(in) :: i
        character(len=:),allocatable :: str
        integer,intent(in),optional :: ipad !! number of digits to pad with zeros
        character(len=100) :: istr
        integer :: istat
        write(istr,'(I100)',iostat=istat) i
        if (istat==0) then
            str = trim(adjustl(istr))
            if (present(ipad)) then
                if (len(str)<ipad) then
                    str = repeat('0',ipad-len(str))//str
                end if
            end if
        else
            str ='****'
        end if
        end function int_to_string
    !*************************************************

    end function get_case_name
!*****************************************************************************************

!*****************************************************************************************
!>
!  Write the current `x` variables to a file (compatible with the Copernicus optvar file).

    subroutine write_optvars_to_file(me,filename,x)

    implicit none

    class(mission_type),intent(inout) :: me
    character(len=*),intent(in)       :: filename !! file name [without extension]
    real(wp),dimension(:),intent(in)  :: x        !! solver opt vars vector

    integer           :: i     !! counter
    type(json_file)   :: json  !! for writing the solution file
    character(len=10) :: istr  !! x index string

    call json%initialize()

    do i = 1, size(x)
        write(istr,'(I10)') i
        istr = adjustl(istr)
        call json%add('xvec('//trim(istr)//').label',trim(me%xname(i)))
        call json%add('xvec('//trim(istr)//').value',x(i)*me%xscale(i))  ! unscaled value
        call json%add('xvec('//trim(istr)//').scale',me%xscale(i))
    end do

    call json%print_file(trim(filename)//'_'//get_case_name()//'.json')

    call json%destroy()

    end subroutine write_optvars_to_file
!*****************************************************************************************

!*****************************************************************************************
    end module nrho_module
!*****************************************************************************************
