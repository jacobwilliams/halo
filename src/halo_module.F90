!*****************************************************************************************
!>
!  Main module with all the stuff to solve the Halo problem.

    module halo_module

    use parameters_module
    use fortran_astrodynamics_toolkit
    use ddeabm_module
    use nlesolver_module
    use numerical_differentiation_module
    use json_module
    use pyplot_module
    use bspline_module
    use config_file_module
    use splined_ephemeris_module
    use halo_utilities_module
    use moon_frame_module

    implicit none

    public  ! just make it all public for now

    type segment_data

        !! data for a segment

        real(wp) :: et_ref = zero  !! ephemeris time reference epoch
                                   !! (all times relative to this)

        real(wp)              :: t0 = zero              !! initial time [days]
        real(wp)              :: t0_scale = one         !! opt var scale for `t0`
        real(wp)              :: t0_dpert = 1.0e-5_wp   !! opt var dpert for `t0` [scaled]

        real(wp),dimension(6) :: x0_rotating = zero !! initial state [km, km/s], rotating frame
        real(wp),dimension(6) :: x0_rotating_scale = one !! opt var scale for `x0_rotating`
        real(wp),dimension(6) :: x0_rotating_dpert = [1.0e-2_wp,&
                                                      1.0e-2_wp,&
                                                      1.0e-2_wp,&
                                                      1.0e-5_wp,&
                                                      1.0e-5_wp,&
                                                      1.0e-5_wp]  !! opt var dperts for `x0_rotating` [scaled]

        real(wp),dimension(6) :: x0 = zero          !! initial state [km, km/s], inertial frame
        real(wp)              :: tf = zero          !! final time [days]
        real(wp),dimension(6) :: xf = zero          !! final state [km, km/s], inertial frame

        real(wp),dimension(6) :: xf_rotating = zero      !! final state [km, km/s], rotating frame
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

        ! These can be pointers that are pointing to the global ones in the mission,
        ! Or, when using OpenMP, they are allocated in each segment.
        type(geopotential_model_pines),pointer :: grav => null() !! central body geopotential model
        class(jpl_ephemeris),pointer :: eph => null()  !! the ephemeris
        type(moon_frame_interpolater),pointer :: moon_pa => null() !! the moon_pa frame interpolater
        logical :: pointmass_central_body = .false.
        logical :: include_pointmass_earth = .true. !! if true, earth is included as a pointmass in the force model
        logical :: include_pointmass_sun = .true. !! if true, sun is included as a pointmass in the force model
        logical :: include_pointmass_jupiter = .false. !! if true, jupiter is included as a pointmass in the force model
        procedure(third_body_grav_f),pointer,nopass :: third_body_gravity => null() !! procedure used to compute 3rd body gravity

        ! for saving the trajectory for plotting:
        type(trajectory) :: traj_inertial  !! in the inertial frame (j2000-moon)
        type(trajectory) :: traj_rotating  !! in the rotating frame (moon-earth, moon-centered)
        type(trajectory) :: traj_se_rotating  !! in the rotating frame (sun-earth, earth-centered)

        contains

        procedure,public :: set_input   => set_segment_inputs
        procedure,public :: get_inputs  => get_segment_inputs
        procedure,public :: get_outputs => get_segment_outputs
        procedure,public :: set_outputs => set_segment_outputs
        procedure,public :: propagate   => propagate_segment

        procedure :: cache
        procedure :: uncache
        procedure :: put_data_in_segment

    end type segment

    type,extends(numdiff_type) :: mission_type

        !! the mission [this is a `numdiff_type` for
        !! organizational purposes.... rethink this maybe ...

        character(len=:),allocatable :: initial_guess_from_file !! to read the initial guess from a JSOn file

        logical :: use_splined_ephemeris = .false. !! if true, the ephemeris is splined
        real(wp) :: dt_spline_sec = 3600.0_wp !! time step in seconds for spline step [1 hr]

        logical :: pointmass_central_body = .false. !! if true, the central body is a pointmass (moon).
                                                   !! otherwise, the `grav` model is used.
        type(geopotential_model_pines),pointer :: grav => null() !! central body geopotential model [global]
        class(jpl_ephemeris),pointer :: eph => null()            !! the ephemeris [global]
        type(moon_frame_interpolater),pointer :: moon_pa => null() !! the moon_pa frame interpolater [global]
        logical :: include_pointmass_earth = .true. !! if true, earth is included as a pointmass in the force model
        logical :: include_pointmass_sun = .true. !! if true, sun is included as a pointmass in the force model
        logical :: include_pointmass_jupiter = .false. !! if true, jupiter is included as a pointmass in the force model

        type(segment),dimension(:),allocatable :: segs  !! the list of segments

        real(wp),dimension(:),allocatable :: xscale  !! opt var scale factors
        real(wp),dimension(:),allocatable :: dpert_  !! opt var dpert factors [rename because already in base class]
        real(wp),dimension(:),allocatable :: fscale  !! function scale factors
        character(len=20),dimension(:),allocatable :: xname !! opt var names

        character(len=:),allocatable :: N_or_S      !! 'N' or 'S'
        character(len=:),allocatable :: L1_or_L2    !! 'L1' or 'L2'

        logical :: solve = .true. !! run the nle solver to solve the problem
        logical :: generate_plots = .true.  !! to generate the python plots
        logical :: generate_trajectory_files = .true.  !! to export the txt trajectory files.
        logical :: generate_json_trajectory_file = .false. !! to export the JSON trajectory file for the solution
                                                           !! (contains inertial and rotating data). This one is used
                                                           !! by the PyVista python script for plotting
        logical :: generate_guess_and_solution_files = .true.  !! to export the json guess and solution files.
        logical :: generate_kernel = .false.  !! to generate a spice kernel (bsp) of the solution
                                              !! [this requires the external mkspk SPICE tool]
        logical :: generate_defect_file = .false. !! generate the file that shows the pos/vel
                                                  !! constraint defects for the solution
        logical :: generate_eclipse_files = .false. !! generate the eclipse data file for the solution
        logical :: run_pyvista_script = .false. !! run the pyvista script to generate the interacdtive 3d plot

        logical :: generate_rp_ra_file = .false. !! to generate the rp_ra file

        real(wp) :: r_eclipse_bubble = 0.0_wp !! radius of the "eclipse bubble" [km]
        real(wp) :: eclipse_dt_step = 3600.0_wp !! dense time step output for eclipse calculations in sec [1 hour]
        integer :: eclipse_filetype = 1 !! type of eclipse file: 1=csv or 2=json

        integer :: epoch_mode = 1 !! 1 - calendar date was specified, 2 - ephemeris time was specified
        integer :: year   = 0 !! epoch of first point (first periapsis crossing)
        integer :: month  = 0 !! epoch of first point (first periapsis crossing)
        integer :: day    = 0 !! epoch of first point (first periapsis crossing)
        integer :: hour   = 0 !! epoch of first point (first periapsis crossing)
        integer :: minute = 0 !! epoch of first point (first periapsis crossing)
        real(wp) :: sec   = zero !! epoch of first point (first periapsis crossing)
        real(wp) :: et_ref = zero  !! ephemeris time reference epoch [sec]
                                   !! (all times relative to this).
                                   !! this is computed from the year,month,day,hour,minute,sec values
                                   !! if `epoch_mode=1`

        integer :: n_revs = 10  !! Number of revs in the Halo.

        real(wp) :: rtol = 1.0e-12_wp !! integrator rtol
        real(wp) :: atol = 1.0e-12_wp !! integrator atol
        real(wp) :: nlesolver_tol = 1.0e-6_wp !! nlesolver tol

        real(wp),dimension(:),allocatable :: initial_r
            !! if `fix_initial_r` is True, this can be used
            !! to specify the initial `r` (moon-centered earth-moon rotating frame)

        character(len=:),allocatable :: ephemeris_file !! the JPL ephemeris file to load
        ! [note: these are build using the get_third_party script in FAT.
        !  the files are platform (and maybe compiler specific)
        !  'data/eph/JPLEPH.421'               !! JPL DE421 ephemeris file [mac,gfortran]
        !  'data/eph/JPLEPH_intel_mac.421'     !! [mac,ifort]
        !  'data/eph/JPLEPH_windows_ifort.421' !! [windows,ifort]

        ! environment and integrator parameters:

        character(len=:),allocatable :: gravfile  !! spherical harmonic gravity coeff file (Moon)
        !! example: 'data/grav/gggrx_0020pm_sha.tab'

        character(len=:),allocatable ::  moon_pa_file !! file for moon_pa frame.
        !! example: `data/moon_pa_2000_2100.csv`

        character(len=:),allocatable :: patch_point_file   !! Halo CR3BP patch point solution file
        !! example: 'data/L2_halos.json'
        logical :: patch_point_file_is_periapsis = .false. !! if the state in the patch point file are periapsis states
                                                           !! (false means they are apoapsis states.
                                                           !! the ones in L2_halos.json are apoapsis states)

        ! the initial guess for the Halo:
        ! This is for the "forward/backward" version of the transcription from the paper.
        type(patch_point) :: periapsis    !! Patchpoint Periapse State
        type(patch_point) :: quarter      !! Patchpoint 1/4 Rev State
        type(patch_point) :: apoapsis     !! Patchpoint 1/2 Rev State
        real(wp) :: period  = 0.0_wp !! Halo period [days]
        real(wp) :: period8 = 0.0_wp !! 1/8 of Halo period [days]

        ! optional problem inputs
        ! [can remove some of the variables from the optimization problem]
        logical :: fix_initial_time = .false. !! to fix the initial epoch in the mission
        logical :: fix_initial_r = .false. !! fix the initial periapsis position vector (x,y,z)
        integer :: fix_ry_at_end_of_rev = 0 !! fix ry at the end of the specified rev (at periapsis)
                                            !! not used if `<= 0`.
        logical :: fix_final_ry_and_vx = .false. !! fix ry and vx at the end of the last rev (at periapsis)

        logical :: constrain_initial_rdot = .false. !! impose an rdot=0 constraint on the initial point

        integer :: solver_mode = 1  !! use sparse or dense solver:
                                    !!
                                    !! * 1 - dense (LAPACK)
                                    !! * 2 - sparse (LSQR)
                                    !! * 3 - sparse (LUSOL)
                                    !! * 4 - sparse (LMSR)
                                    !! * 5 - sparse (qr_mumps)

    contains

        procedure,public :: init => initialize_the_mission
        procedure,public :: plot => plot_trajectory
        procedure,public :: export_trajectory_json_file

        procedure,public :: write_optvars_to_file
        procedure,public :: constraint_violations
        procedure,public :: print_constraint_defects
        procedure,public :: generate_eclipse_data
        procedure :: get_problem_arrays
        procedure :: put_x_in_segments
        procedure :: get_x_from_json_file
        procedure :: get_scales_from_segs
        procedure :: segs_to_propagate
        procedure :: get_sparsity_pattern
        procedure :: define_problem_size
        procedure :: get_case_name
        procedure :: generate_patch_points
        procedure :: update_epoch
        procedure :: export_rp_ra_json_file

    end type mission_type

    type,extends(nlesolver_type) :: my_solver_type

        !! the solver class, which contains
        !! an instance of the mission

        type(mission_type) :: mission !! the Halo mission

    contains

        procedure,public :: init => initialize_the_solver

        procedure :: read_config_file

    end type my_solver_type

    abstract interface
        subroutine third_body_grav_f(r,rb,mu,acc) !! for compute third-body gravity
            import :: wp
            implicit none
            real(wp),dimension(3),intent(in)  :: r   !! satellite position vector [km]
            real(wp),dimension(3),intent(in)  :: rb  !! third-body position vector [km]
            real(wp),intent(in)               :: mu  !! third-body gravitational parameter [km^3/s^2]
            real(wp),dimension(3),intent(out) :: acc !! gravity acceleration vector [km/s^2]
        end subroutine third_body_grav_f
    end interface

    public :: halo_solver_main !! main program
    public :: read_epoch ! also used by optimizer app

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main program to solve the Halo targeting problem.
!
!### Author
!  * Jacob Williams : Sept. 2017

    subroutine halo_solver_main(config_file_name,debug)

    use parameters_module
!$  use omp_lib

    implicit none

    character(len=*),intent(in) :: config_file_name  !! the config file to read
    logical,intent(in) :: debug !! for debugging prints

    type(my_solver_type) :: solver  !! an instance of the solver that we will use
    real(wp),dimension(:),allocatable :: x  !! solver opt vars vector ["forward-backward" formulation]
    integer :: m  !! number of functions
    real(wp),dimension(:),allocatable :: f  !! function vector (constraint violations)
    real(wp) :: tstart, tend  !! for timing
    real(wp) :: tstart_cpu, tend_cpu  !! for timing
    integer :: istat
    character(len=:),allocatable :: message  !! Text status message from solver
    integer :: n_segs, iseg
    real(wp),dimension(6) :: x_rotating
    character(len=:),allocatable :: mkspk_input, bsp_output, mkspk_setup, pck_output  !! filenames for mkspk
!$  integer :: tid, nthreads

    write(*,'(A)') ''
    write(*,'(A)') ' * HALO start'

!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!$
!$  tid = omp_get_thread_num()
!$
!$  if (tid == 0) then
!$      nthreads = omp_get_num_threads()
!$      write(*,'(A,1X,I3)') ' * Number of OMP threads: ', OMP_get_num_threads()
!$  end if
!$
!$OMP END PARALLEL

    if (debug) then
        write(*,*) ''
        write(*,*) '----------------------'
        write(*,*) 'Initializing...'
        write(*,*) '----------------------'
        write(*,*) ''
    end if

    call solver%init(config_file_name,x)  ! initialize the solver & mission (and generate the initial guess)

    if (solver%mission%generate_plots) &
        call solver%mission%plot('guess', draw_trajectory=.true.)    ! plot the initial guess
    if (solver%mission%generate_guess_and_solution_files) &
        call solver%mission%write_optvars_to_file('guess',x)    ! write guess to a file

    !....debugging....
    ! if (solver%mission%generate_trajectory_files) &
    !     call solver%mission%plot('guess',&
    !             draw_trajectory = .false., &
    !             export_trajectory=solver%mission%generate_trajectory_files)
    !....debugging....

    if (debug) then
        call solver%mission%define_problem_size(m=m)
        allocate(f(m))
        call solver%mission%constraint_violations(x,f)
        write(*,*) ''
        write(*,*) '----------------------'
        write(*,*) 'Initial Guess...'
        write(*,*) '----------------------'
        write(*,*) ''
        write(*,'(A/,*(F30.16/))') 'x:      ', x * solver%mission%xscale    ! unscaled values
        write(*,*) ''
        write(*,'(A/,*(F30.16/))') 'f:      ', f    ! scaled values
        write(*,*) ''
    end if

    if (debug) then
        write(*,*) 'INITIAL GUESS:'
        call solver%mission%define_problem_size(n_segs=n_segs)
        do iseg = 1, n_segs
            call solver%mission%segs(iseg)%get_inputs(x0_rotating=x_rotating)
            write(*,'(I5, *(F15.6,1X))') iseg, x_rotating
        end do
    end if

    if (solver%mission%solve) then
        if (debug) then
            write(*,*) ''
            write(*,*) '----------------------'
            write(*,*) 'Solving...'
            write(*,*) '----------------------'
            write(*,*) ''
        end if
        write(*,'(A)') ' * Solving'
        write(*,*) ''

        call cpu_time(tstart_cpu)
!$      tstart = omp_get_wtime()
        call solver%solve(x)  ! call the solver
        call solver%status(istat=istat, message=message)
        call cpu_time(tend_cpu)
!$      tend = omp_get_wtime()
        call solver%mission%put_x_in_segments(x) ! populate segs with solution
    else
        call cpu_time(tstart_cpu)
        call cpu_time(tend_cpu)
        message = 'Not solved'
    end if

    if (debug) then
        write(*,*) ''
        write(*,*) '----------------------'
        write(*,*) 'Solution...'
        write(*,*) '----------------------'
        write(*,*) ''
    end if

    write(*,*) ''
    write(*,'(A)') ' * Status: '//message
    write(*,'(A,1x,F10.3,1x,a)') ' * Elapsed cpu_time: ', (tend_cpu-tstart_cpu), 'sec'
!$  write(*,'(A,1x,F10.3,1x,a)') ' * OMP wall time   : ', (tend-tstart), 'sec'

    if (solver%mission%generate_guess_and_solution_files) then
        write(*,'(A)') ' * Generate solution file'
        call solver%mission%write_optvars_to_file('solution',x) ! write solution to a file
    end if

    ! export solution to plot or trajectory file
    if (solver%mission%generate_trajectory_files .or. solver%mission%generate_plots) then
        write(*,'(A)') ' * Export solution trajectory'
        call solver%mission%plot('solution',&
                draw_trajectory=solver%mission%generate_plots, &
                export_trajectory=solver%mission%generate_trajectory_files)
    end if

    if (solver%mission%generate_kernel) then
        if (.not. solver%mission%generate_trajectory_files) then
            error stop 'error: kernel generation requires the trajectory file to be exported'
        else
            write(*,*) '* Generate BSP kernel'

            mkspk_input = 'solution_'//solver%mission%get_case_name()//'.txt'
            bsp_output  = 'solution_'//solver%mission%get_case_name()//'.bsp'
            mkspk_setup = 'mkspk_'//solver%mission%get_case_name()//'.txt'
            pck_output  = 'solution_'//solver%mission%get_case_name()//'.tpc.json'

            !note: have to delete the kernel if it is already there or mkspk will fail.
            call delete_file(bsp_output)

            call generate_mkspk_setup_file(mkspk_setup)
            call generate_pck_file(pck_output)
            call execute_command_line(mkspk_path//' -setup '//mkspk_setup//' -input '//mkspk_input//' -output '//bsp_output)
        end if
    end if

    if (debug) then
        write(*,*) 'SOLUTION:'
        call solver%mission%define_problem_size(n_segs=n_segs)
        do iseg = 1, n_segs
            call solver%mission%segs(iseg)%get_inputs(x0_rotating=x_rotating)
            write(*,'(I5, *(F15.6,1X))') iseg, x_rotating
        end do
    end if

    if (solver%mission%generate_defect_file) then
        write(*,'(A)') ' * Generate defect file'
        call solver%mission%print_constraint_defects('defects_'//&
                                                     solver%mission%get_case_name()//&
                                                     '.csv')
    end if

    if (solver%mission%generate_json_trajectory_file) then
        call solver%mission%export_trajectory_json_file('traj_'//solver%mission%get_case_name())
    end if

    if (solver%mission%generate_eclipse_files) then
        call solver%mission%generate_eclipse_data('eclipse', &
                                                  filetype = solver%mission%eclipse_filetype)
    end if

    if (solver%mission%run_pyvista_script) then
        write(*,'(A)') ' * Run PyVista script'
        !mkspk_input = 'solution_'//solver%mission%get_case_name()//'.txt'
        mkspk_input = 'traj_'//solver%mission%get_case_name()//'.json'
        call execute_command_line('python ./python/plot_utilities.py '//mkspk_input, wait=.false.)
    end if

    write(*,'(A)') ' * HALO done'
    write(*,'(A)') ''

    end subroutine halo_solver_main
!*****************************************************************************************

!*****************************************************************************************
!>
!  Delete a file if it exists

    subroutine delete_file(name)

    character(len=*),intent(in) :: name

    integer :: istat, iunit

    open(file=name,newunit=iunit,status='OLD',iostat=istat)
    if (istat==0) then
        write(*,'(A)') ' * Deleting existing file: '//trim(name)
        close(iunit,status='DELETE',iostat=istat)
    end if

    end subroutine delete_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Generate a PCK (JSON version) to go with the BSP file.
!  Note: this kernel can be read by `jsonspice` + `SpiceyPy`.
!
!### Reference
!  * [NAIF Integer ID codes](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/naif_ids.html)

    subroutine generate_pck_file(filename)

    character(len=*),intent(in) :: filename

    type(json_file) :: json

    call json%initialize(compress_vectors=.true.)
    call json%add('+NAIF_BODY_NAME', object_name)
    call json%add('+NAIF_BODY_CODE', object_id)
    call json%add('BODY'//int_to_string(object_id)//'_RADII', [0.001_wp, 0.001_wp, 0.001_wp])
    call json%print(filename)
    call json%destroy()

    end subroutine generate_pck_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Generate the mkspk setup file. See `kernel/mkspk` for an example.
!
!  See also: [MKSPK User's Guide](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/ug/mkspk.html)

    subroutine generate_mkspk_setup_file(filename)

    character(len=*),intent(in) :: filename

    integer :: iunit, istat

    open(newunit=iunit, file=trim(filename), status='REPLACE', iostat=istat)
    if (istat/=0) error stop 'error opening file: '//trim(filename)

    write(iunit, '(A)') "\begindata"
    write(iunit, '(A)') ""
    write(iunit, '(A)') "INPUT_DATA_TYPE   = 'STATES'"
    write(iunit, '(A)') "DATA_ORDER        = 'EPOCH X Y Z VX VY VZ'"
    write(iunit, '(A)') "DATA_DELIMITER    = ';'"
    write(iunit, '(A)') "TIME_WRAPPER      = '# ETSECONDS'"
    write(iunit, '(A)') "CENTER_ID         = 301"
    write(iunit, '(A)') "CENTER_NAME       = 'MOON'"
    write(iunit, '(A)') "REF_FRAME_NAME    = 'J2000'"
    write(iunit, '(A)') "INPUT_DATA_UNITS  = ('ANGLES=DEGREES' 'DISTANCES=km')"
    write(iunit, '(A)') "PRODUCER_ID       = 'HALO'"
    write(iunit, '(A)') "OUTPUT_SPK_TYPE   = "//int_to_string(output_spk_type)
    write(iunit, '(A)') "POLYNOM_DEGREE    = "//int_to_string(polynom_degree)
    write(iunit, '(A)') "SEGMENT_ID        = '"//segment_id//"'"
    write(iunit, '(A)') "OBJECT_ID         = "//int_to_string(object_id)
    write(iunit, '(A)') "OBJECT_NAME       = '"//trim(object_name)//"'"
    write(iunit, '(A)') "LEAPSECONDS_FILE  = '"//leapseconds_file//"'"
    close(iunit, iostat=istat)

    end subroutine generate_mkspk_setup_file
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
!
!### History
!  * JW : 8/1/2022 : added option to fix initial time

    pure subroutine define_problem_size(me,n,m,n_segs,n_nonzero,full_problem)

    implicit none

    class(mission_type),intent(in) :: me
    integer,intent(out),optional :: n          !! number of opt vars for the solver
    integer,intent(out),optional :: m          !! number of equality constraints for the solver
    integer,intent(out),optional :: n_segs     !! number of segments
    integer,intent(out),optional :: n_nonzero  !! number of nonzero elements in the Jacobian
    logical,intent(in),optional :: full_problem !! if true (default) always return the results for the full problem
                                                !! (i.e., without fixing any of the variables)

    if (present(n))      n      = 28 * me%n_revs + 7
    if (present(m))      m      = 24 * me%n_revs
    if (present(n_segs)) n_segs = 8 * me%n_revs

    ! there are 4 blocks of nonzeros per rev
    ! (each block contains 84 elements)
    if (present(n_nonzero)) then
        n_nonzero = me%n_revs * (84*4)
        if (me%constrain_initial_rdot) then
            if (me%fix_initial_r) then ! add elements (rdot is function of initial r,v)
                n_nonzero = n_nonzero + 3  ! only v is in the opt var vec
            else
                n_nonzero = n_nonzero + 6  ! r,v in opt var vec
            end if
        end if
    end if
    ! optional rdot constraint
    if (me%constrain_initial_rdot) then
        if (present(m)) m = m + 1  ! add one function (the last row)
    end if

    if (present(full_problem)) then
        if (full_problem) return ! don't do the stuff below
    end if

    if (me%fix_initial_time) then
        if (present(n)) n = n - 1  ! remove the t0 optimization variable
        if (present(n_nonzero)) n_nonzero = n_nonzero - 6 ! remove the first column of the jacobian
    end if

    if (me%fix_initial_r) then
        if (present(n)) n = n - 3  ! remove the three optimization variables
        if (present(n_nonzero)) n_nonzero = n_nonzero - 3*6 ! remove columns 2,3,4 of the jacobian
    end if

    ! note: assuming there are more than 2 revs for these
    if (me%fix_ry_at_end_of_rev > 0) then
        if (me%n_revs<3) error stop 'at least 3 revs are required for fix_ry_at_end_of_rev'
        if (me%fix_ry_at_end_of_rev >= me%n_revs) &
            error stop 'fix_ry_at_end_of_rev must be < number of revs'
        if (present(n)) n = n - 1 ! remove the optimization variable
        if (present(n_nonzero)) n_nonzero = n_nonzero - 12
    end if

    if (me%fix_final_ry_and_vx) then
        if (me%n_revs<3) error stop 'at least 3 revs are required for fix_final_ry_and_vx'
        if (present(n)) n = n - 2 ! remove the two optimization variable
        if (present(n_nonzero)) n_nonzero = n_nonzero - 2*6
    end if

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

    subroutine initialize_the_solver(me,config_file_name,x)

    implicit none

    class(my_solver_type),intent(inout) :: me
    character(len=*),intent(in) :: config_file_name !! the config file to read
    real(wp),dimension(:),allocatable,intent(out),optional :: x !! initial guess

    integer :: n          !! number of opt vars for the solver
    integer :: m          !! number of equality constraints for the solver
    logical :: status_ok  !! status flag for solver initialization
    integer :: istat      !! status code from solver
    integer,dimension(:),allocatable :: irow,icol !! sparsity pattern

    ! note: we don't know the problem size
    ! until we read the config file.

    ! first we have to read the config file:
    ! this will populate some of the mission variables
    call me%read_config_file(config_file_name)

    ! initialize the mission
    call me%mission%init(x)
    call me%mission%define_problem_size(n,m)

    ! initialize the solver:
    select case (me%mission%solver_mode)
    case(1)
        ! dense - uses lapack to solve the linear system
        call me%initialize(     n                = n,            &
                                m                = m,            &
                                max_iter         = 100,          & ! maximum number of iteration
                                func             = halo_func,    &
                                grad             = halo_grad,    &
                                tol              = me%mission%nlesolver_tol,    & ! tolerance
                                ! step_mode        = 4,            & ! 3-point "line search" (2 intervals)
                                ! n_intervals      = 2,            & ! number of intervals for step_mode=4
                                ! alpha_min = 0.2_wp, &
                                ! alpha_max = 0.8_wp, &
                                step_mode = 2,            & ! backtracking "line search"
                                alpha_min = 0.1_wp, &
                                alpha_max = 1.0_wp, &
                                use_broyden      = .false.,      & ! broyden update
                                ! use_broyden=.true.,broyden_update_n=10, & ! ... test ...
                                export_iteration = halo_export   )

    case(5)
        ! this is using the qr_mumps solver as a user-defined solver to nlesolver-fortran.
        ! the solver is defined in qrm_solver. you must use the WITH_QRMUMPS preprocessor
        ! directive to use this method and link the code with the appropriate libraries.

        call me%mission%get_sparsity_pattern(irow,icol) ! it's already been computed, but for now, just compute it again for this call
        call me%initialize(     n                = n,            &
                                m                = m,            &
                                max_iter         = 100,          & ! maximum number of iteration
                                func             = halo_func,    &
                                grad_sparse      = halo_grad_sparse,    &
                                tol              = me%mission%nlesolver_tol,    & ! tolerance
                                step_mode        = 4,            & ! 3-point "line search" (2 intervals)
                                n_intervals      = 2,            & ! number of intervals for step_mode=4
                                alpha_min        = 0.2_wp, &
                                alpha_max        = 1.0_wp, &
                                ! step_mode = 2,  & ! backtracking "line search"  seems to not work as well for large problems. why?
                                ! alpha_min = 0.1_wp, &
                                ! alpha_max = 1.0_wp, &
                                use_broyden      = .false.,      & ! broyden update
                                ! use_broyden=.true.,broyden_update_n=4, & ! ... test ...
                                sparsity_mode    = me%mission%solver_mode, &  ! use a sparse solver
                                custom_solver_sparse = qrm_solver, &  ! the qr_mumps solver wrapper
                                irow             = irow, &  ! sparsity pattern
                                icol             = icol, &
                                export_iteration = halo_export   )
    case (2:4)
        ! varions sparse options available in nlesolver-fortran
        call me%mission%get_sparsity_pattern(irow,icol) ! it's already been computed, but for now, just compute it again for this call
        call me%initialize(     n                = n,            &
                                m                = m,            &
                                max_iter         = 100,          & ! maximum number of iteration
                                func             = halo_func,    &
                                grad_sparse      = halo_grad_sparse,    &
                                tol              = me%mission%nlesolver_tol,    & ! tolerance
                                ! step_mode = 1,& ! TEST TEST TEST
                                ! alpha = 1.0_wp,&
                                step_mode        = 4,            & ! 3-point "line search" (2 intervals)
                                n_intervals      = 2,            & ! number of intervals for step_mode=4
                                use_broyden      = .false.,      & ! broyden update
                                !use_broyden=.true.,broyden_update_n=10, & ! ... test ...
                                export_iteration = halo_export,  &
                                sparsity_mode = me%mission%solver_mode, &  ! use a sparse solver
                                atol          = 1.0e-12_wp,&  ! relative error in definition of `A`
                                btol          = 1.0e-12_wp,&  ! relative error in definition of `b`
!                                damp          = 0.00001_wp, & !  TEST: LSQR damp factor !
                                damp          = 0.0_wp, & !  TEST: LSQR damp factor !
   !                             damp          = 0.1_wp, & !  TEST: LSQR damp factor !
                                itnlim        = 1000000, &  ! max iterations
                                irow          = irow, &  ! sparsity pattern
                                icol          = icol, &
                                lusol_method = 0     ) ! test

    case default
        error stop 'invalid solver_mode'
    end select

    call me%status(istat=istat)
    status_ok = istat == 0
    if (.not. status_ok) then
        write(*,*) 'istat = ', istat
        error stop 'error in initialize_the_solver'
    end if

    if (allocated(me%mission%initial_guess_from_file)) then
        if (me%mission%initial_guess_from_file /= '') then
            !TODO: add some error checking here !
            write(*,'(A)') ' * Reading initial guess from file: '//me%mission%initial_guess_from_file
            call me%mission%get_x_from_json_file(x) ! get solution from the file
            call me%mission%put_x_in_segments(x) ! populate segs with solution
        end if
    end if

    end subroutine initialize_the_solver
!*****************************************************************************************

!*****************************************************************************************
!>
!  Custom solver that uses QR_MUMPS

    subroutine qrm_solver(me,n_cols,n_rows,n_nonzero,irow,icol,a,b,x,istat)
#ifdef WITH_QRMUMPS
        use dqrm_mod
#endif
        implicit none

        class(nlesolver_type),intent(inout) :: me
        integer,intent(in) :: n_cols !! `n`: number of columns in A.
        integer,intent(in) :: n_rows !! `m`: number of rows in A.
        integer,intent(in) :: n_nonzero !! number of nonzero elements of A.
        integer,dimension(n_nonzero),intent(in) :: irow, icol !! sparsity pattern (size is `n_nonzero`)
        real(wp),dimension(n_nonzero),intent(in) :: a !! matrix elements (size is `n_nonzero`)
        real(wp),dimension(n_rows),intent(in) :: b !! right hand side (size is `m`)
        real(wp),dimension(n_cols),intent(out) :: x !! solution (size is `n`)
        integer,intent(out) :: istat !! status code (=0 for success)

#ifdef WITH_QRMUMPS
        type(dqrm_spmat_type) :: qrm_spmat

        ! hack because we have to point to them ! can we avoid this ?? <-----
        integer,dimension(:),allocatable,target :: irow_, icol_
        real(wp),dimension(:),allocatable,target :: a_, r_
        !--------------------------------------------------------------------

        allocate(irow_, source=irow)
        allocate(icol_, source=icol)
        allocate(a_ , source=a)
        allocate(r_ , source=b)

        call qrm_init()

        ! initialize the matrix data structure.
        call qrm_spmat_init(qrm_spmat)

        qrm_spmat%m   =  n_rows
        qrm_spmat%n   =  n_cols
        qrm_spmat%nz  =  n_nonzero
        qrm_spmat%irn => irow_
        qrm_spmat%jcn => icol_
        qrm_spmat%val => a_

        call qrm_spmat_gels(qrm_spmat, r_, x)

        istat = 0 ! how to get status code?

#else
    error stop 'This code was not build with QR_MUMPS'
#endif
    end subroutine qrm_solver

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
    real(wp) :: rdot      !! initial d(rmag)/dt (km/s) in moon-centered rotating frame

    if (present(x)) then

        x = -huge(1.0_wp) ! just in case we miss one

        i = 0 ! initialize the index. will be updated by fill_vector
        iseg = 0
        do irev = 1, me%n_revs

            if (irev==1) then

                ! the first one has an extra opt point at the initial periapsis passage

                if (.not. me%fix_initial_time) &
                    call fill_vector(x, me%segs(iseg+1)%data%t0, i)

                if (me%fix_initial_r) then
                    call fill_vector(x, me%segs(iseg+1)%data%x0_rotating(4:6), i) ! v only
                else
                    call fill_vector(x, me%segs(iseg+1)%data%x0_rotating, i) ! r,v
                end if

                call fill_vector(x, me%segs(iseg+2)%data%t0, i)
                call fill_vector(x, me%segs(iseg+2)%data%x0_rotating, i)

                call fill_vector(x, me%segs(iseg+4)%data%t0, i)
                call fill_vector(x, me%segs(iseg+4)%data%x0_rotating, i)

                call fill_vector(x, me%segs(iseg+6)%data%t0, i)
                call fill_vector(x, me%segs(iseg+6)%data%x0_rotating, i)

                call fill_vector(x, me%segs(iseg+8)%data%t0, i)

                if (me%fix_ry_at_end_of_rev == 1) then
                    call fill_vector(x, me%segs(iseg+8)%data%x0_rotating([1,3,4,5,6]), i)
                else
                    call fill_vector(x, me%segs(iseg+8)%data%x0_rotating, i)
                end if

                ! for next rev:
                iseg = iseg + 8

            else

                call fill_vector(x, me%segs(iseg+2)%data%t0, i)
                call fill_vector(x, me%segs(iseg+2)%data%x0_rotating, i)

                call fill_vector(x, me%segs(iseg+4)%data%t0, i)
                call fill_vector(x, me%segs(iseg+4)%data%x0_rotating, i)

                call fill_vector(x, me%segs(iseg+6)%data%t0, i)
                call fill_vector(x, me%segs(iseg+6)%data%x0_rotating, i)

                call fill_vector(x, me%segs(iseg+8)%data%t0, i)

                if (me%fix_ry_at_end_of_rev == irev) then
                    call fill_vector(x, me%segs(iseg+8)%data%x0_rotating([1,3,4,5,6]), i)
                elseif (irev==me%n_revs .and. me%fix_final_ry_and_vx) then
                    call fill_vector(x, me%segs(iseg+8)%data%x0_rotating([1,3,5,6]), i)
                else
                    call fill_vector(x, me%segs(iseg+8)%data%x0_rotating, i)
                end if

                ! for next rev:
                iseg = iseg + 8

            end if

        end do

        if (any(x == -huge(1.0_wp))) then
            error stop 'error: not all x values were set'
        end if

        ! scale the x vector:
        x = x / me%xscale

    end if

    if (present(f)) then

        ! f = [xf1-xf2, xf3-xf4, xf5-xf6, xf7-xf8, ... ]

        i = 0
        iseg = 0
        do irev = 1, me%n_revs

            call fill_vector(f, -(me%segs(iseg+1)%data%xf_rotating - me%segs(iseg+2)%data%xf_rotating), i)
            call fill_vector(f, -(me%segs(iseg+3)%data%xf_rotating - me%segs(iseg+4)%data%xf_rotating), i)
            call fill_vector(f, -(me%segs(iseg+5)%data%xf_rotating - me%segs(iseg+6)%data%xf_rotating), i)
            call fill_vector(f, -(me%segs(iseg+7)%data%xf_rotating - me%segs(iseg+8)%data%xf_rotating), i)

            iseg = iseg + 8

        end do
        if (me%constrain_initial_rdot) then
            rdot = compute_rdot(me%segs(1)%data%x0)
            call fill_vector(f, rdot, i)  ! last one is the rdot=0 constraint
        end if

        if (size(f) /= size(me%fscale)) error stop 'error: f and fscale are not the same size'

        !scale the f vector:
        f = f / me%fscale

    end if

    end subroutine get_problem_arrays
!*****************************************************************************************

!*****************************************************************************************
!>
!  Print the `r` and `v` constraint defect norms for each segment constraint.

    subroutine print_constraint_defects(me, filename)

        class(mission_type),intent(in) :: me
        character(len=*),intent(in) :: filename !! csv file to write to

        real(wp),dimension(:),allocatable :: f !! constraint violations
        integer :: irev !! rev number
        integer :: i, j
        integer :: m !! number of functions
        integer :: istat , iunit

        open(newunit=iunit, file=filename, status='REPLACE', iostat=istat)
        if (istat/=0) error stop 'error opening '//filename

        call me%define_problem_size(m=m)
        allocate(f(m))
        call me%get_problem_arrays(f=f)

        ! f = [xf1-xf2, xf3-xf4, xf5-xf6, xf7-xf8, ... ] in rotating frame

        i = 0
        write(iunit,'(A26,A1,A26)') 'rerr (km)', ',', 'verr (km/s)'
        do irev = 1, me%n_revs
            do j = 1, 4
                i = i + 1
                write(iunit,'(1P,E26.16,A1,E26.16)') norm2(f(i:i+2)), ',', norm2(f(i+2:i+5))
            end do
        end do

        if (me%constrain_initial_rdot) then  ! also the rdot one
           write(iunit,'(1P,E26.16,A1,E26.16)') f(m), ',', 0.0_wp
        end if

        close(iunit, iostat=istat)

    end subroutine print_constraint_defects
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the JSON solution file and put the `x` vector in the mission.
!  Can be used to restart a solution from a previous run
!  (e.g., with different settings, as long as the fundamental problem isn't changed).
!
!@note No error checking here to make sure the file is consistent with the current mission!

    subroutine get_x_from_json_file(me,x)

    class(mission_type),intent(inout) :: me
    real(wp),dimension(:),intent(out),allocatable :: x !! scaled `x` vector

    real(wp),dimension(:),allocatable :: xscale !! scale factors
    type(json_core) :: json
    type(json_file) :: f
    type(json_value),pointer :: xvec, p_element
    integer :: n_children, i
    logical :: found

    if (allocated(me%initial_guess_from_file)) then

        ! read this file:
        ! {
        ! "xvec": [
        !   {
        !     "i": 1,
        !     "label": "SEG1 T0 (day)",
        !     "value": -0.96866230153337847E-3,
        !     "scale": 0.1E+1
        !   },
        !   ...

        call f%load(me%initial_guess_from_file)
        call f%get('xvec',xvec,found)
        if (.not. found) error stop 'invalid JSON solution file'
        call json%info(xvec,n_children=n_children) ! get size of x
        allocate(x(n_children))
        allocate(xscale(n_children))
        !TODO: should verify size of x compared to current problem.
        !      really should check that all the vars are the same
        ! get each element:
        do i = 1, n_children
            call json%get_child(xvec, i, p_element, found)
            call json%get(p_element,'value',x(i), found)
            if (.not. found) error stop 'could not find value in json file'
            call json%get(p_element,'scale',xscale(i), found)
            if (.not. found) error stop 'could not find scale in json file'
        end do
        x = x / xscale  ! return scaled vector
    else
        error stop 'error: the initial_guess_from_file has not been initialized'
    end if

    end subroutine get_x_from_json_file
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
    real(wp),dimension(:),allocatable :: x !! opt var vector (unscaled)
    integer :: n_segs  !! number of segments

    call me%define_problem_size(n_segs=n_segs)

    ! unscale the x vector:
    allocate(x(size(x_scaled)))
    x = x_scaled * me%xscale

    ! first extract data from the opt var vector and put it into the segments:
    i = 0  ! initialize index, will be updated by extract_vector
    iseg = 0
    do irev = 1, me%n_revs

        ! extract the values:

        if (irev==1) then

            ! the first one has an extra opt point at the initial periapsis passage

            if (me%fix_initial_time) then
                t0(1) = me%segs(1)%data%t0    ! not in x, just keep the current value
            else
                call extract_vector(t0(1), x, i)
            end if
            if (me%fix_initial_r) then
                x0_rotating(1:3,1) = me%segs(1)%data%x0_rotating(1:3) ! not in x, just keep the current r value
                call extract_vector(x0_rotating(4:6,1), x, i) ! v only
            else
                call extract_vector(x0_rotating(:,1), x, i)
            end if
            tf(1)               = t0(1) + me%period8

            call extract_vector(t0(2)               , x, i)
            call extract_vector(x0_rotating(:,2)    , x, i)
            tf(2)               = tf(1)

            call extract_vector(t0(4)               , x, i)
            call extract_vector(x0_rotating(:,4)    , x, i)
            tf(4)               = t0(4) - me%period8

            t0(3)               = t0(2)
            x0_rotating(:,3)    = x0_rotating(:,2)
            tf(3)               = tf(4)

            t0(5)               = t0(4)
            x0_rotating(:,5)    = x0_rotating(:,4)
            tf(5)               = t0(5) + me%period8

            call extract_vector(t0(6)               , x, i)
            call extract_vector(x0_rotating(:,6)    , x, i)
            tf(6)               = tf(5)

            call extract_vector(t0(8)               , x, i)
            if (me%fix_ry_at_end_of_rev == 1) then
                x0_rotating(2,8) = me%segs(8)%data%x0_rotating(2) ! not in x, just keep the current ry value
                call extract_vector(x0_rotating(1,8) , x, i)
                call extract_vector(x0_rotating(3:6,8) , x, i)
            else
                call extract_vector(x0_rotating(:,8)    , x, i)
            end if
            tf(8)               = t0(8) - me%period8

            t0(7)               = t0(6)
            x0_rotating(:,7)    = x0_rotating(:,6)
            tf(7)               = tf(8)

            ! put the data for this rev into the segments:
            do j = 1, 8
                call me%segs(iseg+j)%set_input(t0(j),tf(j),x0_rotating(:,j))
            end do

            ! for next rev:
            iseg = iseg + 8

        else

            t0(1)               =  t0(8)                ! inherits from the previous rev
            x0_rotating(:,1)    =  x0_rotating(:,8)
            tf(1)               =  t0(1) + me%period8

            call extract_vector(t0(2)               , x, i)
            call extract_vector(x0_rotating(:,2)    , x, i)
            tf(2)               =  tf(1)

            call extract_vector(t0(4)               , x, i)
            call extract_vector(x0_rotating(:,4)    , x, i)
            tf(4)               =  t0(4) - me%period8

            t0(3)               =  t0(2)
            x0_rotating(:,3)    =  x0_rotating(:,2)
            tf(3)               =  tf(4)

            t0(5)               =  t0(4)
            x0_rotating(:,5)    =  x0_rotating(:,4)
            tf(5)               =  t0(5) + me%period8

            call extract_vector(t0(6)               , x, i)
            call extract_vector(x0_rotating(:,6)    , x, i)
            tf(6)               =  tf(5)

            call extract_vector(t0(8)               , x, i)

            if (me%fix_ry_at_end_of_rev == irev) then
                x0_rotating(2,8) = me%segs(irev*8)%data%x0_rotating(2) ! not in x, just keep the current ry value
                call extract_vector(x0_rotating(1,8)   , x, i)
                call extract_vector(x0_rotating(3:6,8) , x, i)
            else if (irev==me%n_revs .and. me%fix_final_ry_and_vx) then
                x0_rotating(2,8) = me%segs(n_segs)%data%x0_rotating(2) ! not in x, just keep the current ry value
                x0_rotating(4,8) = me%segs(n_segs)%data%x0_rotating(4) ! not in x, just keep the current vx value
                call extract_vector(x0_rotating(1,8)    , x, i)
                call extract_vector(x0_rotating(3,8)    , x, i)
                call extract_vector(x0_rotating(5:6,8)  , x, i)
            else
                call extract_vector(x0_rotating(:,8)    , x, i)
            end if
            tf(8)               =  t0(8) - me%period8

            t0(7)               =  t0(6)
            x0_rotating(:,7)    =  x0_rotating(:,6)
            tf(7)               =  tf(8)

            ! put the data for this rev into the segments:
            do j = 1, 8
                call me%segs(iseg+j)%set_input(t0(j),tf(j),x0_rotating(:,j))
            end do

            ! for next rev:
            iseg = iseg + 8

        end if

    end do

    end subroutine put_x_in_segments
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the sparsity pattern for the "forward-backward" Halo problem.

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

    integer :: k,ii,jj,icol_start,irow_start
    integer :: n_nonzero  !! number of nonzero elements in the jacobian
    integer :: n  !! number of optimization variables
    integer :: m  !! number of functions
    integer :: iblock
    integer,dimension(:),allocatable :: cols_to_remove !! the columns from the full pattern to remove
    logical,dimension(:),allocatable :: mask
    integer,dimension(:),allocatable :: icol_tmp

    ! get the size of the full problem, we will first construct the full pattern:
    call me%define_problem_size(n=n, m=m, n_nonzero=n_nonzero, full_problem=.true.)

    !-------------------------------------
    ! test... dense pattern
    !-------------------------------------
    ! allocate(irow(n*m))
    ! allocate(icol(n*m))
    ! k = 0
    ! do ii = 1, n
    !     do jj = 1, m
    !         k = k + 1
    !         icol(k) = ii
    !         irow(k) = jj
    !     end do
    ! end do
    ! return
    !-------------------------------------

    allocate(irow(n_nonzero))
    allocate(icol(n_nonzero))

    k = 0
    do iblock = 1, me%n_revs*4 ! block loop

        ! see Figure 22b in Modern Fortran paper
        icol_start = (iblock-1)*7 + 1  ! [1-14, 8-21, 15-29, ...]
        irow_start = (iblock-1)*6 + 1  ! [1-6,  7-12, 13-18, ...]

        do ii = icol_start,icol_start+13
            do jj = irow_start,irow_start+5
                k = k + 1
                irow(k) = jj
                icol(k) = ii
            end do
            ! if necessary, add on the last row for the
            ! initial rdot constraint (depends on initial r,v)
            if (me%constrain_initial_rdot) then
                select case (ii) ! column for initial r,v opt vars
                case(2:7)
                    k = k + 1
                    irow(k) = m  ! last row
                    icol(k) = ii
                end select
            end if
        end do

    end do

    ! the full pattern was generated above.
    ! Here we just remove the columns (optimization variables) we don't need.

    allocate(cols_to_remove(0))
    if (me%fix_initial_time) &
        cols_to_remove = [cols_to_remove, 1]
    if (me%fix_initial_r) &
        cols_to_remove = [cols_to_remove, [2,3,4]]

    if (me%fix_ry_at_end_of_rev > 0) then
        ! remove Ry at end of specified rev
        ! e.g. rev 1 : SEG8 Ry (km) : column 31
        !
        ! patch point: [t rx ry rz vx vy vz]
        !                    ^^
        cols_to_remove = [cols_to_remove, me%fix_ry_at_end_of_rev*7*5 - 4]
    end if

    if (me%fix_final_ry_and_vx) &
        cols_to_remove = [cols_to_remove, [n-4,n-2]] ! last segment Ry & Vx

    if (size(cols_to_remove) > 0) then

        allocate(mask(size(icol)))
        mask = .true. ! keep these
        do k = 1, size(cols_to_remove)
            where(icol == cols_to_remove(k)) mask = .false. ! remove these
        end do

        ! remove the cols:
        irow = pack(irow, mask)
        icol = pack(icol, mask)

        ! decrement ones > the ones removed
        icol_tmp = icol ! so we have the original indices
        do k = 1, size(cols_to_remove)
            where (icol_tmp > cols_to_remove(k)) icol = icol - 1
        end do

        call me%define_problem_size(n_nonzero=n_nonzero)

    end if

    end subroutine get_sparsity_pattern
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the mission.

    subroutine initialize_the_mission(me,x)

    implicit none

    class(mission_type),intent(inout) :: me
    real(wp),dimension(:),allocatable,intent(out),optional :: x !! initial guess

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
    integer :: n_nonzero                        !! number of nonzero elements in the jacobian
    real(wp),dimension(:),allocatable :: dpert  !! perturbation step size array
    real(wp),dimension(:),allocatable :: xlow   !! lower bounds array
    real(wp),dimension(:),allocatable :: xhigh  !! upper bounds
    character(len=10) :: seg_name               !! segment name
    integer,dimension(:),allocatable :: irow    !! sparsity pattern nonzero elements row indices
    integer,dimension(:),allocatable :: icol    !! sparsity pattern nonzero elements column indices
    logical :: use_openmp !! if OpenMP is being used
    real(wp) :: et0 !! initial et for splined ephemeris
    real(wp) :: etf !! final et for splined ephemeirs

    use_openmp = .false.
!$  use_openmp = .true.

    ! problem size:
    call me%define_problem_size(n=n,m=m,n_segs=n_segs,n_nonzero=n_nonzero)

    write(*,*) ''
    write(*,'(A,1X,I10)') '    n:                                     ', n
    write(*,'(A,1X,I10)') '    m:                                     ', m
    write(*,'(A,1X,I10)') '    total number of elements in Jacobian:  ', n*m
    write(*,'(A,1X,I10)') '    number of nonzero elements in Jacobian:', n_nonzero
    write(*,'(A,1X,I10)') '    number of segments:                    ', n_segs
    write(*,'(A,1X,I10)') '    number of revs:                        ', me%n_revs
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
                       info                       = halo_grad_info,&
                       sparsity_mode              = 3,&        ! specified below
                       jacobian_method            = 3,&        ! standard central diff
                       perturb_mode               = 1,&        ! absolute mode
                       partition_sparsity_pattern = .true.,&   ! partition the pattern
                       cache_size                 = 1000 )

    ! generate and set the sparsity pattern for this problem:
    call me%get_sparsity_pattern(irow,icol)
    call me%set_sparsity_pattern(irow,icol)

    if (me%use_splined_ephemeris .or. grav_frame==2) then
        ! for splinting, have to compute et0, etf based on mission range.
        ! use a 2 period buffer.
        et0 = me%et_ref - 2.0*me%period*day2sec
        etf = me%et_ref + (me%period*day2sec * me%n_revs) + 2.0*me%period*day2sec
    end if

    ! set up the ephemeris:
    !write(*,*) 'loading ephemeris file: '//trim(me%ephemeris_file)
    if (me%use_splined_ephemeris) then
        allocate(jpl_ephemeris_splined :: me%eph)
        associate (eph => me%eph)
            select type(eph)
            class is (jpl_ephemeris_splined)
                call eph%initialize_splinded_ephemeris(filename=me%ephemeris_file,&
                                                       status_ok=status_ok,&
                                                       et0=et0,dt=me%dt_spline_sec,etf=etf)
                ! If not using jupiter, destroy it here
                ! (note: inefficient but it's pretty fast so not a big deal)
                if (.not. me%include_pointmass_jupiter) call eph%jupiter_eph_interface%destroy()
            end select
        end associate
        if (.not. status_ok) error stop 'error initializing splined ephemeris'
    else
        allocate(jpl_ephemeris :: me%eph)
        call me%eph%initialize(filename=me%ephemeris_file,status_ok=status_ok)
        if (.not. status_ok) error stop 'error initializing ephemeris'
    end if

    if (.not. me%pointmass_central_body) then
        ! set up the force model [main body is moon]:
        allocate(me%grav)
        call me%grav%initialize(me%gravfile,grav_n,grav_m,status_ok)
        if (.not. status_ok) error stop 'error initializing gravity model'
        if (grav_frame==2) then ! using the splined moon_pa frame for the gravity model
            allocate(me%moon_pa)
            call me%moon_pa%initialize(me%moon_pa_file, et0=et0, etf=etf)
        end if
    end if

    ! now, we set up the segment structure for the problem we are solving:
    ! This is for the "forward-backward" method from the paper (see Figure 2b):
    allocate(me%segs(n_segs))

    ! set up the integrators:
    do i = 1, n_segs

        call me%segs(i)%initialize(n_eoms,maxnum,ballistic_derivs,&
                                   [me%rtol],[me%atol],&
                                   report=trajectory_export_func)

        ! make a copy of this global variable in the segment.
        ! [this was set in read_config_file]
        ! [figure out a way to avoid this..]
        me%segs(i)%data%et_ref = me%et_ref

        if (use_openmp) then
            ! make a copy for each segment, so they can run in parallel
            if (.not. me%pointmass_central_body) then
                allocate(me%segs(i)%grav, source = me%grav)  ! maybe not necessary? (is this threadsafe?)
                if (grav_frame==2) allocate(me%segs(i)%moon_pa, source = me%moon_pa)
            end if
            allocate(me%segs(i)%eph,  source = me%eph)
        else
            ! for serial use, each seg just points to the global ones for the whole mission
            if (.not. me%pointmass_central_body) then
                me%segs(i)%grav => me%grav
                if (grav_frame==2) me%segs(i)%moon_pa  => me%moon_pa
            end if
            me%segs(i)%eph  => me%eph
        end if
        me%segs(i)%pointmass_central_body = me%pointmass_central_body
        me%segs(i)%include_pointmass_earth   =  me%include_pointmass_earth
        me%segs(i)%include_pointmass_sun     =  me%include_pointmass_sun
        me%segs(i)%include_pointmass_jupiter =  me%include_pointmass_jupiter

        ! set 3rd body grav routine:
        if (use_battin_gravity) then
            me%segs(i)%third_body_gravity => third_body_gravity_alt
        else
            me%segs(i)%third_body_gravity => third_body_gravity
        end if

        ! name all the segments:
        write(seg_name,'(I10)') i
        me%segs(i)%name = trim(adjustl(seg_name))

    end do

    ! load the initial guess from the patch point file:
    i = 0
    iseg = 0
    t_periapsis = 0.0_wp
    do irev = 1, me%n_revs

        ! these are all unscaled values:

        t0(1) = t_periapsis + me%periapsis%t

        if (me%fix_initial_r .and. allocated(me%initial_r)) then
            x0_rotating(1:3,1) = me%initial_r(1:3) ! use the user-specified r vector
            x0_rotating(4:6,1) = me%periapsis%rv(4:6)
        else
            x0_rotating(:,1) = me%periapsis%rv
        end if

        tf(1) = t0(1) + me%period8

        t0(2) = t_periapsis + me%quarter%t
        x0_rotating(:,2) = me%quarter%rv
        tf(2) = tf(1)

        t0(4) = t_periapsis + me%apoapsis%t
        x0_rotating(:,4) = me%apoapsis%rv
        tf(4) = t0(4) - me%period8

        t0(3) = t0(2)
        x0_rotating(:,3) = x0_rotating(:,2)
        tf(3) = tf(4)

        t0(5) = t0(4)
        x0_rotating(:,5) = x0_rotating(:,4)
        tf(5) = t0(5) + me%period8

        t0(6) = t_periapsis + me%apoapsis%t + me%quarter%t
        x0_rotating(:,6) = [ me%quarter%rv(1),-me%quarter%rv(2),me%quarter%rv(3),&
                            -me%quarter%rv(4),me%quarter%rv(5),-me%quarter%rv(6)]
        tf(6) = tf(5)

        t0(8) = t_periapsis + me%period
        x0_rotating(:,8) = me%periapsis%rv
        tf(8) = t0(8) - me%period8

        t0(7) = t0(6)
        x0_rotating(:,7) = x0_rotating(:,6)
        tf(7) = tf(8)

        do j = 1, 8

            ! put the data for this rev into the segments:
            call me%segs(iseg+j)%set_input(t0(j),tf(j),x0_rotating(:,j))

            ! also set the scales:
            ! the t0 gradually grows larger so use the initial value for each segment.
            me%segs(iseg+j)%data%t0_scale = abs(magnitude(2.0_wp*me%segs(iseg+j)%data%t0))    ! magic number

            ! do the same for the states, but just in case, specify the min values:
            me%segs(iseg+j)%data%x0_rotating_scale = abs(magnitude(2.0_wp*me%segs(iseg+j)%data%x0_rotating, xscale_x0))

            me%segs(iseg+j)%data%xf_rotating_scale = fscale_xf  ! these are all the same for the constraints

            ! scale the dperts
            me%segs(iseg+j)%data%x0_rotating_dpert = me%segs(iseg+j)%data%x0_rotating_dpert / me%segs(iseg+j)%data%x0_rotating_scale
            me%segs(iseg+j)%data%t0_dpert = me%segs(iseg+j)%data%t0_dpert / me%segs(iseg+j)%data%t0_scale

        end do

        ! for next rev:
        i = i + 35
        iseg = iseg + 8
        t_periapsis = t_periapsis + me%period ! add another period

    end do


    allocate(me%xscale(n))
    allocate(me%xname(n))
    allocate(me%dpert_(n))
    allocate(me%fscale(m))

    ! also set the scale factors:
    call me%get_scales_from_segs()

    ! update with the new dperts:
    call me%set_dpert(me%dpert_)

    if (present(x)) then
        ! get the initial guess from the mission:
        allocate(x(n))    ! size the opt var vector
        call me%get_problem_arrays(x=x)
    end if

    end subroutine initialize_the_mission
!*****************************************************************************************

!*****************************************************************************************
!>
!  Populate the `xscale` and `fscale` problem arrays from the segment data.

    subroutine get_scales_from_segs(me)

    implicit none

    class(mission_type),intent(inout) :: me

    integer :: i,j,ii !! counter
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
    call me%define_problem_size(n_segs=n_segs)

    i = 0 ! for xscale
    ii = 0 ! for dpert
    j = 0 ! for xname
    ! x scales - segment 1:
    if (.not. me%fix_initial_time) then
        call fill_vector(me%xscale, me%segs(1)%data%t0_scale, i)
        call fill_vector(me%dpert_, me%segs(1)%data%t0_dpert, ii)
        call fill_vector(me%xname, 'SEG1 '//t0_label, j)
    end if
    if (me%fix_initial_r) then
        call fill_vector(me%xscale, me%segs(1)%data%x0_rotating_scale(4:6), i)
        call fill_vector(me%dpert_, me%segs(1)%data%x0_rotating_dpert(4:6), ii)
        call fill_vector(me%xname, 'SEG1 '//x0_label(4:6), j)
    else
        call fill_vector(me%xscale, me%segs(1)%data%x0_rotating_scale, i)
        call fill_vector(me%dpert_, me%segs(1)%data%x0_rotating_dpert, ii)
        call fill_vector(me%xname, 'SEG1 '//x0_label, j)
    end if
    ! x scales - the rest:
    do iseg = 2, n_segs, 2
        write(iseg_str,'(I10)') iseg
        call fill_vector(me%xscale, me%segs(iseg)%data%t0_scale, i)
        call fill_vector(me%dpert_, me%segs(iseg)%data%t0_dpert, ii)
        call fill_vector(me%xname, 'SEG'//trim(adjustl(iseg_str))//' '//t0_label, j)
        if (iseg == me%fix_ry_at_end_of_rev*8) then
            call fill_vector(me%xscale, me%segs(iseg)%data%x0_rotating_scale([1,3,4,5,6]), i)
            call fill_vector(me%dpert_, me%segs(iseg)%data%x0_rotating_dpert([1,3,4,5,6]), ii)
            call fill_vector(me%xname, 'SEG'//trim(adjustl(iseg_str))//' '//x0_label([1,3,4,5,6]), j)
            cycle
        else if (me%fix_final_ry_and_vx .and. iseg == n_segs) then ! last state point
            call fill_vector(me%xscale, me%segs(iseg)%data%x0_rotating_scale([1,3,5,6]), i)
            call fill_vector(me%dpert_, me%segs(iseg)%data%x0_rotating_dpert([1,3,5,6]), ii)
            call fill_vector(me%xname, 'SEG'//trim(adjustl(iseg_str))//' '//x0_label([1,3,5,6]), j)
            cycle
        else
            ! otherwise, full state:
            call fill_vector(me%xscale, me%segs(iseg)%data%x0_rotating_scale, i)
            call fill_vector(me%dpert_, me%segs(iseg)%data%x0_rotating_dpert, ii)
            call fill_vector(me%xname, 'SEG'//trim(adjustl(iseg_str))//' '//x0_label, j)
        end if
    end do

    ! f scales:
    i = 0 ! reset for f
    do iseg = 2, n_segs, 2
        call fill_vector(me%fscale, me%segs(iseg)%data%xf_rotating_scale, i)
    end do
    if (me%constrain_initial_rdot) then
        call fill_vector(me%fscale, fscale_rdot, i)  ! scale rdot constraint
    end if
    ! write(*,*) ''
    ! write(*,*) 'xscale: ', me%xscale
    ! write(*,*) ''
    ! write(*,*) 'fscale: ', me%fscale
    ! write(*,*) ''
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
    reaL(wp),dimension(3) :: r_earth_wrt_moon,r_sun_wrt_moon,r_jupiter_wrt_moon
    real(wp),dimension(3,3) :: rotmat
    real(wp),dimension(3) :: a_geopot  !! central body acc
    real(wp),dimension(3) :: a_earth
    real(wp),dimension(3) :: a_sun
    real(wp),dimension(3) :: a_jupiter
    real(wp),dimension(3) :: a_third_body  !! total third-body acc
    real(wp) :: et !! ephemeris time of `t`
    logical :: status_ok

    select type (me)

    class is (segment)

        ! get state:
        r = x(1:3)
        v = x(4:6)

        ! compute ephemeris time [sec]:
        et = me%data%et_ref + t

        if (me%pointmass_central_body) then
            ! pointmass moon gravity model
            a_geopot = -mu_moon / norm2(r)**3 * r
        else
            ! geopotential gravity:
            ! first get the rotation matrix from inertial to body-fixed Moon frame
            select case (grav_frame)
            case (1); rotmat = icrf_to_iau_moon(et) ! iau_moon
            case (2); rotmat = me%moon_pa%j2000_to_frame(et) ! moon_pa (splined)
            case default; error stop 'invalid grav_frame in ballistic_derivs'
            end select
            rb = matmul(rotmat,r)           ! r in body-fixed frame
            call me%grav%get_acc(rb,grav_n,grav_m,a_geopot)  ! get the acc due to the geopotential
            a_geopot = matmul(transpose(rotmat),a_geopot)    ! convert acc back to inertial frame
        end if

        ! third-body state vectors (wrt the central body, which is the moon in this case):
        ! [inertial frame]
        a_third_body  = zero
        if (me%include_pointmass_earth) then
            call me%eph%get_r(et,body_earth,body_moon,r_earth_wrt_moon,status_ok)
            call me%third_body_gravity(r,r_earth_wrt_moon,mu_earth,a_earth)
            a_third_body = a_third_body + a_earth
        end if
        if (me%include_pointmass_sun) then
            call me%eph%get_r(et,body_sun,body_moon,r_sun_wrt_moon,status_ok)
            call me%third_body_gravity(r,r_sun_wrt_moon,mu_sun,a_sun)
            a_third_body = a_third_body + a_sun
        end if
        if (me%include_pointmass_jupiter) then
            call me%eph%get_r(et,body_jupiter,body_moon,r_jupiter_wrt_moon,status_ok)
            call me%third_body_gravity(r,r_jupiter_wrt_moon,mu_jupiter, a_jupiter)
            a_third_body = a_third_body + a_jupiter
        end if

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
    et0 = me%data%et_ref + t0 * day2sec  ! convert to ephemeris time [sec]
    rotating = two_body_rotating_frame(primary_body=body_earth,&
                                        secondary_body=body_moon,&
                                        center=center_at_secondary_body,&
                                        et=et0)
    inertial = icrf_frame(b=body_moon)

    ! from rotating to inertial:
    call rotating%transform(x0_rotating,inertial,et0,me%eph,x0_inertial,status_ok)
    if (.not. status_ok) then
        write(*,*) ''
        write(*,'(A/,*(E30.16/))') 't0:         ', t0
        write(*,'(A/,*(E30.16/))') 'et0:        ', et0
        write(*,'(A/,*(E30.16/))') 'x0_rotating:', x0_rotating
        write(*,*) ''
        error stop 'transformation error in set_segment_inputs'
    end if

    ! set the inputs (needed by the propagator)
    me%data%t0          = t0
    me%data%x0_rotating = x0_rotating
    me%data%x0          = x0_inertial
    me%data%tf          = tf

    end subroutine set_segment_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Gets the initial states of a segment

    subroutine get_segment_inputs(me,t0,x0_rotating)

    implicit none

    class(segment),intent(in) :: me
    real(wp),intent(out),optional :: t0
    real(wp),dimension(6),intent(out),optional :: x0_rotating !! rotating frame

    if (present(t0)) t0 = me%data%t0
    if (present(x0_rotating)) x0_rotating = me%data%x0_rotating

    end subroutine get_segment_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  After propagating a segment, this gets the outputs.

    subroutine get_segment_outputs(me,xf,xf_rotating, x0_rotating)

    implicit none

    class(segment),intent(in) :: me
    real(wp),dimension(6),intent(out),optional :: xf             !!  inertial frame
    real(wp),dimension(6),intent(out),optional :: xf_rotating    !!  rotating frame
    real(wp),dimension(6),intent(out),optional :: x0_rotating

    if (present(xf))          xf = me%data%xf
    if (present(xf_rotating)) xf_rotating = me%data%xf_rotating
    if (present(x0_rotating)) x0_rotating = me%data%x0_rotating

    end subroutine get_segment_outputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the outputs of a segment, assuming it has been propagated elsewhere

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

    subroutine propagate_segment(me,mode,tstep)

    implicit none

    class(segment),intent(inout) :: me
    integer,intent(in),optional :: mode  !! 1 - don't report steps, 2 - report steps (for plotting)
    real(wp),intent(in),optional :: tstep !! fixed time step for mode=2

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
    if (present(tstep)) then
        call me%integrate(t,x,tf,idid=idid,integration_mode=integration_mode,tstep=tstep)
    else
        call me%integrate(t,x,tf,idid=idid,integration_mode=integration_mode)
    end if
    if (idid<0) then
        write(*,'(A,*(I5/))')    'idid: ',idid
        error stop 'error in integrator'
    end if

    xf = x  ! final state [inertial frame]

    ! also save the rotating frame state at tf (to compute the constraint violations):
    etf = me%data%et_ref + tf  ! convert to ephemeris time [sec]
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
    integer :: m    !! number of constraints

    ! In this example, there are no dependencies, so we assume they
    ! can be propagated in any order.

    if (present(funcs_to_compute)) then

        ! this is for the full "forward-backward" problem formulation:

        !                                               optional:
        !       1-6      7-12    13-19     20-26        [m]
        ! f = [xf1-xf2, xf3-xf4, xf5-xf6, xf7-xf8, ..., [rdot] ]
        !
        !      0     1     2
        ! f = [123456123456123456]
        !
        ! 0: 1-2
        ! 1: 3-4
        ! 2: 5-6
        ! 3: 7-8

        if (me%constrain_initial_rdot) call me%define_problem_size(m=m) ! we need m below

        do i = 1, size(funcs_to_compute)
            if (me%constrain_initial_rdot) then
                if (funcs_to_compute(i)==m) cycle  ! for the rdot constraint, handle below separately
            end if
            j = (funcs_to_compute(i)-1) / 6 ! 0,1,2,...
            call add_it(j*2+1,isegs_to_propagate) ! ... is this correct ?
            call add_it(j*2+2,isegs_to_propagate)
        end do

        if (me%constrain_initial_rdot) then  ! for the rdot constraint we have to propagate segment 1
            if (any(m==funcs_to_compute)) then
                call add_it(1,isegs_to_propagate)
            end if
        end if

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

    integer :: i  !! counters
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
!$OMP PARALLEL DO    !...FIRSTPRIVATE(me)
    do i = 1, size(isegs)
        call me%segs(isegs(i))%propagate()
    end do
!$OMP END PARALLEL DO
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

    subroutine halo_func(me,x,f)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:),intent(out)    :: f

    select type (me)
    class is (my_solver_type)
        call me%mission%constraint_violations(x,f)
    class default
        error stop 'invalid class in halo_func'
    end select

    end subroutine halo_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the gradient of the solver function (Jacobian matrix).
!  Dense version.

    subroutine halo_grad(me,x,g)

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
        error stop 'invalid class in halo_grad'
    end select

    end subroutine halo_grad
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the gradient of the solver function (Jacobian matrix).
!  Sparse version.

    subroutine halo_grad_sparse(me,x,g)

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(out)   :: g

    real(wp),dimension(:),allocatable :: jac !! the jacobian matrix returned by `numdiff`
        ! ...note: need to modify so it doesn't
        !          have to return an allocatable array

    integer :: i  !! seg number counter

    select type (me)
    class is (my_solver_type)

        ! first let's cache all the segment data:
        do i=1,size(me%mission%segs)
            call me%mission%segs(i)%cache()
        end do

        ! use numdiff to compute the jacobian matrix (sparse version)
        call me%mission%compute_jacobian(x,jac)
        g = jac

        ! restore data just in case:
        do i=1,size(me%mission%segs)
            call me%mission%segs(i)%uncache()
        end do

    class default
        error stop 'invalid class in halo_grad'
    end select

    end subroutine halo_grad_sparse
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

    subroutine halo_grad_info(me,column,i,x)

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
        error stop 'invalid class in halo_grad_info'
    end select

    end subroutine halo_grad_info
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
        error stop 'invalid class in halo_grad'
    end select

    end subroutine my_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  export an iteration from the solver

    subroutine halo_export(me,x,f,iter)

    use iso_fortran_env, ip => int64

    implicit none

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)     :: x
    real(wp),dimension(:),intent(in)     :: f
    integer,intent(in)                   :: iter !! iteration number

    integer(ip),save :: count_start, count_end, count_rate

    call system_clock(count = count_end, count_rate=count_rate)

    if (iter==1) then
        write(*,'(A4,1X,*(A30,1X))') 'ITER', 'NORM(X)', 'NORM(F)', 'DT(sec)'
    end if

    if (iter==1) then
        write(*,'(I4,1X,2(F30.16,1X),a1)') iter, norm2(x), norm2(f), ''
    else
        write(*,'(I4,1X,2(F30.16,1X),*(I30,1x))') iter, norm2(x), norm2(f), &
                                                  int((count_end-count_start)/real(count_rate,wp), ip)

    end if
    !write(*,'(I4,1X,*(F30.16,1X))') iter, x

    count_start = count_end

    !select type (me)
    !class is (my_solver_type)
    !    write(*,*) 'rdot = ', compute_rdot(me%mission%segs(1)%data%x0)
    !end select

    end subroutine halo_export
!*****************************************************************************************

!*****************************************************************************************
!>
!  Plot the trajectory using matplotlib and/or generate a text file (for MKSPK)
!
!@note It is assumed that all the data is present in the segments needed to propagate.

    subroutine plot_trajectory(me,filename,export_trajectory,draw_trajectory,only_first_rev)

    implicit none

    class(mission_type),intent(inout) :: me
    character(len=*),intent(in) :: filename !! plot file name [without extension]
    logical,intent(in),optional :: export_trajectory    !! if true, a text trajectory
                                                        !! file is also produced
                                                        !! (that can be read by MKSPK)
                                                        !! [default is False]
    logical,intent(in),optional :: draw_trajectory      !! if true [Default], a matplotlib
                                                        !! plot is generated
    logical,intent(in),optional :: only_first_rev       !! to only do the first rev [Default is False]

    type(pyplot) :: plt  !! for generating the plots
    integer :: iseg  !! segment number counter
    integer :: istat !! pyplot status code
    character(len=10) :: iseg_str !! string version of segment number
    logical :: export  !! if the txt file is to be produced
    logical :: plot    !! if the py file is to be produced
    integer :: iunit !! file unit for the txt file
    integer :: i !! counter for txt file write
    integer :: istart,iend,istep  !! index counters
    integer :: nsegs_to_plot !! number of segments to plot
    !integer :: iendprev
    real(wp) :: last_et_written !! to keep track of the last et written to the file
    logical :: first !! the first point has been written to the file

    integer,dimension(2),parameter :: figsize = [20,20] !! figure size for plotting
    logical,parameter :: plot_rotating = .true.  !! if true, the rotating state is plotted.
                                                 !! if false, the inertial state is plotted.

    ! optional arguments:
    if (present(export_trajectory)) then
        export = export_trajectory
    else
        export = .false.
    end if
    if (present(draw_trajectory)) then
        plot = draw_trajectory
    else
        plot = .true.
    end if
    nsegs_to_plot = size(me%segs) ! default export all the segments
    if (present(only_first_rev)) then
        if (only_first_rev) nsegs_to_plot = 8 ! only the first rev (8 segments)
    end if

    if (.not. plot .and. .not. export) return ! nothing to do

    if (plot) then
        ! initialize the plot:
        call plt%initialize(grid=.true.,&
                            xlabel='\n\n\nx (km)',&
                            ylabel='\n\n\ny (km)',&
                            zlabel='\n\n\nz (km)',&
                            figsize         =figsize,&
                            font_size       = 20,&
                            axes_labelsize  = 25,&
                            xtick_labelsize = 20,&
                            ytick_labelsize = 20,&
                            ztick_labelsize = 20,&
                            title='Halo Trajectory: '//trim(filename),&
                            legend=.false.,&
                            axis_equal=.true.,&
                            mplot3d=.true.)
    end if

    if (export) then
        open(newunit=iunit,file=trim(filename)//'_'//me%get_case_name()//&
                                '.txt',status='REPLACE',iostat=istat)
        if (istat/=0) error stop 'error opening trajectory file.'
    end if

    do iseg = 1, size(me%segs)
        ! destroy all trajectories:
        call me%segs(iseg)%traj_inertial%destroy()
        call me%segs(iseg)%traj_rotating%destroy()
        call me%segs(iseg)%traj_se_rotating%destroy()
    end do

    first = .false.  ! to flag when the first point has been written
    do iseg = 1, nsegs_to_plot

        write(iseg_str,'(I10)') iseg

        ! generate the trajectory for this segment:
        call me%segs(iseg)%propagate(mode=2)    ! [export points]

        if (plot) then
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
        end if

        if (export) then

            ! write the trajectory data:     ! Note: THIS HAS TO BE THE INERTIAL FRAME FOR THE SPK FILE

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

            ! note: MKSPK does not allow duplicate time tags and we also want to avoid points
            ! too close together since that may cause interpolation artifacts. so we have the
            ! user-specified min_export_time_step. we don't allow points to be sent to the
            ! bsp file where the time delta is <= this value.

            do i=istart,iend,istep
                associate (et => me%segs(iseg)%traj_inertial%et(i))
                    if (.not. first) then
                        if ((et - last_et_written) <= min_export_time_step) cycle
                    else
                        first = .true.
                    end if
                    last_et_written = et ! save this so we can track this delta
                    ! inertial trajectory for SPK:
                    write(iunit,'(*(E30.16E3,A,1X))',iostat=istat) &
                        et,';',&
                        me%segs(iseg)%traj_inertial%x(i) ,';',&
                        me%segs(iseg)%traj_inertial%y(i) ,';',&
                        me%segs(iseg)%traj_inertial%z(i) ,';',&
                        me%segs(iseg)%traj_inertial%vx(i),';',&
                        me%segs(iseg)%traj_inertial%vy(i),';',&
                        me%segs(iseg)%traj_inertial%vz(i),''
                end associate
            end do
        end if

    end do

    if (plot) then
        ! add the moon as a sphere:
        call plt%add_sphere(r=rad_moon,xc=zero,yc=zero,zc=zero,istat=istat,color='Grey')

        ! save the plot:
        ! call plt%showfig(pyfile=trim(filename)//'.py',istat=istat)
        call plt%savefig(trim(filename)//'_'//me%get_case_name()//'.png',&
                        pyfile=trim(filename)//'_'//me%get_case_name()//'.py',dpi='300',&
                        istat=istat)

        ! cleanup:
        call plt%destroy()
    end if

    if (export) close(iunit,iostat=istat)

    end subroutine plot_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Export the trajectory JSON file.
!
!@note It is assumed that all the data is present in the segments needed to propagate.

    subroutine export_trajectory_json_file(me,filename,only_first_rev)

    implicit none

    class(mission_type),intent(inout) :: me
    character(len=*),intent(in) :: filename !! plot file name [without extension]
    logical,intent(in),optional :: only_first_rev  !! to only do the first rev [Default is False]

    integer :: iseg  !! segment number counter
    integer :: nsegs_to_plot !! number of segments to plot
    type(json_core) :: json
    type(json_value),pointer :: p_root, p_segs, p_seg, p_current
    real(wp),dimension(:),allocatable :: rdot, rmag
    logical :: accumulate_rdot !! to generate the rdot/rp/ra file
    real(wp),dimension(:),allocatable :: cumulative_et,cumulative_rdot,cumulative_r
    integer :: n

    write(*,'(A)') ' * Export the trajectory JSON file.'

    ! optional arguments:
    nsegs_to_plot = size(me%segs) ! default export all the segments
    if (present(only_first_rev)) then
        if (only_first_rev) nsegs_to_plot = 8 ! only the first rev (8 segments)
    end if
    accumulate_rdot = me%generate_rp_ra_file

    ! the JSON file will contain an array of segments:
    call json%initialize(compress_vectors=.true.)
    call json%create_object(p_root, '')
    call json%create_array(p_segs, 'segs')
    call json%add(p_root, p_segs)
    p_current => null()

    do iseg = 1, nsegs_to_plot
        call destroy_traj(iseg)
    end do

    !====================================
    ! now propagate the segments:
!$OMP PARALLEL DO    !...FIRSTPRIVATE(me)
    do iseg = 1, nsegs_to_plot
        call me%segs(iseg)%propagate(mode=2)  ! [export points]
    end do
!$OMP END PARALLEL DO
    !====================================

    if (accumulate_rdot) then
        allocate(cumulative_et(0))
        allocate(cumulative_rdot(0))
        allocate(cumulative_r(0))
    end if

    do iseg = 1, nsegs_to_plot

        associate (seg => me%segs(iseg))
            ! generate the trajectory for this segment:
            !call destroy_traj(iseg)
            !call seg%propagate(mode=2)  ! [export points]

            ! create the segment object for exporting the trajectory:
            call json%create_object(p_seg, '')
            if (associated(p_current)) then
                call json%insert_after(p_current, p_seg) ! next one
            else
                call json%add(p_segs, p_seg) ! first one
            end if
            p_current => p_seg ! update for next seg

            call json%add(p_seg, 'iseg', iseg)
            call json%add(p_seg, 'et', seg%traj_inertial%et)

            call json%add(p_seg, 'x_inertial',  seg%traj_inertial%x)
            call json%add(p_seg, 'y_inertial',  seg%traj_inertial%y)
            call json%add(p_seg, 'z_inertial',  seg%traj_inertial%z)
            call json%add(p_seg, 'vx_inertial', seg%traj_inertial%vx)
            call json%add(p_seg, 'vy_inertial', seg%traj_inertial%vy)
            call json%add(p_seg, 'vz_inertial', seg%traj_inertial%vz)

            if (accumulate_rdot) then
                call compute_rdot_vecs(seg%traj_inertial%x,seg%traj_inertial%y,seg%traj_inertial%z,&
                                       seg%traj_inertial%vx,seg%traj_inertial%vy,seg%traj_inertial%vz,&
                                       rmag, rdot)
                call append_traj_to_arrays(seg)
                !call json%add(p_seg, 'rdot_inertial', rdot)
                ! don't include the last time step since that will overlap the next segment
                ! n = size(seg%traj_inertial%x)
                ! cumulative_et   = [cumulative_et,   seg%traj_inertial%et]
                ! cumulative_rdot = [cumulative_rdot, rdot]
                ! cumulative_r    = [cumulative_r,    rmag]
            end if

            call json%add(p_seg, 'x_rotating',  seg%traj_rotating%x)
            call json%add(p_seg, 'y_rotating',  seg%traj_rotating%y)
            call json%add(p_seg, 'z_rotating',  seg%traj_rotating%z)
            ! call json%add(p_seg, 'vx_rotating', seg%traj_rotating%vx)
            ! call json%add(p_seg, 'vy_rotating', seg%traj_rotating%vy)
            ! call json%add(p_seg, 'vz_rotating', seg%traj_rotating%vz)
            ! call compute_rdot_vecs(seg%traj_rotating%x,seg%traj_rotating%y,seg%traj_rotating%z,&
            !                        seg%traj_rotating%vx,seg%traj_rotating%vy,seg%traj_rotating%vz,&
            !                        rdot)
            ! call json%add(p_seg, 'rdot_rotating', rdot)

            call json%add(p_seg, 'x_se_rotating',  seg%traj_se_rotating%x)
            call json%add(p_seg, 'y_se_rotating',  seg%traj_se_rotating%y)
            call json%add(p_seg, 'z_se_rotating',  seg%traj_se_rotating%z)
            ! call json%add(p_seg, 'vx_se_rotating', seg%traj_se_rotating%vx) ! don't need these
            ! call json%add(p_seg, 'vy_se_rotating', seg%traj_se_rotating%vy)
            ! call json%add(p_seg, 'vz_se_rotating', seg%traj_se_rotating%vz)

            !TODO:
            !  - maybe also earth & sun ephemeris for plotting
            !  - solar fraction to color the trajectory [but really that should go in the eclipse file?]

            !call destroy_traj(iseg) ! keep them for the eclipse file generation ...

        end associate
    end do

    if (accumulate_rdot) then
        call me%export_rp_ra_json_file(cumulative_et, cumulative_r, cumulative_rdot, filename)
    end if

    do iseg = 1, nsegs_to_plot
        call destroy_traj(iseg)
    end do

    call json%print(p_root, trim(filename)//'.json')
    call json%destroy(p_root)

    contains

        subroutine destroy_traj(iseg)
            integer,intent(in) :: iseg !! segment number
            call me%segs(iseg)%traj_inertial%destroy()
            call me%segs(iseg)%traj_rotating%destroy()
            call me%segs(iseg)%traj_se_rotating%destroy()
        end subroutine destroy_traj

        subroutine append_traj_to_arrays(seg)
            type(segment),intent(in) :: seg

            integer :: n, i

            n = size(seg%traj_inertial%et)
            if (seg%traj_inertial%et(1)<=seg%traj_inertial%et(2)) then
                ! forward propagated segment
                cumulative_et   = [cumulative_et,   seg%traj_inertial%et(1:n-1)]
                cumulative_rdot = [cumulative_rdot, rdot(1:n-1)]
                cumulative_r    = [cumulative_r,    rmag(1:n-1)]
            else
                ! some of the segments are propagated backwards, so we have to reverse them
                ! [note this is very inefficient... need to fix this]  -TODO
                cumulative_et   = [cumulative_et,   seg%traj_inertial%et(n:2:-1)]
                cumulative_rdot = [cumulative_rdot, rdot(n:2:-1)]
                cumulative_r    = [cumulative_r,    rmag(n:2:-1)]
            end if

        end subroutine append_traj_to_arrays

    end subroutine export_trajectory_json_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Export the rp/ra file.
!
!@note It is assumed that all the data is present in the segments needed to propagate.

    subroutine export_rp_ra_json_file(me,et,rmag,rdot,filename)

    use root_module
    use bspline_module

    class(mission_type),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: et
    real(wp),dimension(:),intent(in) :: rmag
    real(wp),dimension(:),intent(in) :: rdot
    character(len=*),intent(in) :: filename !! file name [without extension]

    type,extends(chandrupatla_solver) :: my_solver
        !! the solver to use for the root finding
    end type my_solver

    real(wp) :: et0, etf, et_root, rdot_root, initial_et, final_et, rdot0, rdotf, rmag_root
    integer :: i, n
    real(wp),dimension(:),allocatable :: et_for_rp, et_for_ra, rp_vec, ra_vec
    real(wp) :: x, f
    integer :: iflag
    type(my_solver) :: solver
    type(bspline_1d) :: rdot_spline, rmag_spline
    type(json_file) :: json
    logical :: finished !! to indicate the arrays are done
    integer :: n_rp, n_ra, n_rp_et, n_ra_et !! counters for array sizes

    integer,parameter :: kx = 4  !! spline order (cubic bspline)
    integer,parameter :: idx = 0 !! interpolate value only
    real(wp),parameter :: et_step = 3600.0_wp * 24.0_wp !! et step size [should be an input] (sec)
    integer,parameter :: chunk_size = 1000 !! chunk size for expanding the arrays

    write(*,'(A)') ' * Generating ra/rp file.'

    ! error check to make sure et is strictly increasing.
    do i = 2, size(et)
        if (et(i)<=et(i-1)) then
            write(*,*) 'error in et point ', i
            write(*,*) 'et(i)  ', et(i)
            write(*,*) 'et(i-1)', et(i-1)
            error stop
        end if
    end do
    ! first spline rmag and rdot as a functions of et:
    call rdot_spline%initialize(et,rdot,kx,iflag,.false.); call check_spline('rdot')
    call rmag_spline%initialize(et,rmag,kx,iflag,.false.); call check_spline('rmag')

    allocate(et_for_rp(0), et_for_ra(0), rp_vec(0), ra_vec(0)) ! arrays to hold results

    ! step through the trajectory and find all the rdot roots (periapsis and apoapsis):
    call solver%initialize(rdot_func)
    n = size(et)
    initial_et = et(1)
    final_et = et(n)
    et0 = initial_et
    n_rp = 0; n_ra = 0; n_rp_et = 0; n_ra_et = 0
    do
        etf = min(final_et, et0+et_step)
        ! if there is a root on this interval (change of sign of rdot):
        rdot0 = rdot_func(solver, et0)
        rdotf = rdot_func(solver, etf)
        finished = etf==et(n) ! the last step
        if (rdot0*rdotf<=zero) then  ! there is a root on this interval
            ! find the root:
            call solver%solve(et0,etf,et_root,rdot_root,iflag)
            rmag_root = rmag_func(et_root)  ! get corresponding rmag at the root
            if (rdotf>rdot0) then
                ! increasing - periapsis
                ! et_for_rp = [et_for_rp, et_root]
                ! rp_vec    = [rp_vec, rmag_root]
                call add_to(et_for_rp, et_root, n_rp_et, chunk_size, finished)
                call add_to(rp_vec, rmag_root, n_rp, chunk_size, finished)
            else
                ! decreasing - apoapsis
                !et_for_ra = [et_for_ra, et_root]
                !ra_vec    = [ra_vec, rmag_root]
                call add_to(et_for_ra, et_root, n_ra_et, chunk_size, finished)
                call add_to(ra_vec, rmag_root, n_ra, chunk_size, finished)
            end if
        end if
        if (finished) exit
        et0 = etf ! set up for next step
    end do

    ! output results to file:
    call json%initialize(compress_vectors=.true.)
    call json%add('et_for_rp', et_for_rp)
    call json%add('rp_vec',    rp_vec)
    call json%add('et_for_ra', et_for_ra)
    call json%add('ra_vec',    ra_vec)
    call json%print(trim(filename)//'_rp_ra.json')
    call json%destroy()

    contains

        pure subroutine add_to(vec, val, n, chunk_size, finished)
            !! resize an array in chunks

            real(wp), dimension(:), allocatable, intent(inout) :: vec
                !! the vector to add to
            real(wp), intent(in) :: val
                !! the value to add
            integer, intent(inout) :: n
                !! counter for last element added to vec.
                !! must be initialized to size(vec)
                !! (or 0 if not allocated) before first call
            integer, intent(in) :: chunk_size
                !! allocate vec in blocks of this size (>0)
            logical, intent(in) :: finished
                !! set to true to return vec
                !! as its correct size (n)

            real(wp), dimension(:), allocatable :: tmp

            if (allocated(vec)) then
                if (n == size(vec)) then
                    ! have to add another chunk:
                    allocate (tmp(size(vec) + chunk_size))
                    tmp(1:size(vec)) = vec
                    call move_alloc(tmp, vec)
                end if
                n = n + 1
            else
                ! the first element:
                allocate (vec(chunk_size))
                n = 1
            end if

            vec(n) = val

            if (finished) then
                ! set vec to actual size (n):
                if (allocated(tmp)) deallocate (tmp)
                allocate (tmp(n))
                tmp = vec(1:n)
                call move_alloc(tmp, vec)
            end if

        end subroutine add_to

        subroutine check_spline(case) !! error checking for splines
            character(len=*),intent(in) :: case !! case name
            if (iflag/=0) then
                select case (iflag)
                case(2); error stop case//' : bspline error: iknot out of range.'
                case(3); error stop case//' : bspline error: nx out of range.'
                case(4); error stop case//' : bspline error: kx out of range.'
                case(5); error stop case//' : bspline error: x not strictly increasing.'
                end select
            end if
        end subroutine check_spline

        function rdot_func(me,x)  !! compute rdot from spline
            class(root_solver),intent(inout) :: me
            real(wp),intent(in) :: x
            real(wp) :: rdot_func
            integer :: iflag
            call rdot_spline%evaluate(x,idx,rdot_func,iflag)
            if (iflag/=0) error stop 'error evaluating bspline for rdot.'
        end function rdot_func
        function rmag_func(x)  !! compute rmag from spline
            real(wp),intent(in) :: x
            real(wp) :: rmag_func
            integer :: iflag
            call rmag_spline%evaluate(x,idx,rmag_func,iflag)
            if (iflag/=0) error stop 'error evaluating bspline for rmag.'
        end function rmag_func

    end subroutine export_rp_ra_json_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  compute rdot, which is d(rmag)/dt.
!
!  See also: [[compute_rdot_vecs]]

    function compute_rdot(rv) result(rdot)
        real(wp),dimension(6),intent(in) :: rv !! rv vector (km, km/s)
        real(wp) :: rdot !! [km/s]
        associate(r => rv(1:3), v => rv(4:6))
            rdot = dot_product(r,v) / norm2(r)
        end associate
    end function compute_rdot
!*****************************************************************************************

!*****************************************************************************************
!>
!  compute the rmag and rdot vectors, given the state vectors

    subroutine compute_rdot_vecs(x,y,z,vx,vy,vz,rmag,rdot)
        real(wp),dimension(:),intent(in) :: x !! x-position component
        real(wp),dimension(:),intent(in) :: y !! y-position component
        real(wp),dimension(:),intent(in) :: z !! z-position component
        real(wp),dimension(:),intent(in) :: vx !! x-velocity component
        real(wp),dimension(:),intent(in) :: vy !! y-velocity component
        real(wp),dimension(:),intent(in) :: vz !! z-velocity component
        real(wp),dimension(:),intent(out),allocatable :: rmag !! magnitude of position vector [km]
        real(wp),dimension(:),intent(out),allocatable :: rdot !! derivative of radius magniutde [km/s]
        integer :: i !! counter
        real(wp),dimension(3) :: r, v
        allocate(rdot(size(x)),rmag(size(x)))
        do i = 1, size(x)
            r = [x(i), y(i), z(i)]
            v = [vx(i), vy(i), vz(i)]
            rmag(i) = norm2(r)
            rdot(i) = dot_product(r,v) / rmag(i)
        end do
    end subroutine compute_rdot_vecs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Propagate all the segments with a dense time step and export the
!  Eclipsing data w.r.t. Earth. In this file, any point <0 is in eclipse.
!
!  Based on [[plot_trajectory]]
!
!  TODO: use the data from [[export_trajectory_json_file]], so that should be run first

    subroutine generate_eclipse_data(me,fileprefix,filetype)

        class(mission_type),intent(inout) :: me
        character(len=*),intent(in) :: fileprefix !! file prefix for the csv file (case name will be added)
        integer,intent(in),optional :: filetype !! type of output file: 1=csv [default] or 2=json

        integer :: i,iunit,iseg,istat,istart,istep,iend
        real(wp) :: phi !! sunfrac value
        real(wp) :: p
        integer :: ifiletype
        type(json_file) :: json
        character(len=10) :: iseg_str
        character(len=:),allocatable :: full_filename
        real(wp),dimension(:),allocatable :: phi_vec !! vector of `phi` values from a segment (for JSON file)
        real(wp),dimension(:),allocatable :: et_vec !! vector of `et` values from a segment (for JSON file)
        real(wp),dimension(:),allocatable :: phi_vec_total, et_vec_total !! cumulative for all segments

        real(wp),dimension(:),allocatable :: x_se, y_se, z_se, vx_se, vy_se, vz_se, et_se, phi_se
        type(icrf_frame) :: inertial
        type(two_body_rotating_frame) :: se_rotating
        real(wp),dimension(6) :: x_se_rotating
        logical :: status_ok

        integer,parameter :: FILETYPE_CSV = 1
        integer,parameter :: FILETYPE_JSON = 2
        character(len=*),dimension(*),parameter :: file_ext = ['.csv ', &
                                                               '.json']

        if (present(filetype)) then
            ifiletype = filetype
        else
            ifiletype = FILETYPE_JSON
        end if
        full_filename = trim(fileprefix)//'_'//me%get_case_name()//&
                        trim(file_ext(ifiletype))

        write(*,'(A)') ' * Generating eclipse file: '//full_filename

        select case (ifiletype)
        case (FILETYPE_CSV)
            open(newunit=iunit,file=full_filename,status='REPLACE',iostat=istat)
            if (istat/=0) error stop 'error opening eclipse csv file.'
            write(iunit,'(A30,A,1X,A30)',iostat=istat) 'ET (sec)', ',', 'PHI'
        case (FILETYPE_JSON)
            call json%initialize(compress_vectors=.true.)
        case default
            error stop 'invalid filetype for eclipse file'
        end select

        allocate(x_se(0));allocate(y_se(0));allocate(z_se(0))
        allocate(vx_se(0));allocate(vy_se(0));allocate(vz_se(0))
        allocate(et_se(0))
        allocate(phi_se(0))
        allocate(phi_vec_total(0))
        allocate(et_vec_total(0))

        ! destroy all trajectories first:
        do iseg = 1, size(me%segs)
            call me%segs(iseg)%traj_inertial%destroy()
            call me%segs(iseg)%traj_rotating%destroy()
            call me%segs(iseg)%traj_se_rotating%destroy()
        end do

        !====================================
        ! now propagate the segments:
        ! [export points at a fixed time step]
!$OMP PARALLEL DO    !...FIRSTPRIVATE(me)
        do iseg = 1, size(me%segs)
            call me%segs(iseg)%propagate(mode=2,tstep=me%eclipse_dt_step)
        end do
!$OMP END PARALLEL DO
        !====================================

        ! TODO: could also try to do this stuff in parallel too...
        do iseg = 1, size(me%segs)

            if (ifiletype==FILETYPE_JSON) then
                if (allocated(phi_vec)) deallocate(phi_vec); allocate(phi_vec(0))
                if (allocated(et_vec))  deallocate(et_vec);  allocate(et_vec(0))
            end if
            write(iseg_str,'(I10)') iseg    ! segment number as string
            iseg_str = adjustl(iseg_str)

            ! destroy all trajectories first:
            ! call me%segs(iseg)%traj_inertial%destroy()
            ! call me%segs(iseg)%traj_rotating%destroy()
            ! call me%segs(iseg)%traj_se_rotating%destroy()
            ! generate the trajectory for this segment:
            ! [export points at a fixed time step]
            !call me%segs(iseg)%propagate(mode=2,tstep=me%eclipse_dt_step)

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

            do i=istart,iend,istep
                associate (et => me%segs(iseg)%traj_inertial%et(i),&
                            rv_moon => [me%segs(iseg)%traj_inertial%x(i),&
                                        me%segs(iseg)%traj_inertial%y(i),&
                                        me%segs(iseg)%traj_inertial%z(i),&
                                        me%segs(iseg)%traj_inertial%vx(i),&
                                        me%segs(iseg)%traj_inertial%vy(i),&
                                        me%segs(iseg)%traj_inertial%vz(i) ])

                    p = get_sun_fraction(me%eph, et, rv_moon, me%r_eclipse_bubble)
                    phi = min(0.0_wp, p) ! phi >0 mean the spacecraft is in sunlight

                    select case (ifiletype)
                    case (FILETYPE_CSV)
                        write(iunit,'(*(E30.16E3,A,1X))',iostat=istat) et, ',', phi
                    case (FILETYPE_JSON)
                        ! for this we accumulate the data for the whole
                        ! segment and write all at once at the end
                        phi_vec = [phi_vec, phi] ! TODO: could use expand routines to make this more efficient
                        phi_se = [phi_se, p] ! save the actual value
                        et_se = [et_se, et]

                        ! also convert to sun-earth rotating
                        se_rotating = two_body_rotating_frame(primary_body=body_sun,&
                                                              secondary_body=body_earth,&
                                                              center=center_at_secondary_body,&
                                                              et=et)
                        ! from inertial to rotating:
                        inertial = icrf_frame(b=body_moon)
                        call inertial%transform(rv_moon,se_rotating,et,me%segs(iseg)%eph,x_se_rotating,status_ok)
                        if (.not. status_ok) error stop 'transformation error in propagate_segment'
                        x_se = [x_se, x_se_rotating(1)]
                        y_se = [y_se, x_se_rotating(2)]
                        z_se = [z_se, x_se_rotating(3)]

                    end select

                end associate
            end do
            if (ifiletype==FILETYPE_JSON) then
                ! accumulate only the points where phi < 0 (in shadow)
                et_vec  = pack(me%segs(iseg)%traj_inertial%et, mask=phi_vec<0.0_wp)
                phi_vec = pack(phi_vec, mask=phi_vec<0.0_wp)
                if (size(et_vec)>0) then
                    et_vec_total = [et_vec_total, et_vec]
                    phi_vec_total = [phi_vec_total, phi_vec]
                end if
            end if

        end do

        ! save file:
        select case (ifiletype)
        case (FILETYPE_CSV)
            close(iunit,iostat=istat)
        case (FILETYPE_JSON)
            call json%add('et',  et_vec_total)     ! TODO: could remove any duplicate time entries (e.g., at the segment interfaces)
            call json%add('phi', phi_vec_total)

            call json%add('phi_se', phi_se) ! also the full S-E trajectory data for plotting
            call json%add('et_se', et_se)
            call json%add('x_se_rotating', x_se)
            call json%add('y_se_rotating', y_se)
            call json%add('z_se_rotating', z_se)
            call json%print(full_filename)
        end select
        call json%destroy()

    end subroutine generate_eclipse_data
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
    real(wp),dimension(6) :: x_se_rotating !! the state to export
    type(icrf_frame) :: inertial
    type(two_body_rotating_frame) :: rotating
    type(two_body_rotating_frame) :: se_rotating
    logical :: status_ok !! transformation status flag

    select type (me)
    class is (segment)

        et = me%data%et_ref + t  ! convert to ephemeris time [sec]

        ! convert state to the moon-centered rotating frame
        inertial = icrf_frame(b=body_moon)
        rotating = two_body_rotating_frame(primary_body=body_earth,&
                                           secondary_body=body_moon,&
                                           center=center_at_secondary_body,&
                                           et=et)
        ! from inertial to rotating:
        call inertial%transform(x,rotating,et,me%eph,x_rotating,status_ok)
        if (.not. status_ok) error stop 'transformation error in propagate_segment'

        ! also convert to sun-earth rotating for plotting: -- really should only do when when exporting the json traj file
        se_rotating = two_body_rotating_frame(primary_body=body_sun,&
                                              secondary_body=body_earth,&
                                              center=center_at_secondary_body,&
                                              et=et)
        ! from inertial to rotating:
        call inertial%transform(x,se_rotating,et,me%eph,x_se_rotating,status_ok)
        if (.not. status_ok) error stop 'transformation error in propagate_segment'

        ! save both the inertial and rotating trajectories:
        call me%traj_inertial%add(et,x)
        call me%traj_rotating%add(et,x_rotating)
        call me%traj_se_rotating%add(et,x_se_rotating)

    class default
        error stop 'invalid class in trajectory_export_func'
    end select

    end subroutine trajectory_export_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the epoch info from a config file.

    subroutine read_epoch(f, epoch_mode, year, month, day, hour, minute, sec, et)

    type(config_file),intent(inout) :: f
    integer,intent(out) :: year, month, day, hour, minute !! only output if found (epoch_mode=1)
    real(wp),intent(out) :: sec                           !! only output if found (epoch_mode=1)
    integer,intent(out) :: epoch_mode !! 1 : calendar date was specified, 2: ephemeris time was specified
    real(wp),intent(out) :: et !! the ephemeris time [always output]

    logical :: found, found_et
    logical,dimension(6) :: found_calendar

    ! can specify either calendar date or et:
    call f%get('year',   year   , found_calendar(1) )
    call f%get('month',  month  , found_calendar(2) )
    call f%get('day',    day    , found_calendar(3) )
    call f%get('hour',   hour   , found_calendar(4) )
    call f%get('minute', minute , found_calendar(5) )
    call f%get('sec',    sec    , found_calendar(6) )
    call f%get('et_ref', et     , found_et )

    found = all(found_calendar) .or. found_et ! at least one
    if (found) found = .not. (all(found_calendar) .and. found_et) ! only one
    if (.not. found) error stop 'error: just specify epoch as year,month,day,hour,minute,sec or et_ref'
    if (found_et) then
        epoch_mode = 2 ! ephemeris time was specified
    else
        epoch_mode = 1 ! calendar date was specified
    end if

    if (epoch_mode==1) then
        ! have to convert to et [see update_epoch]
        et = jd_to_et(julian_date(year,&
                                  month,&
                                  day,&
                                  hour,&
                                  minute,&
                                  sec))
    end if

    end subroutine read_epoch
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the config file that defines the problem
!  to be solved and set all the global variables.
!
!### Notes
!  * "period" is the normalized period ("jc", the jacobi constant, is also allowed)
!  * The states in the patch point file are assumed to be "S" family

    subroutine read_config_file(me, filename)

    implicit none

    class(my_solver_type),intent(inout) :: me
    character(len=*),intent(in) :: filename  !! the JSON config file to read

    type(config_file) :: f
    logical :: found, found_jc, found_period, found_initial_r
    real(wp):: jc !! jacobii constant from the file
    real(wp):: p !! period from file (normalized)
    type(patch_point),dimension(3) :: pp !! patch points
    logical :: use_json_file  !! if true, use the JSON file.
                              !! if false, use the CSV file.
    real(wp) :: lstar, tstar
    real(wp),dimension(:),allocatable :: jcvec, normalized_period
    real(wp),dimension(:),allocatable :: x0, z0, ydot0
    real(wp) :: pp_period, pp_x0, pp_z0, pp_ydot0

    write(*,*) '* Loading config file: '//trim(filename)
    call f%open(filename)

    write(*,*) '* Reading config file'

    call f%get('jc',                        jc, found_jc)          ! one of these two must be present
    call f%get('period',                    p, found_period)       !
    call f%get('solve',                             me%mission%solve,                     found)
    call f%get('fix_initial_time',                  me%mission%fix_initial_time,          found)
    call f%get('fix_initial_r',                     me%mission%fix_initial_r,             found)
    call f%get('fix_ry_at_end_of_rev',              me%mission%fix_ry_at_end_of_rev,      found)
    call f%get('fix_final_ry_and_vx',               me%mission%fix_final_ry_and_vx,       found)
    call f%get('generate_plots',                    me%mission%generate_plots,            found)
    call f%get('generate_trajectory_files',         me%mission%generate_trajectory_files, found)
    call f%get('generate_json_trajectory_file',     me%mission%generate_json_trajectory_file,     found)
    call f%get('generate_guess_and_solution_files', me%mission%generate_guess_and_solution_files, found)
    call f%get('generate_kernel',                   me%mission%generate_kernel,           found)
    call f%get('generate_defect_file',              me%mission%generate_defect_file,      found)
    call f%get('generate_eclipse_files',            me%mission%generate_eclipse_files,    found)
    call f%get('run_pyvista_script',                me%mission%run_pyvista_script,        found)
    call f%get('generate_rp_ra_file',               me%mission%generate_rp_ra_file,       found)

    ! if generating the kernel: have to generate the trajectory files
    if (me%mission%generate_kernel) me%mission%generate_trajectory_files = .true.

    call f%get('constrain_initial_rdot', me%mission%constrain_initial_rdot, found)

    call f%get('r_eclipse_bubble',    me%mission%r_eclipse_bubble, found)
    call f%get('eclipse_dt_step',     me%mission%eclipse_dt_step,  found)
    call f%get('eclipse_filetype',    me%mission%eclipse_filetype, found)

    call f%get('use_splined_ephemeris',  me%mission%use_splined_ephemeris,  found)
    call f%get('dt_spline_sec',          me%mission%dt_spline_sec,          found)

    call f%get('pointmass_central_body',  me%mission%pointmass_central_body,  found)
    call f%get('include_pointmass_earth',    me%mission%include_pointmass_earth,    found)
    call f%get('include_pointmass_sun',      me%mission%include_pointmass_sun,      found)
    call f%get('include_pointmass_jupiter',  me%mission%include_pointmass_jupiter,  found)

    call f%get('initial_guess_from_file', me%mission%initial_guess_from_file, found)

    call f%get('solver_mode', me%mission%solver_mode, found)
    if (.not. found) me%mission%solver_mode = 1

    !optional global parameter inputs:
    call f%get('mu_earth'  , mu_earth,  found)
    call f%get('mu_moon'   , mu_moon,   found)
    call f%get('mu_sun'    , mu_sun,    found)
    call f%get('mu_jupiter', mu_jupiter, found)
    call f%get('rad_moon'  , rad_moon,  found)
    call f%get('rad_sun'   , rad_sun,   found)
    call f%get('rad_earth' , rad_earth, found)
    call f%get('maxnum'    , maxnum,    found)
    call f%get('grav_n'    , grav_n,    found)
    call f%get('grav_m'    , grav_m,    found)
    call f%get('grav_frame', grav_frame,found)
    call f%get('xscale_x0' , xscale_x0, found) ! scale values (defaults are good values for the 4500 Rp case)
    if (.not. found) xscale_x0 = [1.0e+05_wp,1.0e+05_wp,1.0e+05_wp,2.0e+00_wp,2.0e+00_wp,2.0e+00_wp]
    if (size(xscale_x0) /= 6) error stop 'error: xscale_x0 must be a 6 element vector'
    call f%get('fscale_xf', fscale_xf, found)
    if (.not. found) fscale_xf = [1.0e+04_wp,1.0e+04_wp,1.0e+04_wp,1.0e+02_wp,1.0e+02_wp,1.0e+02_wp]
    if (size(fscale_xf) /= 6) error stop 'error: fscale_xf must be a 6 element vector'
    call f%get('fscale_rdot' , fscale_rdot, found) ! rdot function scale value

    call f%get('use_battin_gravity', use_battin_gravity, found)

    call f%get('object_id', object_id, found)
    call f%get('object_name', object_name, found); if (.not. found) object_name = 'HALO'
    call f%get('leapseconds_file', leapseconds_file, found); if (.not. found) leapseconds_file = 'kernel/naif0012.tls'
    call f%get('mkspk_path', mkspk_path, found); if (.not. found) mkspk_path = 'kernel/mkspk'
    call f%get('polynom_degree', polynom_degree, found)
    call f%get('output_spk_type', output_spk_type, found)
    call f%get('segment_id', segment_id, found); if (.not. found) segment_id = 'SPK_STATES_09'
    call f%get('min_export_time_step', min_export_time_step, found)
    min_export_time_step = abs(min_export_time_step) ! only non-negative values allowed

    ! required inputs:
    call f%get('N_or_S',   me%mission%N_or_S   )
    call f%get('L1_or_L2', me%mission%L1_or_L2 )

    ! reach the epoch. can specify either calendar date or et:
    call read_epoch(f,  me%mission%epoch_mode, &
                        me%mission%year, &
                        me%mission%month, &
                        me%mission%day, &
                        me%mission%hour, &
                        me%mission%minute, &
                        me%mission%sec, &
                        me%mission%et_ref)

    ! tolerances:
    call f%get('rtol',          me%mission%rtol          , found )
    call f%get('atol',          me%mission%atol          , found )
    call f%get('nlesolver_tol', me%mission%nlesolver_tol , found )

    ! if fix_initial_r is true, can specify the r to use
    ! [otherwise, the initial guess is used from the patch points]
    call f%get('initial_r', me%mission%initial_r , found_initial_r )

    call f%get('n_revs',          me%mission%n_revs           )
    call f%get('ephemeris_file',  me%mission%ephemeris_file   )
    call f%get('gravfile',        me%mission%gravfile         )
    call f%get('patch_point_file',me%mission%patch_point_file )
    call f%get('patch_point_file_is_periapsis',me%mission%patch_point_file_is_periapsis, found )

    call f%get('moon_pa_file', me%mission%moon_pa_file, found)
    if (.not. found) me%mission%moon_pa_file = 'data/moon_pa_2000_2100.csv'

    call f%close() ! cleanup ---------------------------------------------------

    use_json_file = index(me%mission%patch_point_file, '.json') > 0
    if (.not. use_json_file) error stop 'error: patch point file must be a JSON file.'
    if (use_json_file .and. (found_jc .eqv. found_period)) &
        error stop 'error: must either specify period or jc in the config file to use JSON patch point file.'
        ! csv file no longer supported

    ! set the epoch:
    call me%mission%update_epoch()

    ! read the patch point file and load it:

    ! read the JSON file:
    write(*,*) '* Reading JSON patch point file: '//trim(me%mission%patch_point_file)
    call f%open(me%mission%patch_point_file)

    ! all inputs are required:
    write(*,*) '* getting data... '
    call f%get('lstar', lstar            )
    call f%get('tstar', tstar            )
    call f%get('jc',    jcvec            )
    call f%get('period',normalized_period)
    call f%get('x0',    x0               )
    call f%get('z0',    z0               )
    call f%get('ydot0', ydot0            )
    call f%close()

    ! write(*,'(A, *(F6.2,",",1X))') 'period range in database (days): ', &
    !       normalized_period * tstar * sec2day

    ! which is the independant variable:
    if (found_jc) then
        pp_period = interpolate_point(jcvec, normalized_period, jc)
        pp_x0     = interpolate_point(jcvec, x0, jc)
        pp_z0     = interpolate_point(jcvec, z0, jc)
        pp_ydot0  = interpolate_point(jcvec, ydot0, jc)
    else
        pp_period = p
        pp_x0     = interpolate_point(normalized_period, x0, p)
        pp_z0     = interpolate_point(normalized_period, z0, p)
        pp_ydot0  = interpolate_point(normalized_period, ydot0, p)
    end if

    call me%mission%generate_patch_points(lstar, tstar, pp_period, &
                                          pp_x0, pp_z0, pp_ydot0, pp)
    me%mission%periapsis = pp(1)
    me%mission%quarter   = pp(2)
    me%mission%apoapsis  = pp(3)

    ! write(*,'(a15,1x,F10.6,1x,a)') 'periapsis t: ' , me%mission%periapsis%t, 'days'
    ! write(*,'(a15,1x,F10.6,1x,a)') 'quarter t:   ' , me%mission%quarter%t,   'days'
    ! write(*,'(a15,1x,F10.6,1x,a)') 'apoapsis t:  ' , me%mission%apoapsis%t,  'days'

    ! compute some time variables:
    me%mission%period  = me%mission%apoapsis%t * 2.0_wp  ! Halo period [days]
    me%mission%period8 = me%mission%period / 8.0_wp      ! 1/8 of Halo period [days]

    write(*,*) ''
    write(*,'(A30,1x,f16.6,a,i4,a,i2,a,i2,a,i0.2,a,i0.2,a,f9.6,a)') '   Reference ephemeris time: ', &
                                                                          me%mission%et_ref, &
                                                                    ' (', me%mission%year,'-',&
                                                                          me%mission%month,'-',&
                                                                          me%mission%day,' ',&
                                                                          me%mission%hour,':',&
                                                                          me%mission%minute,':',&
                                                                          me%mission%sec, ')'
    write(*,'(a30,1x,f21.17)') '   Orbit period (days): ', me%mission%period

    contains
!*****************************************************************************************

    !**********************************************
    !>
    !  Iterpolate if necessary to compute `fvec(x)`.
    !  We will allow data sets that are not strictly increasing/decreasing
    !  (`jc` isn't, but `period` is for L2 file),
    !  so, first the input `x` is checked to see if that is a point in `xvec`.
    !  If not, the bspline interpolation is used.

        function interpolate_point(xvec, fvec, x) result(y)

        implicit none

        real(wp),dimension(:),intent(in) :: xvec
        real(wp),dimension(:),intent(in) :: fvec
        real(wp),intent(in) :: x
        real(wp) :: y !! `fvec(x)`

        type(bspline_1d) :: spline
        integer :: i, iflag, irow

        integer,parameter :: kx = 4  !! spline order (cubic bspline)
        integer,parameter :: idx = 0 !! interpolate value only

        ! so, first just look for the exact point:
        irow = 0
        do i=1,size(xvec)
            if (xvec(i)==x) then
                irow = i
                exit
            end if
        end do

        if (irow/=0) then
            ! if found, then use it:
            y = fvec(irow)
        else
            call spline%initialize(xvec,fvec,kx,iflag, extrap=.true.)
            if (iflag/=0) then
                select case (iflag)
                case(2); error stop 'bspline error: iknot out of range.'
                case(3); error stop 'bspline error: nx out of range.'
                case(4); error stop 'bspline error: kx out of range.'
                case(5); error stop 'bspline error: x not strictly increasing.'
                end select
            end if
            call spline%evaluate(x,idx,y,iflag)
            if (iflag/=0) error stop 'error evaluating bspline.'
        end if

        end function interpolate_point
    !**********************************************

    end subroutine read_config_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Update the variables for the reference epoch in the mission class.
!  Either computes the et from the calendar state, or vice versa.

    subroutine update_epoch(me)

    implicit none

    class(mission_type),intent(inout) :: me

    select case (me%epoch_mode)
    case(1)
        ! compute reference epoch from the date (TDB):
        ! save this in the mission class
        me%et_ref = calendar_date_to_et(me%year,&
                                        me%month,&
                                        me%day,&
                                        me%hour,&
                                        me%minute,&
                                        me%sec)
    case(2)
        ! then compute calendar from et
        call julian_date_to_calendar_date(et_to_jd(me%et_ref),&
                                          me%year,&
                                          me%month,&
                                          me%day,&
                                          me%hour,&
                                          me%minute,&
                                          me%sec)
    case default
        error stop 'error: epoch_mode must be 1 or 2.'
    end select

    end subroutine update_epoch
!*****************************************************************************************

!*****************************************************************************************
!>
!  Generates the patch points (unnormalized, moon-centered) from the CR3BP normalized guess.

    subroutine generate_patch_points(me, lstar, tstar, period, x0, z0, ydot0, pp)

    implicit none

    class(mission_type),intent(inout) :: me
    real(wp),intent(in) :: lstar
    real(wp),intent(in) :: tstar
    real(wp),intent(in) :: period
    real(wp),intent(in) :: x0
    real(wp),intent(in) :: z0
    real(wp),intent(in) :: ydot0
    type(patch_point),dimension(3), intent(out) :: pp !! periapsis, quarter, apoapsis patch points

    integer,parameter  :: n = 6 !! number of state variables

    real(wp),dimension(n) :: x_wrt_moon_normalized   !! normalized state wrt to moon
    real(wp),dimension(n) :: x_wrt_moon_unnormalized !! unnormalized state wrt to moon
    real(wp),dimension(n) :: x                       !! unnormalized state wrt barycenter, for integration
    real(wp),dimension(:),allocatable :: t_crtbp, x_crtbp,y_crtbp,z_crtbp,vx_crtbp,vy_crtbp,vz_crtbp
    type(ddeabm_class) :: prop  !! integrator
    real(wp) :: t_unnormalized ! unnormalized time
    real(wp) :: mu    !! CRTPB parameter
    real(wp) :: t    !! normalized time for integration
    real(wp) :: tf   !! final time (normalized) for integration
    real(wp) :: dt   !! time step (normalized) for integration
    integer :: idid  !! ddeabm output flag
    integer :: i     !! counter
    real(wp) :: t_, tf_ !! temp time vars for advancing guess by 1/2 period

    ! preallocate the arrays used in the report function:
    allocate(t_crtbp(0),x_crtbp(0),y_crtbp(0),z_crtbp(0),&
             vx_crtbp(0),vy_crtbp(0),vz_crtbp(0))

    ! first generate the patch points in the normalized system wrt the barycenter.
    ! do this with a numerical integration of the CR3BP equations of motion.

    mu = compute_crtpb_parameter(mu_earth,mu_moon)
    t  = zero          ! start at periapsis
    tf = period / two  ! propagate to apoapsis
    dt = period / four ! time step of 1/4 period
    x = [x0, zero, z0, zero, ydot0, zero] ! initial state: normalized wrt barycenter

    ! integrate and report the points at dt steps (periapsis, 1/4, and apoapsis)
    call prop%initialize(n,maxnum=maxnum,df=func,&
                         rtol=[me%rtol],atol=[me%atol],&
                         report=report)
    if (.not. me%patch_point_file_is_periapsis) then
        ! advance state by 1/2 period to put x at periapsis
        t_ = t ! make temp copies so the originals aren't changed
        tf_ = tf
        call prop%first_call()
        call prop%integrate(t_,x,tf_,idid=idid,integration_mode=1)
    end if
    call prop%first_call()
    call prop%integrate(t,x,tf,idid=idid,integration_mode=2,tstep=dt)

    ! normalized patch points are now in x_crtbp,y_crtbp,z_crtbp,vx_crtbp,vy_crtbp,vz_crtbp
    ! write(*,*) ''
    ! write(*,*) 'generate_patch_points'
    ! write(*,'(A,*(F10.6,1X))') '  t_crtbp  = ', t_crtbp
    ! write(*,'(A,*(F10.6,1X))') '  x_crtbp  = ', x_crtbp
    ! write(*,'(A,*(F10.6,1X))') '  y_crtbp  = ', y_crtbp
    ! write(*,'(A,*(F10.6,1X))') '  z_crtbp  = ', z_crtbp
    ! write(*,'(A,*(F10.6,1X))') '  vx_crtbp = ', vx_crtbp
    ! write(*,'(A,*(F10.6,1X))') '  vy_crtbp = ', vy_crtbp
    ! write(*,'(A,*(F10.6,1X))') '  vz_crtbp = ', vz_crtbp
    ! write(*,*) ''

    ! transform state to moon-centered and unnormalize
    do i = 1, 3 ! periapsis, quarter, apoapsis

        ! transform state to moon-centered rotating frame

        !            y
        !            ^
        !            |
        !  (M1)------+----------------(M2)  --> x
        !       -mu           1-mu
        !

        x_wrt_moon_normalized = [x_crtbp(i) - (1.0_wp - mu), &
                                 y_crtbp(i), &
                                 z_crtbp(i), &
                                 vx_crtbp(i), &
                                 vy_crtbp(i), &
                                 vz_crtbp(i)]

        ! convert to km, km/s, days (see also unnormalize_variables)
        x_wrt_moon_unnormalized(1:3) = x_wrt_moon_normalized(1:3) * lstar         ! unscale distance
        x_wrt_moon_unnormalized(4:6) = x_wrt_moon_normalized(4:6) * (lstar/tstar) ! unscale velocity
        t_unnormalized = t_crtbp(i) * tstar ! unscale time

        ! results:
        select case (me%N_or_S)
        case ('S')
            ! the data in the file is for the South family
            pp(i) = patch_point(t = t_unnormalized*sec2day,&
                                rv = x_wrt_moon_unnormalized)
        case ('N')
            x_wrt_moon_unnormalized(3) = -x_wrt_moon_unnormalized(3)
            x_wrt_moon_unnormalized(6) = -x_wrt_moon_unnormalized(6)
            pp(i) = patch_point(t = t_unnormalized*sec2day,&
                                rv = x_wrt_moon_unnormalized)
        case default
            error stop 'invalid value for N_or_S: '//trim(me%N_or_S)
        end select

    end do

    contains

    subroutine func(me,t,x,xdot)  !! CRTBP derivative function
        implicit none
        class(ddeabm_class),intent(inout) :: me
        real(wp),intent(in)               :: t
        real(wp),dimension(:),intent(in)  :: x
        real(wp),dimension(:),intent(out) :: xdot

        call crtbp_derivs(mu,x,xdot)

    end subroutine func

    subroutine report(me,t,x)  !! report function
        implicit none
        class(ddeabm_class),intent(inout)    :: me
        real(wp),intent(in)                  :: t
        real(wp),dimension(:),intent(in)     :: x

        t_crtbp   = [t_crtbp,  t]
        x_crtbp   = [x_crtbp,  x(1)]
        y_crtbp   = [y_crtbp,  x(2)]
        z_crtbp   = [z_crtbp,  x(3)]
        vx_crtbp  = [vx_crtbp, x(4)]
        vy_crtbp  = [vy_crtbp, x(5)]
        vz_crtbp  = [vz_crtbp, x(6)]

    end subroutine report

    end subroutine generate_patch_points
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a string that can be used for file names, etc.
!  (read the values of the global variables).
!
!@note We are only saving `sec` as an integer.

    function get_case_name(me) result(case_name)

    implicit none

    class(mission_type),intent(in) :: me
    character(len=:),allocatable :: case_name

    case_name = int_to_string(me%year,4)//int_to_string(me%month,2)//&
                int_to_string(me%day,2)//int_to_string(me%hour,2)//&
                int_to_string(me%minute,2)//int_to_string(int(me%sec),2)//&
                '_'//me%L1_or_L2//'_'//me%N_or_S//'_NREVS='//int_to_string(me%n_revs)

    end function get_case_name
!*****************************************************************************************

!*****************************************************************************************
!>
!  integer to string with optional zero padding

    function int_to_string(i,ipad) result(str)
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
!*****************************************************************************************

!*****************************************************************************************
!>
!  Write the current `x` variables to a file (compatible with the Copernicus optvar file).

    subroutine write_optvars_to_file(me,filename,x)

    implicit none

    class(mission_type),intent(inout) :: me
    character(len=*),intent(in)       :: filename !! file name [without extension]
    real(wp),dimension(:),intent(in)  :: x        !! solver opt vars vector

    integer :: i  !! counter
    type(json_core) :: json
    type(json_value),pointer :: p_root, p_xvec, p_element, p_last

    call json%initialize()
    call json%create_object(p_root, '')
    call json%create_array(p_xvec, 'xvec')
    call json%add(p_root, p_xvec)
    p_last => null() ! for the first element
    do i = 1, size(x)
        nullify(p_element)
        call json%create_object(p_element, '') ! allocate
        call json%add(p_element, 'i',     i)
        call json%add(p_element, 'label', trim(me%xname(i)))
        call json%add(p_element, 'value', x(i)*me%xscale(i))
        call json%add(p_element, 'scale', me%xscale(i))
        ! to speed this up, we keep track of the
        ! last element and insert after it.
        if (associated(p_last)) then
            ! insert after the last one
            call json%insert_after(p_last, p_element)
        else
            ! first element
            call json%add(p_xvec, p_element)
        end if
        p_last => p_element
    end do

    call json%print(p_root, trim(filename)//'_'//me%get_case_name()//'.json')

    call json%destroy(p_root)

    end subroutine write_optvars_to_file
!*****************************************************************************************

!*****************************************************************************************
    end module halo_module
!*****************************************************************************************
