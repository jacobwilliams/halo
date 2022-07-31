!*****************************************************************************************
!>
!  Main program to solve the NRHO targeting problem.
!
!### Usage
!
!```
!    > halo-solver config.json
!```
!
!### Example config file
!
!```json
! {
!     "period": 1.5872714606,  // normalized period ("jc" is also allowed)
!     "N_or_S": "S",
!     "L1_or_L2": "L2",
!     "year": 2000,
!     "month": 1,
!     "day": 1,
!     "hour": 12,
!     "minute": 0,
!     "sec": 0,
!     "n_revs": 10,
!     "ephemeris_file": "../data/eph/JPLEPH_windows_ifort.421"
! }
!
!### Author
!  * Jacob Williams : Sept. 2017

    program halo_solver

    use parameters_module
    use argv_module,       only: argv
    use halo_module,       only: my_solver_type,nrho_func,&
                                 nrho_grad,nrho_export,&
                                 define_problem_size

    implicit none

    logical,parameter :: debug = .false.

    character(len=:),allocatable :: config_file_name  !! the config file to read
    type(my_solver_type) :: solver  !! an instance of the solver that we will use
    real(wp),dimension(:),allocatable :: x  !! solver opt vars vector ["forward-backward" formulation]
    integer :: m  !! number of functions
    real(wp),dimension(:),allocatable :: f  !! function vector (constraint violations)
    real :: tstart, tend  !! for timing
    integer :: istat
    character(len=:),allocatable :: message  !! Text status message from solver

    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) 'Initializing...'
    write(*,*) '----------------------'
    write(*,*) ''

    if (command_argument_count() /=1 ) &
        error stop 'Error: the first command line argument must be the config file name'

    config_file_name = argv(1)

    call solver%init(config_file_name,x)  ! initialize the solver & mission (and generate the initial guess)
    if (generate_plots) &
        call solver%mission%plot('guess')  ! plot the initial guess
    if (generate_trajectory_files) &
        call solver%mission%write_optvars_to_file('guess',x) ! write guess to a file
    call define_problem_size(m=m)
    allocate(f(m))
    call solver%mission%constraint_violations(x,f)

    if (debug) then
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

    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) 'Solving...'
    write(*,*) '----------------------'
    write(*,*) ''

    call cpu_time(tstart)
    call solver%solve(x)  ! call the solver
    call solver%status(istat=istat, message=message)
    call cpu_time(tend)

    call solver%mission%put_x_in_segments(x) ! populate segs with solution

    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) 'Solution...'
    write(*,*) '----------------------'
    write(*,*) ''

    write(*,'(A)')             'status: '//message
    write(*,*) 'elapsed time: ', (tend-tstart), 'sec'
    write(*,*) ''

    if (generate_trajectory_files) &
        call solver%mission%write_optvars_to_file('solution',x)       ! write solution to a file:
    if (generate_plots) &
        call solver%mission%plot('solution',export_trajectory=.true.) ! plot solution

!*****************************************************************************************
    end program halo_solver
!*****************************************************************************************
