!*****************************************************************************************
!>
!  Main program to solve the Halo targeting problem.
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
!     "period": 1.5872714606,
!     "N_or_S": "S",
!     "L1_or_L2": "L2",
!     "year": 2000,
!     "month": 1,
!     "day": 1,
!     "hour": 12,
!     "minute": 0,
!     "sec": 0,
!     "n_revs": 10,
!     "ephemeris_file":   "data/eph/JPLEPH.421",
!     "gravfile":         "data/grav/gggrx_0020pm_sha.tab",
!     "patch_point_file": "data/L2_halos.json"
! }
!```
!
!  Notes:
!
!  * "period" is the normalized period ("jc", the jacobi constant, is also allowed)
!  * The states in the patch point file are assumed to be "S" family
!
!### Author
!  * Jacob Williams : Sept. 2017

    program halo_solver

    use parameters_module
    use argv_module,       only: argv
    use halo_module,       only: my_solver_type,define_problem_size
!$  use omp_lib

    implicit none

    logical,parameter :: debug = .false. !! for debugging prints

    character(len=:),allocatable :: config_file_name  !! the config file to read
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
    character(len=:),allocatable :: mkspk_input, bsp_output  !! filenames for mkspk
!$  integer :: tid, nthreads

!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!$
!$  tid = omp_get_thread_num()
!$
!$  if (tid == 0) then
!$      nthreads = omp_get_num_threads()
!$      write(*,*) '----------------------'
!$      write(*,'(A,1X,I5)') 'number of OMP threads: ', OMP_get_num_threads()
!$      write(*,*) '----------------------'
!$  end if
!$
!$OMP END PARALLEL

    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) 'Initializing...'
    write(*,*) '----------------------'
    write(*,*) ''

    if (command_argument_count() /=1 ) &
        error stop 'Error: the first command line argument must be the config file name'

    config_file_name = argv(1)

    call solver%init(config_file_name,x)  ! initialize the solver & mission (and generate the initial guess)

    if (allocated(solver%mission%initial_guess_from_file)) then
        !TODO: add some error checking here !
        write(*,*) 'Reading initial guess from file: '//solver%mission%initial_guess_from_file
        call solver%mission%get_x_from_json_file(x) ! get solution from the file
        call solver%mission%put_x_in_segments(x) ! populate segs with solution
    end if

    if (solver%mission%generate_plots) &
        call solver%mission%plot('guess', draw_trajectory=.true.)    ! plot the initial guess
    if (solver%mission%generate_guess_and_solution_files) &
        call solver%mission%write_optvars_to_file('guess',x)    ! write guess to a file

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

    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) 'Solving...'
    write(*,*) '----------------------'
    write(*,*) ''

    call cpu_time(tstart_cpu)
!$  tstart = omp_get_wtime()
    call solver%solve(x)  ! call the solver
    call solver%status(istat=istat, message=message)
    call cpu_time(tend_cpu)
!$  tend = omp_get_wtime()

    call solver%mission%put_x_in_segments(x) ! populate segs with solution

    write(*,*) ''
    write(*,*) '----------------------'
    write(*,*) 'Solution...'
    write(*,*) '----------------------'
    write(*,*) ''

    write(*,'(A)') 'status: '//message
    write(*,*) 'elapsed cpu_time: ', (tend_cpu-tstart_cpu), 'sec'
!$  write(*,*) 'OMP wall time   : ', (tend-tstart), 'sec'
    write(*,*) ''

    if (solver%mission%generate_guess_and_solution_files) &
        call solver%mission%write_optvars_to_file('solution',x) ! write solution to a file
    ! export solution to plot or trajectory file
    if (solver%mission%generate_trajectory_files .or. solver%mission%generate_plots) &
        call solver%mission%plot('solution',&
                draw_trajectory=solver%mission%generate_plots, &
                export_trajectory=solver%mission%generate_trajectory_files)
    if (solver%mission%generate_kernel) then
        if (.not. solver%mission%generate_trajectory_files) then
            write(*,*) 'error: kernel generation requires the trajectory file to be exported'
        else
            mkspk_input = 'solution_'//solver%mission%get_case_name()//'.txt'
            bsp_output  = 'solution_'//solver%mission%get_case_name()//'.bsp'
            call execute_command_line('kernel/mkspk -setup kernel/setup.txt -input '//mkspk_input//' -output '//bsp_output)
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
        call solver%mission%print_constraint_defects('solution_defects_'//&
                                                     solver%mission%get_case_name()//&
                                                     '.csv')
    end if

    if (solver%mission%generate_eclipse_files) then
        call solver%mission%generate_eclipse_data('eclipse', &
                                                  filetype = solver%mission%eclipse_filetype)
    end if

!*****************************************************************************************
    end program halo_solver
!*****************************************************************************************
