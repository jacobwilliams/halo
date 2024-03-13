!*****************************************************************************************
!>
!  Main program to find an eclipse-free halo orbit.
!  This program calls the solver as a command line call and reads the outputs.
!  It uses a simulated annealing algorithm to minimize the time in earth eclipse.
!
!  See: `run-sa-main.sh`
!
! TODO: replace with a single program.
!
!### Notes
!  * could enable the initial_guess_from_file output for iterations after the first one [not currently done]
!  * could switch to real128 for the last iteration? or after the final trajectory is found?

    program halo_optimizer

    use halo_kinds_module,          only: wp
    use argv_module,                only: argv
    use simulated_annealing_module, only: simulated_annealing_type
    use json_module,                only: json_file
    use config_file_module,         only: config_file
    use halo_module

    implicit none

    type(simulated_annealing_type) :: sa
    real(wp), dimension(2) :: x !! [period, et_ref]
    real(wp) :: f
    integer :: iunit, istat

    integer,parameter :: n = 2 !! the 2 optimization variables are `et_ref` and `period` in the config file
    character(len=*),parameter :: base_config_file = './examples/example_sparse.json' ! TODO: should be an input
    character(len=*),parameter :: base_script_file = './run-sa.sh' ! the script to run for the function evaluation

    integer :: ns   !! number of cycles.
    integer :: nt  !! number of iterations before temperature reduction.
    integer :: neps  !! number of final function values used to decide upon termination.
    integer :: n_resets  !! number of times to run the main loop
    real(wp),dimension(:),allocatable :: lb !! lower bounds (size `n`)
    real(wp),dimension(:),allocatable :: ub !! upper bounds (size `n`)
    real(wp),dimension(:),allocatable :: vm !! size `n`
    real(wp) :: fopt, t, rt, period
    real(wp),dimension(n) :: xopt
    integer  :: nacc, nfcnev, ier
    logical :: found
    type(config_file) :: config
    integer :: epoch_mode, year, month, day, hour, minute
    real(wp) :: sec, et_ref

    call config%open(base_config_file)
    call config%json%print()

    open(newunit = iunit, file = 'function_evals.csv', status='REPLACE', iostat=istat)
    write(iunit,'(A20,1X,2(A30,1X),A8)') 'CONFIG', 'PERIOD', 'ET', 'OBJFCN'

    ! initial guess:
    !x = [1.5872714606_wp, 757339200.0_wp] ! period, et [2024-01-01]
    call config%get('period', period) ! required
    call read_epoch(config, epoch_mode, year, month, day, hour, minute, sec, et_ref)
    x = [period, et_ref]

    ! defaults for optional arguments:
    t        = 5.0_wp
    lb       = [1.40_wp, x(2)-7*day2sec]  ! * vary epoch by +/- 7 days
    ub       = [1.60_wp, x(2)+7*day2sec]
    vm       = [0.2_wp, 7*day2sec]
    ! set some small values for testing...
    ns       = 3      ! number of cycles.
    nt       = 3      ! number of iterations before temperature reduction.
    neps     = 2      ! number of final function values used to decide upon termination.
    n_resets = 1      ! number of times to run the main loop
    rt       = 0.5_wp ! temperature reduction factor

    ! if these are specified in the config file then use those instead:
    call config%get('sa_t',       t,        found)
    call config%get('sa_lb',      lb,       found)
    call config%get('sa_ub',      ub,       found)
    call config%get('sa_vm',      vm,       found)
    call config%get('sa_ns',      ns,       found)
    call config%get('sa_nt',      nt,       found)
    call config%get('sa_neps',    neps,     found)
    call config%get('sa_n_resets',n_resets, found)
    call config%get('sa_rt',      rt,       found)

    write(*,*) 'Running optimizer'
    call sa%initialize(fcn,n,lb,ub,&
                       n_resets = n_resets, &
                       ns       = ns, &
                       nt       = nt, &
                       neps     = neps, &
                       optimal_f_specified = .true.,&
                       optimal_f     = 0.0_wp,& ! 0 is the desired solution (no eclipses) if it's found, we can stop.
                       optimal_f_tol = epsilon(1.0_wp)) ! since this is really an int, this doesn't matter much
    call sa%optimize(x, rt, t, vm, xopt, fopt, nacc, nfcnev, ier)

    write(*,*) ''
    write(*,*) '=============='
    write(*,*) '   SOLUTION   '
    write(*,*) '=============='
    write(*,*) ''
    write(*,*) 'xopt = ', xopt
    write(*,*) 'fopt = ', fopt
    write(*,*) 'ier  = ', ier
    write(*,*) ''

    !TODO: maybe make a new config for the solution ?

    close(iunit, iostat=istat) ! function eval file
    call config%close()

    contains

        subroutine fcn(me, x, f, istat)

            !! function for the SA algorithm.

            class(simulated_annealing_type),intent(inout) :: me
            real(wp), dimension(:), intent(in) :: x !! [period, et_ref]
            real(wp), intent(out) :: f
            integer,intent(out) :: istat

            integer,save :: istep = 0 !! number of function evaluations

            character(len=5) :: istep_str
            character(len=:),allocatable :: script !! script to run
            character(len=:),allocatable :: step_config !! config file for this step
            character(len=:),allocatable :: eclipse_file !! eclipse file
            type(json_file) :: json !! for loading the eclipse json file
            real(wp),dimension(:),allocatable :: et_vec !! ephemeris time of eclipse violations
            type(my_solver_type) :: solver !! just so we can read the base config and generate the case name

            istep = istep + 1
            write(istep_str,'(I5)') istep

            ! create the script
            ! read the original script and modify with the inputs, then save the new one
            step_config = './STEP='//trim(adjustl(istep_str))//'.json'
            call config%json%update('period', x(1), found)
            call config%json%update('et_ref', x(2), found)
            ! make sure the output settings are right:
            call config%json%update('initial_guess_from_file', '', found)
            call config%json%update('generate_plots',                    .false., found)
            call config%json%update('generate_trajectory_files',         .false., found)
            call config%json%update('generate_guess_and_solution_files', .false., found)
            call config%json%update('generate_kernel',                   .false., found)
            call config%json%update('generate_defect_file',              .false., found)
            call config%json%update('run_pyvista_script',                .false., found)
            call config%json%update('generate_eclipse_files', .true., found)  ! only need these
            call config%json%update('eclipse_filetype', 2, found)             !
            call config%json%remove('year')
            call config%json%remove('month')
            call config%json%remove('day')
            call config%json%remove('hour')
            call config%json%remove('minute')
            call config%json%remove('sec')
            call config%json%print(step_config)

            ! have to read the config so that get_case_name will return the right string:
            write(*,*) 'reading base config file: '//base_config_file
            call solver%read_config_file(step_config)

            ! the env var is used by the script to indicate the config file to use:
            call execute_command_line('export HALO_CONFIG_FILE='//step_config//'; '//base_script_file)

            ! read the eclipse function output file
            eclipse_file = './eclipse_'//solver%mission%get_case_name()//'.json'
            call json%load(eclipse_file)
            call json%get('et', et_vec, found)
            if (.not. found) error stop 'error reading eclipse file. et not found.'
            call json%destroy()

            ! f is just the count the 1-hr epochs where the phi is < 0
            !TODO:: could use unique() to remove duplidate et's. just count the unique ets (can be zero)
            f = real(size(et_vec), wp)

            istat = 0 ! no problems

            write(*,*) '======================'
            write(*,*) ' Latest f: ', f
            write(*,*) '======================'
            write(iunit,'(A20,1X,2(E30.16,1X),I8)') step_config, x, int(f)

        end subroutine fcn

    end program halo_optimizer