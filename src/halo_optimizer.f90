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
    use halo_module

    implicit none

    type(simulated_annealing_type) :: sa
    real(wp), dimension(2) :: x !! [period, et_ref]
    real(wp) :: f
    integer :: iunit, istat

    integer,parameter :: n = 2 !! the 2 optimization variables are `et_ref` and `period` in the config file

    ! simulated annealing parameters -- set some small values for testing...
    integer,parameter  :: ns = 3  !! number of cycles.
    integer,parameter  :: nt = 3 !! number of iterations before temperature reduction.
    integer,parameter  :: neps = 2 !! number of final function values used to decide upon termination.
    integer,parameter  :: n_resets = 1 !! number of times to run the main loop
    real(wp),parameter :: rt = 0.5_wp

    real(wp),dimension(n) :: lb !! lower bounds
    real(wp),dimension(n) :: ub !! upper bounds
    real(wp) :: fopt, t
    real(wp),dimension(n) :: xopt
    real(wp),dimension(n) :: vm
    integer  :: nacc, nfcnev, ier

    character(len=*),parameter :: base_config_file = './examples/example_sparse.json'
    character(len=*),parameter :: base_script_file = './run-sa.sh'

    open(newunit = iunit, file = 'function_evals.csv', status='REPLACE', iostat=istat)

    ! initial guess:
    x = [1.5872714606_wp, 710212455.0_wp] ! period, et
    ! lower and upper bounds
    !   * vary epoch by +/- 7 days
    lb = [x(1)-0.2_wp, x(2)-7*day2sec]
    ub = [x(1)+0.2_wp, x(2)+7*day2sec]
    t  = 5.0_wp
    vm = 1.0_wp

    write(*,*) 'Running optimizer'
    call sa%initialize(fcn,n,lb,ub,&
                       n_resets = n_resets, &
                       ns       = ns, &
                       nt       = nt, &
                       neps     = neps, &
                       optimal_f_specified = .true.,&
                       optimal_f           = 0.0_wp,& ! 0 is the desired solution (no eclipses) if it's found, we can stop.
                       optimal_f_tol       = epsilon(1.0_wp))
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

    close(iunit, iostat=istat) ! function eval file

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
            character(len=:),allocatable :: config !! config file for this step
            character(len=:),allocatable :: eclipse_file !! eclipse file
            type(json_file) :: json
            real(wp),dimension(:),allocatable :: et_vec !! ephemeris time of eclipse violations
            logical :: found
            type(my_solver_type) :: solver !! just so we can read the base config and generate the case name

            istep = istep + 1
            write(istep_str,'(I5)') istep

            ! create the script
            ! read the original script and modify with the inputs, then save the new one
            config = './STEP='//trim(adjustl(istep_str))//'.json'
            call json%load(base_config_file)
            call json%update('period', x(1), found)   ! should we scale these?
            call json%update('et_ref', x(2), found)   !
            !TODO: make sure the output settings are right (have to generate the eclipse json files ....)
            call json%print(config)
            call json%destroy()

            ! have to read the config so that get_case_name will return the right string:
            write(*,*) 'reading base config file: '//base_config_file
            call solver%read_config_file(config)

            ! the env var is used by the script to indicate the config file to use:
            call execute_command_line('export HALO_CONFIG_FILE='//config//'; '//base_script_file)

            ! read the eclipse function output file
            eclipse_file = './eclipse_'//solver%mission%get_case_name()//'.json'
            call json%load(eclipse_file)
            call json%get('et', et_vec)
            call json%destroy()

            ! f is just the count the 1-hr epochs where the phi is < 0
            !TODO:: could use unique() to remove duplidate et's. just count the unique ets (can be zero)
            f = real(size(et_vec), wp)

            istat = 0 ! no problems

            write(*,*) '======================'
            write(*,*) ' Latest f: ', f
            write(*,*) '======================'
            write(iunit,'(A,1X,*(E30.16,1X))') config, x, f

        end subroutine fcn

    end program halo_optimizer