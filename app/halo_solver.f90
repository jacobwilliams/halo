!*****************************************************************************************
!>
!  Main program to solve the Halo targeting problem.
!
!### Usage
!
!```
!    > halo_solver config.json
!```
!
!### Author
!  * Jacob Williams : Sept. 2017

    program halo_solver

    use argv_module,       only: argv
    use halo_module,       only: halo_solver_main

    implicit none

    logical,parameter :: debug = .false. !! for debugging prints
    character(len=:),allocatable :: config_file_name  !! the config file to read

    ! get the config file name for the command line argument:
    if (command_argument_count() /=1 ) &
        error stop 'Error: the first command line argument must be the config file name'
    config_file_name = argv(1)

    ! the main routine:
    call halo_solver_main(config_file_name,debug)

!*****************************************************************************************
    end program halo_solver
!*****************************************************************************************
