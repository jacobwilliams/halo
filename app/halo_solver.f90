!*****************************************************************************************
!>
!  Main program to solve the Halo targeting problem.
!
!### Usage
!
!```
!    > halo_solver config.json [--gui]
!```
!
!  `--gui` is used by `halo_gui`: the JSON trajectory file is always
!  exported (and its name written to the output), and the pyvista
!  script is not run.
!
!### Author
!  * Jacob Williams : Sept. 2017

    program halo_solver

    use argv_module,       only: argv
    use halo_module,       only: halo_solver_main

    implicit none

    logical,parameter :: debug = .false. !! for debugging prints
    character(len=:),allocatable :: config_file_name  !! the config file to read
    logical :: gui !! if run from the GUI

    ! get the config file name for the command line argument:
    gui = .false.
    select case (command_argument_count())
    case (1)
    case (2)
        if (argv(2) /= '--gui') error stop 'Error: the only option after the config file name is --gui'
        gui = .true.
    case default
        error stop 'Error: the first command line argument must be the config file name'
    end select
    config_file_name = argv(1)

    ! the main routine:
    call halo_solver_main(config_file_name,debug,gui)

!*****************************************************************************************
    end program halo_solver
!*****************************************************************************************
