!*****************************************************************************************
!>
!  For handling JSON config files

 module config_file_module
    use json_module,       only: json_file
    use parameters_module, only: wp

    implicit none

    private

    type,public :: config_file
        !! a class for handling JSON config files.
        private
        character(len=:),allocatable :: name !! the file name
        type(json_file),public :: json !! made it public so we can use the methods
        contains
        private
        procedure,public :: open  => open_config_file
        procedure,public :: close => close_config_file
        generic,public :: get => get_int, get_real, get_logical, get_char,&
                                 get_real_vec
        procedure,private :: get_int, get_real, get_logical, get_char, get_real_vec
    end type config_file

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Open the file.

    subroutine open_config_file(me, filename)

    class(config_file),intent(inout) :: me
    character(len=*),intent(in) :: filename !! the name of the JSON file

    call me%json%load(filename=filename)
    if (me%json%failed()) error stop 'error reading json file '//trim(filename)

    end subroutine open_config_file
!*****************************************************************************************

!*****************************************************************************************
!>
!  Open the file.

    subroutine close_config_file(me)

    class(config_file),intent(inout) :: me

    call me%json%destroy()

    end subroutine close_config_file
!*****************************************************************************************

!**********************************************
!>
!  Get an integer variable from the config file.

    subroutine get_int(me, name, variable, found)

    class(config_file),intent(inout) :: me
    character(len=*),intent(in) :: name !! name of the variable
    integer,intent(inout) :: variable  !! the variable to populate with the value
    logical,intent(out),optional :: found !! if the variable was in the file

    integer :: tmp
    logical :: required  !! if the variable is required
                         !! (if required and not found, `error stop` is called)
    logical :: was_found !! local copy of `found`

    required = .not. present(found)

    call me%json%get(name, tmp, was_found); if (was_found) variable = tmp

    if (required .and. .not. was_found) &
    error stop trim(name)//' integer variable not found in config file: '//trim(me%name)
    if (present(found)) found = was_found

    end subroutine get_int
!**********************************************

!**********************************************
!>
!  Get a real variable from the config file.

    subroutine get_real(me, name, variable, found)

    class(config_file),intent(inout) :: me
    character(len=*),intent(in) :: name !! name of the variable
    real(wp),intent(inout) :: variable  !! the variable to populate with the value
    logical,intent(out),optional :: found !! if the variable was in the file

    real(wp) :: tmp
    logical :: required  !! if the variable is required
                         !! (if required and not found, `error stop` is called)
    logical :: was_found !! local copy of `found`

    required = .not. present(found)

    call me%json%get(name, tmp, was_found); if (was_found) variable = tmp

    if (required .and. .not. was_found) &
    error stop trim(name)//' integer variable not found in config file: '//trim(me%name)
    if (present(found)) found = was_found

    end subroutine get_real
!**********************************************

!**********************************************
!>
!  Get a logical variable from the config file.

    subroutine get_logical(me, name, variable, found)

    class(config_file),intent(inout) :: me
    character(len=*),intent(in) :: name !! name of the variable
    logical,intent(inout) :: variable  !! the variable to populate with the value
    logical,intent(out),optional :: found !! if the variable was in the file

    logical :: tmp
    logical :: required  !! if the variable is required
                         !! (if required and not found, `error stop` is called)
    logical :: was_found !! local copy of `found`

    required = .not. present(found)

    call me%json%get(name, tmp, was_found); if (was_found) variable = tmp

    if (required .and. .not. was_found) &
    error stop trim(name)//' integer variable not found in config file: '//trim(me%name)
    if (present(found)) found = was_found

    end subroutine get_logical
!**********************************************

!**********************************************
!>
!  Get a character variable from the config file.

    subroutine get_char(me, name, variable, found)

    class(config_file),intent(inout) :: me
    character(len=*),intent(in) :: name !! name of the variable
    character(len=:),allocatable,intent(inout) :: variable  !! the variable to populate with the value
    logical,intent(out),optional :: found !! if the variable was in the file

    character(len=:),allocatable :: tmp
    logical :: required  !! if the variable is required
                         !! (if required and not found, `error stop` is called)
    logical :: was_found !! local copy of `found`

    required = .not. present(found)

    call me%json%get(name, tmp, was_found); if (was_found) variable = tmp

    if (required .and. .not. was_found) &
    error stop trim(name)//' integer variable not found in config file: '//trim(me%name)
    if (present(found)) found = was_found

    end subroutine get_char
!**********************************************

!**********************************************
!>
!  Get a real vector from the config file.

    subroutine get_real_vec(me, name, variable, found)

    class(config_file),intent(inout) :: me
    character(len=*),intent(in) :: name !! name of the variable
    real(wp),dimension(:),allocatable,intent(inout) :: variable  !! the variable to populate with the value
    logical,intent(out),optional :: found !! if the variable was in the file

    real(wp),dimension(:),allocatable :: tmp
    logical :: required  !! if the variable is required
                         !! (if required and not found, `error stop` is called)
    logical :: was_found !! local copy of `found`

    required = .not. present(found)

    call me%json%get(name, tmp, was_found); if (was_found) variable = tmp

    if (required .and. .not. was_found) &
    error stop trim(name)//' integer variable not found in config file: '//trim(me%name)
    if (present(found)) found = was_found

    end subroutine get_real_vec
!**********************************************

!*****************************************************************************************
    end module config_file_module
!*****************************************************************************************