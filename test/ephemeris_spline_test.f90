!*****************************************************************************************
!>
!  Time the setup of the splined ephemeris (and of the splined MOON_PA frame) for
!  the problem sizes of the technical note, with the ephemeris options of a config
!  file (default `examples/example_sparse.json`: `ephemeris_file`, `dt_spline_sec`,
!  `moon_pa_file`, the epoch and the period). The time span is the solver's
!  (see [[initialize_the_mission]]): from two periods before the reference epoch to
!  two periods after the last rev.
!
!  Timed, for each number of revs (the minimum of `n_repeat` runs):
!
!  * opening the JPL ephemeris file (`jpl_ephemeris%initialize`),
!  * the splined ephemeris (`initialize_splinded_ephemeris`: opening the file, and
!    evaluating DE421 at each spline point for the Earth, Sun, SSB and Jupiter
!    w.r.t. the Moon, and fitting the splines),
!  * the splined MOON_PA frame (`moon_pa_frame%initialize`: reading the CSV
!    file and fitting the splines; only used with `grav_frame=2`),
!  * the whole problem setup (`mission%init`, with the config's options: it includes
!    the above, and the per-segment copies made for OpenMP).
!
!  Usage (from the halo root directory):
!
!    fpm test ephemeris_spline_test --profile release --flag "-fopenmp" -- [config.json] [n_revs ...]
!
!  The default numbers of revs are 100 and 1120.

    program ephemeris_spline_test

    use halo_module
    use halo_kinds_module, only: wp
    use splined_ephemeris_module, only: jpl_ephemeris_splined
    use splined_moon_frame_module, only: moon_pa_frame
    use fortran_astrodynamics_toolkit, only: jpl_ephemeris
    use iso_fortran_env, only: output_unit
!$  use omp_lib

    implicit none

    integer,parameter :: n_repeat = 3  !! number of runs of each (the minimum time is reported)

    type(my_solver_type) :: solver
    character(len=1024) :: arg
    character(len=:),allocatable :: config
    integer,dimension(:),allocatable :: revs
    real(wp),dimension(:),allocatable :: x
    real(wp) :: et0, etf, t_file, t_spline, t_moon, t_init, t0
    integer :: i, j, n_points
    logical :: status_ok

    config = 'examples/example_sparse.json'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        config = trim(arg)
    end if
    if (command_argument_count() >= 2) then
        allocate(revs(command_argument_count()-1))
        do i = 2, command_argument_count()
            call get_command_argument(i, arg)
            read(arg,*) revs(i-1)
        end do
    else
        revs = [100, 1120]
    end if

    call solver%read_config_file(config)  ! (also sets the reference epoch)

    write(output_unit,'(A)') ''
    write(output_unit,'(A)') ' config:         '//config
    write(output_unit,'(A)') ' ephemeris file: '//solver%mission%ephemeris_file
    write(output_unit,'(A,F8.1,A)') ' spline step:    ', solver%mission%dt_spline_sec, ' sec'
    write(output_unit,'(A)') ' moon_pa file:   '//solver%mission%moon_pa_file
    write(output_unit,'(A,I0,A)') ' (minimum of ', n_repeat, ' runs)'

    do j = 1, size(revs)

        ! the solver's time span for this number of revs:
        associate (m => solver%mission)
            m%n_revs = revs(j)
            et0 = m%et_ref - 2.0_wp*m%period*day2sec
            etf = m%et_ref + m%period*day2sec*m%n_revs + 2.0_wp*m%period*day2sec
            n_points = floor((etf-et0)/m%dt_spline_sec) + 1
        end associate

        t_file = huge(1.0_wp); t_spline = huge(1.0_wp); t_moon = huge(1.0_wp); t_init = huge(1.0_wp)
        do i = 1, n_repeat
            t_file   = min(t_file,   time_file())
            t_spline = min(t_spline, time_spline())
            t_moon   = min(t_moon,   time_moon_pa())
            t_init   = min(t_init,   time_init())
        end do

        write(output_unit,'(A)') ''
        write(output_unit,'(A,I0,A,F9.1,A,I0,A)') ' revs: ', revs(j), '   span: ', (etf-et0)/day2sec, &
                                                  ' days   spline points: ', n_points, ' (per coordinate)'
        write(output_unit,'(A,F9.4,A)') '   open the ephemeris file:     ', t_file,   ' s'
        write(output_unit,'(A,F9.4,A)') '   splined ephemeris (total):   ', t_spline, ' s'
        write(output_unit,'(A,F9.4,A)') '   splined MOON_PA frame:       ', t_moon,   ' s'
        write(output_unit,'(A,F9.4,A)') '   whole setup (mission%init):  ', t_init,   ' s'

    end do

    contains

        function time_file() result(t)
            !! open the JPL ephemeris file (not splined)
            real(wp) :: t
            type(jpl_ephemeris),allocatable :: eph
            allocate(eph)
            t0 = wall_time()
            call eph%initialize(filename=solver%mission%ephemeris_file, status_ok=status_ok)
            t = wall_time() - t0
            if (.not. status_ok) error stop 'error opening the ephemeris file'
            call eph%close()
        end function time_file

        function time_spline() result(t)
            !! the splined ephemeris, as in the solver
            real(wp) :: t
            type(jpl_ephemeris_splined),allocatable :: eph
            allocate(eph)
            t0 = wall_time()
            call eph%initialize_splinded_ephemeris(filename=solver%mission%ephemeris_file, status_ok=status_ok, &
                                                   et0=et0, dt=solver%mission%dt_spline_sec, etf=etf)
            t = wall_time() - t0
            if (.not. status_ok) error stop 'error initializing the splined ephemeris'
            call eph%close()
        end function time_spline

        function time_moon_pa() result(t)
            !! the splined MOON_PA frame, as in the solver
            real(wp) :: t
            type(moon_pa_frame),allocatable :: moon_pa
            allocate(moon_pa)
            t0 = wall_time()
            call moon_pa%initialize(solver%mission%moon_pa_file, et0=et0, etf=etf)
            t = wall_time() - t0
            call moon_pa%destroy()
        end function time_moon_pa

        function time_init() result(t)
            !! the whole problem setup, with the config's options
            real(wp) :: t
            type(my_solver_type) :: s
            call s%read_config_file(config)
            s%mission%n_revs = solver%mission%n_revs
            t0 = wall_time()
            call s%mission%init(x)
            t = wall_time() - t0
        end function time_init

        function wall_time() result(t)
            real(wp) :: t
            integer(8) :: count, rate
!$          t = omp_get_wtime()
!$          return
            call system_clock(count, rate)
            t = real(count,wp) / real(rate,wp)
        end function wall_time

    end program ephemeris_spline_test
!*****************************************************************************************
