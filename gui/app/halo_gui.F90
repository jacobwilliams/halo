!*****************************************************************************************
!>
!  A GUI for the Halo solver: choose a config file (e.g. `examples/example_sparse.json`),
!  run the solver, and plot the resulting trajectory in 3D. Optionally, each
!  iteration of the solver is plotted as it goes.
!
!  The window is made with [qt-fortran](https://github.com/jacobwilliams/qt-fortran),
!  and the 3D view with [vkframes](https://github.com/jacobwilliams/vkframes).
!
!  * The solver runs in this program, in a worker thread (`halo_solver_main`, with
!    `gui=.true.`: the JSON trajectory file is always written, and the pyvista script
!    is not run), in the current directory (the config files' paths are relative to
!    it). Its output (standard output, redirected to a pipe while it runs) is shown
!    in the "Solver output" panel (a reader thread reads the pipe), and in the
!    terminal. Stop asks the solver to stop at its next iteration.
!  * "Plot each iteration": the solver's iteration callback gives each iteration's
!    trajectory (in memory) to the GUI thread, which puts it straight into the
!    vkframes trajectories. When the solver is done, the solution's `traj_*.json`
!    file is plotted.
!  * File > Open Trajectory plots an existing `traj_*.json` file without solving.
!  * The trajectory is drawn in the Moon-centered Earth-Moon rotating frame or the
!    Moon-centered inertial (J2000) frame (as in `python/plot_utilities.py`), colored by
!    the distance from the Moon, with a body following it in time.
!
!  Notes:
!  * the solver runs in this process, so if it stops the program (`error stop`,
!    e.g. for a bad config file), the GUI closes too.
!  * large arrays here are allocatable (on the heap), not automatic or
!    temporaries: with OpenMP, local arrays are on the stack, and worker threads'
!    stacks are small (macOS: 512 kB; Windows: 1 MB).
!  * the POSIX calls for the pipe are the C runtime's `_pipe`, `_dup`, ... on Windows.
!
!### Usage
!
!```
!    > pixi run gui -- [config.json] [--run] [--iterations] [--traj traj_file.json]
!```
!
!### Author
!  * Jacob Williams : Nov. 2026

    module halo_gui_module

    use qt
    use vkframes
    use vkframes_qt,        only: connect_view, ViewLink, DemoOptions
    use vkf_text,           only: get_font
    use json_module,        only: json_core, json_value
    use halo_module,        only: halo_solver_main
    use halo_kinds_module,  only: hwp => wp  ! (the solver's real kind)
    use iso_c_binding
    use iso_fortran_env,    only: output_unit, int32

    implicit none

    private

    public :: app, build_scene, build_ui, start_up, change_directory

    integer,parameter :: slider_steps = 1000  !! the time slider's range
    real(wp),parameter :: moon_radius = 1737.4_wp  !! km
    real(wp),parameter :: day = 86400.0_wp  !! sec

    !> the frames the trajectory can be shown in
    integer,parameter :: ROTATING = 1, INERTIAL = 2
    character(len=*),parameter :: frame_names(2) = [character(len=32) :: &
                                    'Earth-Moon rotating (EMR)', 'Inertial (J2000)']

    !> colormap for the distance from the Moon: near (red) to far (blue)
    real(wp),parameter :: cmap(3,5) = reshape([ 1.00_wp, 0.25_wp, 0.20_wp, &
                                                1.00_wp, 0.75_wp, 0.15_wp, &
                                                0.55_wp, 1.00_wp, 0.35_wp, &
                                                0.20_wp, 0.85_wp, 1.00_wp, &
                                                0.35_wp, 0.45_wp, 1.00_wp ], [3,5])

    !> The trajectory in one frame, and what draws it.
    type :: frame_t
        type(Trajectory) :: traj
        type(DrawableTrajectory) :: drawable
        type(CurveArtist) :: curve
        real(wp) :: center(3) = 0.0_wp  !! the center of the trajectory's bounding box (km)
        real(wp) :: radius = 1.0e5_wp  !! half its diagonal (km)
    end type frame_t

    !> One iteration's trajectory, made by the worker thread and handed to the
    !> GUI thread (which deletes it) through `solver_job%iteration`.
    type :: iteration_points
        integer :: iter = 0
        real(wp),dimension(:),allocatable :: t  !! (n) ET (sec)
        real(wp),dimension(:,:),allocatable :: rot, inr  !! (3,n) km
    end type iteration_points

    !> A solver run (in the worker thread). The worker reads the inputs, sets the
    !> results, and emits `iteration`; it touches nothing else of the GUI.
    type :: solver_job
        ! inputs
        character(len=:),allocatable :: config  !! the config file
        logical :: iterations = .false.  !! plot each iteration
        ! results
        character(len=:),allocatable :: traj_file  !! the solution's trajectory file (if written)
        ! iteration points: (the address of an iteration_points, in two halves)
        type(QSignalInt2) :: iteration
    end type solver_job

    !> The reader thread's job: read the pipe `fd` until it ends, sending what
    !> it reads to the GUI thread through `text`.
    type :: reader_job
        integer(c_int) :: fd = -1  !! the pipe's read end
        type(QSignalString) :: text
    end type reader_job

    !> The GUI, the scene and the solver.
    type :: app_t
        ! window
        type(QMainWindow) :: main
        type(QStackedWidget) :: stack  !! the center: a startup message, then the view
        type(QLabel) :: startup_message
        type(QWindow) :: window
        type(QWidget) :: container
        type(QtfEventFilter) :: events
        type(ViewLink) :: view_link
        type(WindowProxy) :: proxy
        integer :: frame = 0  !! frames drawn
        ! controls
        type(QLineEdit) :: config_edit
        type(QPushButton) :: run_button, stop_button, browse_button
        type(QCheckBox) :: plot_iterations
        type(QLabel) :: status_label, info_label
        type(QComboBox) :: frame_choice
        type(QDoubleSpinBox) :: time_scale
        type(QCheckBox) :: pause_time
        type(QLabel) :: time_label
        type(QSlider) :: time_slider
        logical :: slider_held = .false.
        logical :: slider_setting = .false.
        type(QPlainTextEdit) :: log
        ! the solver (worker thread)
        type(QThread) :: worker
        logical :: running = .false.
        logical :: plotted_this_run = .false.  !! something was plotted in this run (keep the camera)
        type(solver_job) :: job
        ! its output: standard output redirected to a pipe while it runs
        integer(c_int) :: saved_stdout = -1  !! the original standard output (while redirected)
        integer(c_int) :: tee_fd = -1  !! a copy of it, the output is also written to (until the pipe ends)
        type(QThread) :: reader  !! reads the pipe (`read_pipe`)
        logical :: reading = .false.  !! the reader is running
        type(reader_job) :: rjob
        character(len=:),allocatable :: out_buffer  !! output not yet shown (no newline yet)
        ! the scene
        type(FrameManager) :: scene
        type(ReferenceFrame) :: root  !! Moon-centered
        type(Sphere) :: moon
        type(Sphere) :: body  !! follows the trajectory
        type(frame_t) :: frames(2)  !! rotating, inertial
        integer :: current_frame = ROTATING
        type(View) :: view
        type(SkySphere) :: sky
        type(HUDText) :: hud_title, hud_frame, hud_time, hud_near, hud_far
        logical :: loaded = .false.  !! a trajectory is loaded
        real(wp) :: t0 = 0.0_wp, tf = 1.0_wp  !! the trajectory's times (ET sec)
        ! startup options
        character(len=:),allocatable :: config_file, initial_traj
        logical :: run_at_start = .false.  !! --run: run the solver at startup
        logical :: iterations_at_start = .false.  !! --iterations: "Plot each iteration" checked
        type(DemoOptions) :: options  !! --capture, --quit-after (for testing)
    end type app_t

    type(app_t),target,save :: app

    integer,parameter :: pipe_size = 65536  !! the pipe's (and the reader's buffer) size

    ! the C runtime, for redirecting standard output to a pipe:
    interface
#ifdef _WIN32
        function c_chdir(path) bind(C, name='_chdir') result(r)
            import :: c_int, c_char
            character(kind=c_char),dimension(*),intent(in) :: path
            integer(c_int) :: r
        end function c_chdir
        function c_pipe_win(fds, psize, textmode) bind(C, name='_pipe') result(r)
            import :: c_int
            integer(c_int),dimension(2),intent(out) :: fds
            integer(c_int),value :: psize, textmode
            integer(c_int) :: r
        end function c_pipe_win
        function c_dup(fd) bind(C, name='_dup') result(r)
            import :: c_int
            integer(c_int),value :: fd
            integer(c_int) :: r
        end function c_dup
        function c_dup2(fd, fd2) bind(C, name='_dup2') result(r)
            import :: c_int
            integer(c_int),value :: fd, fd2
            integer(c_int) :: r
        end function c_dup2
        function c_close(fd) bind(C, name='_close') result(r)
            import :: c_int
            integer(c_int),value :: fd
            integer(c_int) :: r
        end function c_close
        function c_read(fd, buf, n) bind(C, name='_read') result(r)
            import :: c_int, c_char
            integer(c_int),value :: fd
            character(kind=c_char),dimension(*),intent(out) :: buf
            integer(c_int),value :: n
            integer(c_int) :: r
        end function c_read
        function c_write(fd, buf, n) bind(C, name='_write') result(r)
            import :: c_int, c_char
            integer(c_int),value :: fd
            character(kind=c_char),dimension(*),intent(in) :: buf
            integer(c_int),value :: n
            integer(c_int) :: r
        end function c_write
#else
        function c_chdir(path) bind(C, name='chdir') result(r)
            import :: c_int, c_char
            character(kind=c_char),dimension(*),intent(in) :: path
            integer(c_int) :: r
        end function c_chdir
        function c_pipe_posix(fds) bind(C, name='pipe') result(r)
            import :: c_int
            integer(c_int),dimension(2),intent(out) :: fds
            integer(c_int) :: r
        end function c_pipe_posix
        function c_dup(fd) bind(C, name='dup') result(r)
            import :: c_int
            integer(c_int),value :: fd
            integer(c_int) :: r
        end function c_dup
        function c_dup2(fd, fd2) bind(C, name='dup2') result(r)
            import :: c_int
            integer(c_int),value :: fd, fd2
            integer(c_int) :: r
        end function c_dup2
        function c_close(fd) bind(C, name='close') result(r)
            import :: c_int
            integer(c_int),value :: fd
            integer(c_int) :: r
        end function c_close
        function c_read_posix(fd, buf, n) bind(C, name='read') result(r)
            import :: c_int, c_size_t, c_char, c_long
            integer(c_int),value :: fd
            character(kind=c_char),dimension(*),intent(out) :: buf
            integer(c_size_t),value :: n
            integer(c_long) :: r
        end function c_read_posix
        function c_write_posix(fd, buf, n) bind(C, name='write') result(r)
            import :: c_int, c_size_t, c_char, c_long
            integer(c_int),value :: fd
            character(kind=c_char),dimension(*),intent(in) :: buf
            integer(c_size_t),value :: n
            integer(c_long) :: r
        end function c_write_posix
#endif
    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  The scene: the Moon at the center, and the trajectory in each frame (empty
!  until one is loaded), with a body following the one shown.

    subroutine build_scene()

    integer :: i

    app%root = ReferenceFrame('Moon-centered')
    call app%root%show_decorations(.false.)

    app%moon = Sphere('Moon', radius=moon_radius, color=[0.6_wp, 0.6_wp, 0.6_wp])
    call app%moon%set_texture_map('lroc_color_2k.jpg')  ! (vkframes' data/, on VKFRAMES_DATA_PATH)
    call app%moon%show_axes(NO_AXES)
    call app%moon%show_axes_labels(NO_AXES)
    if (.not. app%root%add_child(app%moon)) write(*,'(A)') 'halo_gui: add_child failed'

    do i = 1, size(app%frames)
        associate (f => app%frames(i))
        f%traj = Trajectory()
        f%curve = CurveArtist(f%traj, color=[1.0_wp, 1.0_wp, 1.0_wp], width=2.0_wp)
        call f%curve%set_point_colors(.true.)  ! colored by the distance from the Moon
        f%drawable = DrawableTrajectory(trim(frame_names(i)))
        call f%drawable%show_decorations(.false.)
        call f%drawable%add_artist(f%curve)
        if (.not. app%root%add_child(f%drawable)) write(*,'(A)') 'halo_gui: add_child failed'
        end associate
    end do

    app%body = Sphere('HALO', radius=500.0_wp, color=[1.0_wp, 1.0_wp, 1.0_wp])
    call app%body%show_axes(NO_AXES)
    call app%body%show_axes_labels(NO_AXES)
    if (.not. app%root%add_child(app%body)) write(*,'(A)') 'halo_gui: add_child failed'

    call app%scene%set_frame(app%root)

    app%view = View(app%root, app%moon, frame_type=ABSOLUTE_FRAME)
    call app%view%set_default_view_distance(3.0e5_wp)

    app%hud_title = HUDText('No trajectory loaded', position=[0.01_wp, 0.99_wp], &
                            alignment=ALIGN_LEFT_TOP, size=16.0_wp, color=[0.85_wp, 0.85_wp, 0.85_wp])
    app%hud_frame = HUDText('', position=[0.01_wp, 0.955_wp], &
                            alignment=ALIGN_LEFT_TOP, size=14.0_wp, color=[0.6_wp, 0.8_wp, 1.0_wp])
    app%hud_time = HUDText('', position=[0.01_wp, 0.01_wp], &
                           alignment=ALIGN_LEFT_BOTTOM, size=18.0_wp, color=[1.0_wp, 0.85_wp, 0.2_wp])
    app%hud_near = HUDText('', position=[0.99_wp, 0.045_wp], &
                           alignment=ALIGN_RIGHT_BOTTOM, size=14.0_wp, color=cmap(:,1))
    app%hud_far = HUDText('', position=[0.99_wp, 0.01_wp], &
                          alignment=ALIGN_RIGHT_BOTTOM, size=14.0_wp, color=cmap(:,size(cmap,2)))

    if (.not. app%sky%set_star_data('Stars_HYGv3.txt', min_mag=-2.0_wp, max_mag=6.0_wp, max_stars=40000)) &
        write(*,'(A)') 'halo_gui: no star catalog'

    call show_frame(ROTATING)

    end subroutine build_scene
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read a trajectory JSON file written by the solver (`traj_*.json`: an array
!  `segs` of segments, each with `et` and the positions in each frame), and plot
!  it. The segments are joined in time order (backward-propagated ones
!  reversed), as in `python/plot_utilities.py`.

    logical function load_trajectory(file, keep_view) result(ok)

    character(len=*),intent(in) :: file  !! the file to read
    logical,intent(in) :: keep_view  !! keep the camera and the time (else reset them)

    type(json_core) :: json
    type(json_value),pointer :: root, segs, seg
    real(wp),dimension(:),allocatable :: et, xr, yr, zr, xi, yi, zi
    real(wp),dimension(:),allocatable :: t_all
    real(wp),dimension(:,:),allocatable :: rot_all, inr_all
    integer :: nsegs, iseg, n, ntot, k, kstart, kend, kstep
    logical :: found, status_ok
    character(len=:),allocatable :: msg

    ok = .false.
    call json%initialize()
    call json%load(file=file, p=root)
    if (json%failed()) then
        call json%check_for_errors(status_ok, msg)
        call error_message('Could not read '//file//': '//msg)
        call json%destroy()
        return
    end if
    call json%get(root, 'segs', segs, found)
    if (.not. found) then
        call error_message(file//' has no "segs" array (is it a traj_*.json file?)')
        call json%destroy(root)
        return
    end if
    nsegs = json%count(segs)

    ! join the segments:
    allocate(t_all(0), rot_all(3,0), inr_all(3,0))
    ntot = 0
    do iseg = 1, nsegs
        call json%get_child(segs, iseg, seg)
        call json%get(seg, 'et', et, found);  if (.not. found) cycle
        call json%get(seg, 'x_rotating', xr, found);  if (.not. found) cycle
        call json%get(seg, 'y_rotating', yr, found);  if (.not. found) cycle
        call json%get(seg, 'z_rotating', zr, found);  if (.not. found) cycle
        call json%get(seg, 'x_inertial', xi, found);  if (.not. found) cycle
        call json%get(seg, 'y_inertial', yi, found);  if (.not. found) cycle
        call json%get(seg, 'z_inertial', zi, found);  if (.not. found) cycle
        n = min(size(et), size(xr), size(yr), size(zr), size(xi), size(yi), size(zi))
        if (n == 0) cycle
        if (et(n) >= et(1)) then  ! forward propagated
            kstart = 1;  kend = n;  kstep = 1
        else                      ! backward propagated: reversed
            kstart = n;  kend = 1;  kstep = -1
        end if
        call grow(ntot + n)
        do k = kstart, kend, kstep
            ! skip points not after the last one (the segments' shared ends)
            if (ntot > 0) then
                if (et(k) <= t_all(ntot)) cycle
            end if
            ntot = ntot + 1
            t_all(ntot) = et(k)
            rot_all(1,ntot) = xr(k);  rot_all(2,ntot) = yr(k);  rot_all(3,ntot) = zr(k)
            inr_all(1,ntot) = xi(k);  inr_all(2,ntot) = yi(k);  inr_all(3,ntot) = zi(k)
        end do
    end do
    call json%destroy(root)

    if (ntot < 2) then
        call error_message('No trajectory points in '//file)
        return
    end if

    call plot_points(t_all(1:ntot), rot_all(:,1:ntot), inr_all(:,1:ntot), basename(file), &
                     int_text(nsegs)//' segments', keep_view)
    call app%main%setWindowTitle('HALO - '//basename(file))
    ok = .true.

    contains

        subroutine grow(nmax)  !! make room for nmax points (copies on the heap)
            integer,intent(in) :: nmax
            integer :: m
            real(wp),dimension(:),allocatable :: t_new
            real(wp),dimension(:,:),allocatable :: rot_new, inr_new
            if (size(t_all) >= nmax) return
            m = max(nmax, 2 * size(t_all))
            allocate(t_new(m), rot_new(3,m), inr_new(3,m))
            if (ntot > 0) then
                t_new(1:ntot) = t_all(1:ntot)
                rot_new(:,1:ntot) = rot_all(:,1:ntot)
                inr_new(:,1:ntot) = inr_all(:,1:ntot)
            end if
            call move_alloc(t_new, t_all)
            call move_alloc(rot_new, rot_all)
            call move_alloc(inr_new, inr_all)
        end subroutine grow

    end function load_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Plot a trajectory (in both frames): its points put into the vkframes
!  trajectories, colored by the distance from the Moon.

    subroutine plot_points(t, rot, inr, title, details, keep_view)

    real(wp),dimension(:),intent(in) :: t  !! (n) ET (sec)
    real(wp),dimension(:,:),intent(in) :: rot  !! (3,n) Moon-centered rotating frame (km)
    real(wp),dimension(:,:),intent(in) :: inr  !! (3,n) Moon-centered inertial frame (km)
    character(len=*),intent(in) :: title  !! for the HUD and the panel
    character(len=*),intent(in) :: details  !! for the panel
    logical,intent(in) :: keep_view  !! keep the camera and the time (else reset them)

    real(wp),dimension(:),allocatable :: r
    real(wp),dimension(:,:),allocatable :: colors
    real(wp) :: rmin, rmax, time, rmax_inertial
    integer :: j, n

    n = size(t)
    allocate(r(n), colors(3,n))
    rmax_inertial = 0.0_wp
    do j = 1, n
        r(j) = norm2(rot(:,j))  ! (the same in both frames)
        rmax_inertial = max(rmax_inertial, norm2(inr(:,j)))
    end do
    rmin = minval(r)
    rmax = maxval(r)
    do j = 1, n
        colors(:,j) = colormap((r(j) - rmin) / max(rmax - rmin, 1.0e-9_wp))
    end do

    call app%frames(ROTATING)%traj%clear()
    call app%frames(ROTATING)%traj%add(t, pos=rot, color=colors)
    call app%frames(INERTIAL)%traj%clear()
    call app%frames(INERTIAL)%traj%add(t, pos=inr, color=colors)
    call set_bounds(app%frames(ROTATING), rot)
    call set_bounds(app%frames(INERTIAL), inr)
    call app%body%set_radius(0.01_wp * max(rmax_inertial, rmax))

    app%t0 = t(1)
    app%tf = t(n)
    app%loaded = .true.

    call app%hud_title%set_text(title)
    call app%hud_near%set_text('|r| = '//real_text(rmin, '(F0.0)')//' km')
    call app%hud_far%set_text('|r| = '//real_text(rmax, '(F0.0)')//' km')
    call app%info_label%setText(title//'<br>'// &
                                int_text(n)//' points, '//details//'<br>'// &
                                real_text((app%tf - app%t0)/day, '(F0.2)')//' days<br>'// &
                                '|r|: '//real_text(rmin, '(F0.0)')//' to '//real_text(rmax, '(F0.0)')//' km')
    call show_frame(app%current_frame)

    if (keep_view) then
        ! the same time (within the new range), the same camera:
        if (app%proxy%is_open()) then
            time = max(app%t0, min(app%tf, app%proxy%time()))
            call app%proxy%set_time_range(app%t0, app%tf)
            call app%proxy%set_time(time)
        end if
    else
        call app%time_scale%setValue(real(max(1.0_wp, (app%tf - app%t0) / 60.0_wp), c_double))  ! all in ~1 min
        call update_time_range()
        call reset_view()
    end if

    contains

        subroutine set_bounds(f, pos)  !! the bounding box of the trajectory (including the Moon)
            type(frame_t),intent(inout) :: f
            real(wp),dimension(:,:),intent(in) :: pos
            real(wp),dimension(3) :: lo, hi
            integer :: i
            do i = 1, 3
                lo(i) = min(minval(pos(i,:)), -moon_radius)
                hi(i) = max(maxval(pos(i,:)),  moon_radius)
            end do
            f%center = 0.5_wp * (lo + hi)
            f%radius = 0.5_wp * norm2(hi - lo)
        end subroutine set_bounds

    end subroutine plot_points
!*****************************************************************************************

!*****************************************************************************************
!>
!  The color for a fraction `x` (0 to 1) along the colormap.

    pure function colormap(x) result(c)

    real(wp),intent(in) :: x
    real(wp),dimension(3) :: c

    real(wp) :: s, w
    integer :: i

    s = max(0.0_wp, min(1.0_wp, x)) * (size(cmap,2) - 1)
    i = min(int(s) + 1, size(cmap,2) - 1)
    w = s - (i - 1)
    c = (1.0_wp - w) * cmap(:,i) + w * cmap(:,i+1)

    end function colormap
!*****************************************************************************************

!*****************************************************************************************
!>
!  Show the trajectory in frame `iframe` (ROTATING or INERTIAL), with the body following it.

    subroutine show_frame(iframe)

    integer,intent(in) :: iframe

    integer :: i

    app%current_frame = iframe
    do i = 1, size(app%frames)
        call app%frames(i)%drawable%show_contents(i == iframe .and. app%loaded)
    end do
    call app%body%follow(app%frames(iframe)%traj, data=FOLLOW_POSITION, mode=FOLLOW_LIMIT)
    call app%body%show_contents(app%loaded)
    call app%body%show_name_label(app%loaded)
    if (app%loaded) then
        call app%hud_frame%set_text(trim(frame_names(iframe))//', Moon-centered')
    else
        call app%hud_frame%set_text('')
    end if

    end subroutine show_frame
!*****************************************************************************************

!*****************************************************************************************
!>
!  The camera put back to see the whole trajectory in the frame shown: looking
!  at the center of its bounding box from an isometric direction (+x, +y, +z),
!  z up.

    subroutine reset_view()

    type(View),pointer :: v
    real(wp),dimension(3) :: eye, center
    real(wp),dimension(3),parameter :: direction = [1.0_wp, 1.0_wp, 1.0_wp] / sqrt(3.0_wp)
    real(wp),dimension(3),parameter :: up = [0.0_wp, 0.0_wp, 1.0_wp]

    if (.not. app%loaded) return
    center = app%frames(app%current_frame)%center
    eye = center + 3.0_wp * app%frames(app%current_frame)%radius * direction
    call app%view%set_default_view_parameters(eye, center, up)
    if (.not. app%proxy%is_open()) return
    v => app%proxy%current_view()
    if (.not. associated(v)) return
    call v%set_default_view_parameters(eye, center, up)
    call v%reset_view()

    end subroutine reset_view
!*****************************************************************************************

!*****************************************************************************************
!>
!  Give the trajectory's times to the view (if open), starting at the first.

    subroutine update_time_range()

    if (.not. app%proxy%is_open()) return
    call app%proxy%set_time_range(app%t0, app%tf)
    call app%proxy%set_time_scale(app%time_scale%value())
    call app%proxy%set_time(app%t0)

    end subroutine update_time_range
!*****************************************************************************************

!*****************************************************************************************
!>
!  Open the vkframes view in the window, with the scene.

    subroutine open_view()

    if (.not. app%proxy%open(app%window%width(), app%window%height(), 'halo', &
                             native_window=app%window%winId(), &
                             pixel_ratio=app%window%devicePixelRatio())) then
        call error_message('Could not open the 3D view (Vulkan).')
        return
    end if
    if (app%proxy%num_views() > 0) return  ! (opened again, on a new native window)
    call app%proxy%set_scene(app%scene)
    call app%proxy%set_background_color([0.0_wp, 0.0_wp, 0.0_wp])
    call app%proxy%set_sky(app%sky)
    call app%proxy%add_view(app%view)
    call app%proxy%select_view(1)
    call app%proxy%add_hud_text(app%hud_title)
    call app%proxy%add_hud_text(app%hud_frame)
    call app%proxy%add_hud_text(app%hud_time)
    call app%proxy%add_hud_text(app%hud_near)
    call app%proxy%add_hud_text(app%hud_far)
    call update_time_range()

    end subroutine open_view
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exposed: start the frame loop; on a new native window, open the view again.

    subroutine on_expose(data)

    class(*),intent(inout),optional :: data

    if (.not. app%window%isExposed()) return
    if (.not. app%proxy%is_open() .and. app%proxy%num_views() > 0) then
        call open_view()
        if (app%proxy%is_open()) call app%window%requestUpdate()
        return
    end if
    if (app%proxy%is_open() .and. app%frame == 0) call app%window%requestUpdate()

    end subroutine on_expose
!*****************************************************************************************

!*****************************************************************************************
!>
!  Render a frame, update the time readouts, and ask for the next frame.

    subroutine on_update(data)

    class(*),intent(inout),optional :: data

    real(wp) :: t, dpr
    character(len=:),allocatable :: text

    dpr = app%window%devicePixelRatio()
    if (dpr /= app%proxy%pixel_ratio()) &  ! moved to a screen of another scale
        call app%proxy%resize(nint(app%window%width() * dpr), nint(app%window%height() * dpr), pixel_ratio=dpr)
    if (.not. app%proxy%frame()) return
    app%frame = app%frame + 1

    if (app%loaded) then
        t = app%proxy%time()
        text = 't = '//real_text((t - app%t0)/day, '(F0.3)')//' days  (ET '//real_text(t, '(F0.1)')//')'
        call app%hud_time%set_text(text)
        call app%time_label%setText(text)
        if (.not. app%slider_held) then  ! the slider follows the time (unless dragged)
            app%slider_setting = .true.
            call app%time_slider%setValue(nint(slider_steps * (t - app%t0) / max(app%tf - app%t0, 1.0_wp)))
            app%slider_setting = .false.
        end if
    end if
    if (app%options%after_frame(app%proxy, app%window, app%frame)) return  ! (--capture)
    call app%window%requestUpdate()

    end subroutine on_update
!*****************************************************************************************

!*****************************************************************************************
!>
!  After the main window is first painted: the slow parts of the startup
!  (fonts, the Vulkan device), then the view in place of the message; then the
!  trajectory given on the command line, and the solver run, if asked for.

    subroutine on_startup(data)

    class(*),intent(inout),optional :: data

    type(c_ptr) :: font

    call show_startup_message('Loading fonts...')
    font = get_font('')
    call show_startup_message('Starting Vulkan...')
    call open_view()
    call app%stack%setCurrentIndex(1)
    if (app%proxy%is_open()) call app%window%requestUpdate()
    if (allocated(app%initial_traj)) then
        if (load_trajectory(app%initial_traj, keep_view=.false.)) call set_status('Loaded '//basename(app%initial_traj))
    end if
    if (app%run_at_start) call on_run()

    end subroutine on_startup
!*****************************************************************************************

!*****************************************************************************************
    subroutine show_startup_message(text)
    character(len=*),intent(in) :: text
    call app%startup_message%setText(text)
    call QApplication_processEvents()
    end subroutine show_startup_message
!*****************************************************************************************

!*****************************************************************************************
!>
!  Browse for a config file.

    subroutine on_browse(data)

    class(*),intent(inout),optional :: data

    character(len=:),allocatable :: file

    if (app%running) return
    file = QFileDialog_getOpenFileName('Open HALO config file', dir_of(app%config_edit%text()), &
                                       'Config files (*.json);;All files (*)', app%main)
    if (len(file) > 0) call app%config_edit%setText(relative_path(file))

    end subroutine on_browse
!*****************************************************************************************

!*****************************************************************************************
!>
!  Open an existing trajectory file (`traj_*.json`) and plot it.

    subroutine on_open_trajectory(data)

    class(*),intent(inout),optional :: data

    character(len=:),allocatable :: file

    if (app%running) return
    file = QFileDialog_getOpenFileName('Open HALO trajectory file', QDir_currentPath(), &
                                       'Trajectory files (traj_*.json);;JSON files (*.json);;All files (*)', &
                                       app%main)
    if (len(file) == 0) return
    if (load_trajectory(file, keep_view=.false.)) call set_status('Loaded '//basename(file))

    end subroutine on_open_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Run the solver on the config file, in a worker thread (`solve`), with its
!  output shown in the log.

    subroutine on_run(data)

    class(*),intent(inout),optional :: data

    character(len=:),allocatable :: config
    logical :: exists

    if (app%running .or. app%reading) return
    config = trim(app%config_edit%text())
    if (len(config) == 0) then
        call error_message('Choose a config file first.')
        return
    end if
    inquire(file=config, exist=exists)
    if (.not. exists) then
        call error_message('Config file not found: '//config)
        return
    end if

    call app%log%clear()
    call app%log%appendPlainText('Solving '//config//' (in '//QDir_currentPath()//')')
    call capture_output()

    app%job%config = config
    app%job%iterations = app%plot_iterations%isChecked()
    if (allocated(app%job%traj_file)) deallocate(app%job%traj_file)
    app%plotted_this_run = .false.
    app%running = .true.

    app%worker = QThread_create(solve, app%job, app%main)
    call app%worker%finished%connect(on_solver_finished)
    call app%worker%start()

    call set_controls_enabled(.false.)
    call set_status('Solving '//basename(config)//'...')

    end subroutine on_run
!*****************************************************************************************

!*****************************************************************************************
!>
!  The worker thread: run the solver (only the job is touched here).

    subroutine solve(data)

    class(*),intent(inout),optional :: data

    select type (job => data)
    type is (solver_job)
        if (job%iterations) then
            call halo_solver_main(job%config, .false., gui=.true., &
                                  iteration_callback=send_iteration, stop_check=stop_requested, &
                                  traj_file=job%traj_file)
        else
            call halo_solver_main(job%config, .false., gui=.true., &
                                  stop_check=stop_requested, traj_file=job%traj_file)
        end if
        flush(output_unit)
    end select

    end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver's iteration callback (worker thread): copy the points to a new
!  `iteration_points`, and hand it to the GUI thread (`on_iteration`), which
!  plots and deletes it. Its address goes in the signal's two integers (the
!  signal is queued, so nothing here is shared with the GUI thread).

    subroutine send_iteration(iter, et, r_rotating, r_inertial)

    integer,intent(in) :: iter
    real(hwp),dimension(:),intent(in) :: et
    real(hwp),dimension(:,:),intent(in) :: r_rotating, r_inertial

    type(iteration_points),pointer :: p
    integer(c_intptr_t) :: address
    integer(int32),dimension(2) :: halves
    integer :: i, n

    flush(output_unit)  ! (the iteration's output shown now)
    n = size(et)
    if (n < 2) return
    allocate(p)
    allocate(p%t(n), p%rot(3,n), p%inr(3,n))
    p%iter = iter
    do i = 1, n
        p%t(i)     = real(et(i), wp)
        p%rot(:,i) = real(r_rotating(:,i), wp)
        p%inr(:,i) = real(r_inertial(:,i), wp)
    end do
    address = transfer(c_loc(p), address)
    halves = transfer(address, halves)
    call app%job%iteration%emit(int(halves(1)), int(halves(2)))  ! (delivered on the GUI thread)

    end subroutine send_iteration
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver's stop check (worker thread): Stop was clicked.

    logical function stop_requested()
    stop_requested = QThread_isInterruptionRequested()
    end function stop_requested
!*****************************************************************************************

!*****************************************************************************************
!>
!  An iteration's trajectory from the worker (GUI thread): plot it, then delete it.

    subroutine on_iteration(half1, half2, data)

    integer,intent(in) :: half1, half2
    class(*),intent(inout),optional :: data

    type(iteration_points),pointer :: p
    integer(c_intptr_t) :: address
    type(c_ptr) :: cp

    address = transfer([int(half1,int32), int(half2,int32)], address)
    cp = transfer(address, cp)
    call c_f_pointer(cp, p)

    call plot_points(p%t, p%rot, p%inr, 'Iteration '//int_text(p%iter), &
                     'solver iteration', keep_view=app%plotted_this_run)
    app%plotted_this_run = .true.
    call set_status('Solving: iteration '//int_text(p%iter)//'...')
    deallocate(p)

    end subroutine on_iteration
!*****************************************************************************************

!*****************************************************************************************
    subroutine on_stop(data)
    class(*),intent(inout),optional :: data
    if (.not. app%running) return
    call app%worker%requestInterruption()
    call app%log%appendPlainText('*** stop requested (at the next iteration)')
    call set_status('Stopping at the next iteration...')
    call app%stop_button%setEnabled(.false.)
    end subroutine on_stop
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver is done (GUI thread): output back to the terminal, and the
!  solution's trajectory plotted.

    subroutine on_solver_finished(data)

    class(*),intent(inout),optional :: data

    call restore_output()  ! (the reader then sees the end of the output)
    app%running = .false.
    call app%worker%deleteLater()
    call set_controls_enabled(.true.)

    if (.not. allocated(app%job%traj_file)) then
        call set_status('Solver done, but it wrote no trajectory file.')
        return
    end if
    if (load_trajectory(app%job%traj_file, keep_view=app%plotted_this_run)) &
        call set_status('Done: plotted '//app%job%traj_file)

    end subroutine on_solver_finished
!*****************************************************************************************

!*****************************************************************************************
    subroutine set_controls_enabled(enabled)
    logical,intent(in) :: enabled
    call app%run_button%setEnabled(enabled)
    call app%stop_button%setEnabled(.not. enabled)
    call app%config_edit%setEnabled(enabled)
    call app%browse_button%setEnabled(enabled)
    call app%plot_iterations%setEnabled(enabled)
    end subroutine set_controls_enabled
!*****************************************************************************************

!*****************************************************************************************
!>
!  Standard output (the solver's `write(*,...)`, and the programs it runs) to a
!  pipe, read by a reader thread (`read_pipe`) and shown in the log (and copied
!  to the terminal) on the GUI thread (`on_output`).

    subroutine capture_output()

    integer(c_int),dimension(2) :: fds

    if (app%saved_stdout >= 0 .or. app%reading) return  ! (already)
    if (make_pipe(fds) /= 0) return  ! (then the output just goes to the terminal)
    flush(output_unit)
    app%saved_stdout = c_dup(1_c_int)  ! (-1 if there is none, e.g. no console)
    if (c_dup2(fds(2), 1_c_int) < 0) then
        if (c_close(fds(1)) /= 0 .or. c_close(fds(2)) /= 0) continue
        if (app%saved_stdout >= 0) then
            if (c_close(app%saved_stdout) /= 0) continue
        end if
        app%saved_stdout = -1
        return
    end if
    if (c_close(fds(2)) /= 0) continue  ! (standard output is now the only write end)
    app%tee_fd = -1
    if (app%saved_stdout >= 0) app%tee_fd = c_dup(app%saved_stdout)
    app%out_buffer = ''

    app%rjob%fd = fds(1)
    app%reader = QThread_create(read_pipe, app%rjob, app%main)
    call app%reader%finished%connect(on_reader_finished)
    app%reading = .true.
    call app%reader%start()

    end subroutine capture_output
!*****************************************************************************************

!*****************************************************************************************
!>
!  Standard output back to where it was. That closes the pipe's write end, so
!  the reader gets the rest of the output, then the end.

    subroutine restore_output()

    flush(output_unit)
    if (app%saved_stdout < 0) return
    if (c_dup2(app%saved_stdout, 1_c_int) < 0) continue
    if (c_close(app%saved_stdout) /= 0) continue
    app%saved_stdout = -1

    end subroutine restore_output
!*****************************************************************************************

!*****************************************************************************************
!>
!  The reader thread: read the pipe until it ends (all its write ends closed),
!  sending each chunk to the GUI thread. Blocking reads, so it works with pipes
!  on Windows too (where QSocketNotifier doesn't).

    subroutine read_pipe(data)

    class(*),intent(inout),optional :: data

    character(kind=c_char),dimension(:),allocatable :: buf  ! (heap: a worker's stack is small)
    character(len=:),allocatable :: chunk
    integer :: n, i

    select type (job => data)
    type is (reader_job)
        allocate(buf(pipe_size))
        do
            n = read_fd(job%fd, buf)
            if (n <= 0) exit  ! the end (or an error)
            if (allocated(chunk)) deallocate(chunk)
            allocate(character(len=n) :: chunk)
            do i = 1, n
                chunk(i:i) = buf(i)
            end do
            call job%text%emit(chunk)  ! (a copy, delivered on the GUI thread)
        end do
    end select

    end subroutine read_pipe
!*****************************************************************************************

!*****************************************************************************************
!>
!  The reader is done (GUI thread): the last partial line shown, the pipe closed.

    subroutine on_reader_finished(data)

    class(*),intent(inout),optional :: data

    if (len(app%out_buffer) > 0) call app%log%appendPlainText(app%out_buffer)
    app%out_buffer = ''
    if (c_close(app%rjob%fd) /= 0) continue
    app%rjob%fd = -1
    if (app%tee_fd >= 0) then
        if (c_close(app%tee_fd) /= 0) continue
        app%tee_fd = -1
    end if
    app%reading = .false.
    call app%reader%deleteLater()

    end subroutine on_reader_finished
!*****************************************************************************************

!*****************************************************************************************
!>
!  Output from the reader (GUI thread): copied to the terminal, and its complete
!  lines shown in the log.

    subroutine on_output(text, data)

    character(len=*),intent(in) :: text
    class(*),intent(inout),optional :: data

    integer :: i

    if (app%tee_fd >= 0 .and. len(text) > 0) i = write_fd(app%tee_fd, text)

    app%out_buffer = app%out_buffer//text
    do
        i = index(app%out_buffer, new_line('a'))
        if (i == 0) exit
        call add_line(app%out_buffer(1:i-1))
        app%out_buffer = app%out_buffer(i+1:)
    end do

    contains

        subroutine add_line(line)
            character(len=*),intent(in) :: line
            if (len(line) > 0) then
                if (line(len(line):len(line)) == achar(13)) then  ! (CRLF)
                    call app%log%appendPlainText(line(1:len(line)-1))
                    return
                end if
            end if
            call app%log%appendPlainText(line)
        end subroutine add_line

    end subroutine on_output
!*****************************************************************************************

!*****************************************************************************************
!>
!  The C runtime's pipe, read and write, portably.

    integer function make_pipe(fds) result(r)
    integer(c_int),dimension(2),intent(out) :: fds  !! read end, write end
#ifdef _WIN32
    integer(c_int),parameter :: O_BINARY = int(z'8000', c_int)
    r = c_pipe_win(fds, int(pipe_size, c_int), O_BINARY)
#else
    r = c_pipe_posix(fds)
#endif
    end function make_pipe

    integer function read_fd(fd, buf) result(n)  !! bytes read (0: the end; < 0: an error)
    integer(c_int),intent(in) :: fd
    character(kind=c_char),dimension(:),intent(out) :: buf
#ifdef _WIN32
    n = c_read(fd, buf, int(size(buf), c_int))
#else
    n = int(c_read_posix(fd, buf, int(size(buf), c_size_t)))
#endif
    end function read_fd

    logical function change_directory(dir) result(ok)  !! the current directory set to `dir`
    character(len=*),intent(in) :: dir
    ok = c_chdir(trim(dir)//c_null_char) == 0
    end function change_directory

    integer function write_fd(fd, text) result(n)  !! bytes written (< 0: an error)
    integer(c_int),intent(in) :: fd
    character(len=*),intent(in) :: text
#ifdef _WIN32
    n = c_write(fd, text, int(len(text), c_int))
#else
    n = int(c_write_posix(fd, text, int(len(text), c_size_t)))
#endif
    end function write_fd
!*****************************************************************************************

!*****************************************************************************************
!  view and time controls

    subroutine on_frame_choice(index, data)
    integer,intent(in) :: index
    class(*),intent(inout),optional :: data
    call show_frame(index + 1)
    call reset_view()
    end subroutine on_frame_choice

    subroutine on_reset_view(data)
    class(*),intent(inout),optional :: data
    call reset_view()
    end subroutine on_reset_view

    subroutine on_time_scale(value, data)
    real(c_double),intent(in) :: value
    class(*),intent(inout),optional :: data
    if (app%proxy%is_open()) call app%proxy%set_time_scale(value)
    end subroutine on_time_scale

    subroutine on_pause_time(checked, data)
    logical,intent(in) :: checked
    class(*),intent(inout),optional :: data
    if (app%proxy%is_open()) call app%proxy%pause_time(checked)
    end subroutine on_pause_time

    subroutine on_reset_time(data)
    class(*),intent(inout),optional :: data
    if (app%proxy%is_open()) call app%proxy%set_time(app%t0)
    end subroutine on_reset_time

    subroutine on_time_slider(value, data)
    integer,intent(in) :: value
    class(*),intent(inout),optional :: data
    if (app%slider_setting .or. .not. app%proxy%is_open()) return
    call app%proxy%set_time(app%t0 + (app%tf - app%t0) * real(value, wp) / slider_steps)
    end subroutine on_time_slider

    subroutine on_slider_pressed(data)
    class(*),intent(inout),optional :: data
    app%slider_held = .true.
    end subroutine on_slider_pressed

    subroutine on_slider_released(data)
    class(*),intent(inout),optional :: data
    app%slider_held = .false.
    end subroutine on_slider_released

    subroutine on_quit(data)
    class(*),intent(inout),optional :: data
    call app%main%close()
    end subroutine on_quit
!*****************************************************************************************

!*****************************************************************************************
!>
!  Closing: not while the solver runs (it can only be stopped between
!  iterations); otherwise close the view.

    subroutine on_close(event, data)

    type(QCloseEvent),intent(in) :: event
    class(*),intent(inout),optional :: data

    integer :: answer

    if (app%running) then
        answer = QMessageBox_question('Solver running', &
                    'The solver is still running. Stop it at its next iteration, then quit?', &
                    ior(QMessageBox_Yes, QMessageBox_No), QMessageBox_No, app%main)
        call event%ignore()
        if (answer == QMessageBox_Yes) then
            call on_stop()
            call app%worker%finished%connect(on_quit)  ! (closes once it's done)
        end if
        return
    end if
    call app%proxy%close()

    end subroutine on_close
!*****************************************************************************************

!*****************************************************************************************
!>
!  The main window: the 3D view in the center, the controls in a dock on the
!  left, and the solver's output in a dock at the bottom.

    subroutine build_ui()

    type(QDockWidget) :: dock
    type(QScrollArea) :: scroll
    type(QWidget) :: panel, row
    type(QVBoxLayout) :: column
    type(QHBoxLayout) :: hbox
    type(QGroupBox) :: group
    type(QFormLayout) :: form
    type(QPushButton) :: reset, reset_view_button
    type(QMenuBar) :: bar
    type(QMenu) :: menu
    type(QAction) :: action

    app%main = QMainWindow()
    call app%main%setWindowTitle('HALO')
    call app%main%closeEvent%connect(on_close)

    ! menu: File
    bar = app%main%menuBar()
    menu = bar%addMenu('&File')
    action = QAction('&Open Config...', app%main)
    call action%setShortcut('Ctrl+O')
    call action%triggered%connect(on_browse)
    call menu%addAction(action)
    action = QAction('Open &Trajectory...', app%main)
    call action%setShortcut('Ctrl+T')
    call action%triggered%connect(on_open_trajectory)
    call menu%addAction(action)
    call menu%addSeparator()
    action = QAction('&Quit', app%main)
    call action%setShortcut('Ctrl+Q')
    call action%triggered%connect(on_quit)
    call menu%addAction(action)

    ! the 3D view's window
    app%window = QWindow()
    app%container = QWidget_createWindowContainer(app%window)
    call app%container%setMinimumSize(320, 240)
    app%events = QtfEventFilter(app%window)
    call app%events%exposeEvent%connect(on_expose)
    call app%events%updateRequest%connect(on_update)
    call connect_view(app%view_link, app%events, app%proxy, app%window)  ! (mouse, keys, resizes, the surface)
    app%stack = QStackedWidget()
    app%startup_message = QLabel('Starting...')
    call app%startup_message%setAlignment(Qt_AlignCenter)
    call app%startup_message%setStyleSheet('background-color: black; color: #cccccc; font-size: 16pt;')
    call app%stack%addWidget(app%startup_message)
    call app%stack%addWidget(app%container)
    call app%main%setCentralWidget(app%stack)

    panel = QWidget()
    column = QVBoxLayout(panel)

    ! the solver
    group = QGroupBox('Solver')
    form = QFormLayout(group)
    row = QWidget()
    hbox = QHBoxLayout(row)
    call hbox%setContentsMargins(0, 0, 0, 0)
    app%config_edit = QLineEdit(app%config_file)
    call app%config_edit%setPlaceholderText('config file (.json)')
    call hbox%addWidget(app%config_edit)
    app%browse_button = QPushButton('Browse...')
    call app%browse_button%clicked%connect(on_browse)
    call hbox%addWidget(app%browse_button)
    call form%addRow('Config file', row)
    app%plot_iterations = QCheckBox('Plot each iteration')
    call app%plot_iterations%setChecked(app%iterations_at_start)
    call app%plot_iterations%setToolTip('Show the trajectory of each solver iteration as it solves '// &
                                        '(propagates all the segments at each iteration).')
    call form%addRow('', app%plot_iterations)
    row = QWidget()
    hbox = QHBoxLayout(row)
    call hbox%setContentsMargins(0, 0, 0, 0)
    app%run_button = QPushButton('Run solver')
    call app%run_button%clicked%connect(on_run)
    call hbox%addWidget(app%run_button)
    app%stop_button = QPushButton('Stop')
    call app%stop_button%setEnabled(.false.)
    call app%stop_button%clicked%connect(on_stop)
    call hbox%addWidget(app%stop_button)
    call form%addRow('', row)
    app%status_label = QLabel('Choose a config file and run the solver.')
    call app%status_label%setWordWrap(.true.)
    call form%addRow(app%status_label)
    call column%addWidget(group)

    ! the trajectory
    group = QGroupBox('Trajectory')
    form = QFormLayout(group)
    app%info_label = QLabel('None loaded')
    call app%info_label%setWordWrap(.true.)
    call form%addRow(app%info_label)
    app%frame_choice = QComboBox()
    call app%frame_choice%addItems(frame_names)
    call app%frame_choice%currentIndexChanged%connect(on_frame_choice)
    call form%addRow('Frame', app%frame_choice)
    reset_view_button = QPushButton('Reset view')
    call reset_view_button%clicked%connect(on_reset_view)
    call form%addRow('', reset_view_button)
    call column%addWidget(group)

    ! time
    group = QGroupBox('Time')
    form = QFormLayout(group)
    app%time_scale = QDoubleSpinBox()
    call app%time_scale%setRange(0.0_c_double, 1.0e9_c_double)
    call app%time_scale%setValue(10000.0_c_double)
    call app%time_scale%setSuffix(' x')
    call app%time_scale%valueChanged%connect(on_time_scale)
    call form%addRow('Time scale', app%time_scale)
    app%pause_time = QCheckBox('Pause')
    call app%pause_time%toggled%connect(on_pause_time)
    call form%addRow('', app%pause_time)
    reset = QPushButton('Reset time')
    call reset%clicked%connect(on_reset_time)
    call form%addRow('', reset)
    app%time_label = QLabel('')
    call form%addRow(app%time_label)
    app%time_slider = QSlider(Qt_Horizontal)
    call app%time_slider%setRange(0, slider_steps)
    call app%time_slider%setSingleStep(slider_steps / 100)
    call app%time_slider%setPageStep(slider_steps / 10)
    call app%time_slider%valueChanged%connect(on_time_slider)
    call app%time_slider%sliderPressed%connect(on_slider_pressed)
    call app%time_slider%sliderReleased%connect(on_slider_released)
    call form%addRow(app%time_slider)
    call column%addWidget(group)

    call column%addWidget(QLabel('<small>Left drag: rotate<br>Right drag: zoom<br>Middle drag: pan<br>' // &
                                 'Wheel: zoom</small>'))
    call column%addStretch()
    scroll = QScrollArea()
    call scroll%setWidgetResizable(.true.)
    call scroll%setWidget(panel)
    dock = QDockWidget('Controls', app%main)
    call dock%setWidget(scroll)
    call app%main%addDockWidget(Qt_LeftDockWidgetArea, dock)

    ! the solver's output
    app%log = QPlainTextEdit()
    call app%log%setReadOnly(.true.)
    call app%log%setMaximumBlockCount(20000)
    call app%log%setStyleSheet('font-family: Menlo, Monaco, monospace; font-size: 11pt;')
    dock = QDockWidget('Solver output', app%main)
    call dock%setWidget(app%log)
    call app%main%addDockWidget(Qt_BottomDockWidgetArea, dock)

    ! the iterations' trajectories (from the solver's thread), and its output (from the reader's)
    app%job%iteration = QSignalInt2(app%main)
    call app%job%iteration%connect(on_iteration)
    app%rjob%text = QSignalString(app%main)
    call app%rjob%text%connect(on_output)

    end subroutine build_ui
!*****************************************************************************************

!*****************************************************************************************
!>
!  Show the main window; the view opens once it is painted.

    subroutine start_up()

    call app%main%resize(1400, 900)
    call app%main%show()
    call QTimer_singleShot(0, on_startup)
    call app%options%start_quit_timer(app%main)

    end subroutine start_up
!*****************************************************************************************

!*****************************************************************************************
!  small utilities

    subroutine set_status(text)
    character(len=*),intent(in) :: text
    call app%status_label%setText(text)
    end subroutine set_status

    subroutine error_message(text)
    character(len=*),intent(in) :: text
    integer :: button
    call set_status(text)
    button = QMessageBox_warning('HALO', text, parent=app%main)
    end subroutine error_message

    function basename(path) result(name)
    character(len=*),intent(in) :: path
    character(len=:),allocatable :: name
    name = path(index(path, '/', back=.true.)+1:)
    end function basename

    function dir_of(path) result(dir)  !! the folder of `path` (the current one if none)
    character(len=*),intent(in) :: path
    character(len=:),allocatable :: dir
    integer :: i
    i = index(path, '/', back=.true.)
    if (i > 0) then
        dir = path(1:i-1)
    else
        dir = QDir_currentPath()
    end if
    end function dir_of

    function relative_path(path) result(rel)  !! `path` relative to the current folder, if in it
    character(len=*),intent(in) :: path
    character(len=:),allocatable :: rel, cwd
    cwd = QDir_currentPath()//'/'
    if (len(path) > len(cwd)) then
        if (path(1:len(cwd)) == cwd) then
            rel = path(len(cwd)+1:)
            return
        end if
    end if
    rel = path
    end function relative_path

    function int_text(i) result(s)
    integer,intent(in) :: i
    character(len=:),allocatable :: s
    character(len=32) :: buf
    write(buf,'(I0)') i
    s = trim(buf)
    end function int_text

    function real_text(x, fmt) result(s)
    real(wp),intent(in) :: x
    character(len=*),intent(in) :: fmt
    character(len=:),allocatable :: s
    character(len=64) :: buf
    write(buf, fmt) x
    s = trim(adjustl(buf))
    if (len(s) > 0) then
        if (s(len(s):len(s)) == '.') s = s(1:len(s)-1)  ! (F0.0: "3210.")
        if (s(1:1) == '.') s = '0'//s
    end if
    end function real_text
!*****************************************************************************************

    end module halo_gui_module
!*****************************************************************************************

!*****************************************************************************************
!>
!  The HALO GUI program (see [[halo_gui_module]]).
!
!  Arguments: `[--dir DIR] [config.json] [--run] [--iterations] [--traj traj_file.json]`.
!  `--dir` sets the current directory (the halo root, which the config files'
!  paths are relative to; give it first: the other paths are relative to it).
!  The config file defaults to `examples/example_sparse.json` (if it exists);
!  `--run` runs the solver on it at startup; `--iterations` checks "Plot each
!  iteration"; `--traj` plots a trajectory file at startup.
!  For testing: `--capture file.png` saves frame `--capture-frame n` (default 120),
!  and `--quit-after sec` closes the program.

    program halo_gui

    use qt
    use halo_gui_module

    implicit none

    type(QApplication) :: qapp
    character(len=4096) :: arg
    integer :: i, ierr
    logical :: exists

    qapp = QApplication()

    app%config_file = ''
    i = 1
    do while (i <= command_argument_count())
        if (app%options%take(i)) cycle  ! (--capture file, --capture-frame n, --quit-after sec)
        call get_command_argument(i, arg)
        select case (arg)
        case ('--dir')  ! (first: the other arguments' paths are relative to it)
            call get_command_argument(i + 1, arg)
            if (.not. change_directory(trim(arg))) then
                write(*,'(A)') 'halo_gui: cannot change to directory '//trim(arg)
                stop 1
            end if
            i = i + 1
        case ('--run')
            app%run_at_start = .true.
        case ('--iterations')
            app%iterations_at_start = .true.
        case ('--traj')
            call get_command_argument(i + 1, arg)
            app%initial_traj = trim(arg)
            i = i + 1
        case default
            app%config_file = trim(arg)
        end select
        i = i + 1
    end do
    if (len(app%config_file) == 0) then
        inquire(file='examples/example_sparse.json', exist=exists)
        if (exists) app%config_file = 'examples/example_sparse.json'
    end if

    call build_scene()
    call build_ui()
    call start_up()
    ierr = qapp%exec()

    end program halo_gui
!*****************************************************************************************
