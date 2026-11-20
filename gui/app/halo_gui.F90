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
!  * File > Open Trajectory plots an existing `traj_*.json` file without solving (and
!    its `eclipse_*.json` file, if there is one).
!  * Plots (tabs next to the solver output): the constraint defects at the segments'
!    interfaces (|dr| and |dv| where segments 2k-1 and 2k end; updated at each
!    iteration while solving), and the eclipses (`phi` < 0: in the Earth's shadow)
!    over time, with the current time.
!  * The trajectory is drawn in the Moon-centered Earth-Moon rotating frame, the
!    Moon-centered inertial (J2000) frame, or the Earth-centered Sun-Earth rotating
!    frame (with the Earth's umbra and penumbra, along +x, away from the Sun), as in
!    `python/plot_utilities.py`. It is colored ("Color by") by the distance from the
!    Moon, by eclipses (in sunlight or in the Earth's shadow: `phi_se` from the
!    `eclipse_*.json` file written with the solution), or cyan. With a marker for the current state (following it
!    in time), and green and red markers at each segment's start and end (where its
!    propagation starts and ends: during the iterations, the red ones show the gaps
!    the solver is closing).
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
    public :: COLOR_ALTITUDE, COLOR_ECLIPSE, COLOR_NONE  ! (for --color)

    integer,parameter :: slider_steps = 1000  !! the time slider's range
    real(wp),parameter :: moon_radius = 1737.4_wp  !! km
    real(wp),parameter :: day = 86400.0_wp  !! sec

    real(wp),parameter :: earth_radius = 6378.137_wp  !! km
    real(wp),parameter :: sun_radius = 696000.0_wp  !! km
    real(wp),parameter :: au = 149597870.7_wp  !! km (the shadow's geometry: at 1 au)
    real(wp),parameter :: shadow_length = 6.0e5_wp  !! km (the umbra and penumbra drawn: beyond the Moon)

    !> the frames the trajectory can be shown in
    integer,parameter :: ROTATING = 1, INERTIAL = 2, SUN_EARTH = 3
    character(len=*),parameter :: frame_names(3) = [character(len=32) :: &
                                    'Earth-Moon rotating (EMR)', 'Inertial (J2000)', 'Sun-Earth rotating (SER)']
    character(len=*),parameter :: frame_centers(3) = [character(len=16) :: &
                                    'Moon-centered', 'Moon-centered', 'Earth-centered']

    !> how the trajectory is colored ("Color by": the combo box's index)
    integer,parameter :: COLOR_ALTITUDE = 0, COLOR_ECLIPSE = 1, COLOR_NONE = 2
    character(len=*),parameter :: color_names(0:2) = [character(len=36) :: &
                                    'Altitude (distance from the Moon)', 'Eclipses (Earth''s shadow)', 'None (cyan)']
    real(wp),parameter :: lit_color(3) = [1.0_wp, 0.85_wp, 0.25_wp]  !! in sunlight (yellow)
    real(wp),parameter :: shadow_color_edge(3) = [0.3_wp, 0.8_wp, 1.0_wp]  !! just in the shadow (cyan)
    real(wp),parameter :: shadow_color_deep(3) = [0.5_wp, 0.25_wp, 1.0_wp]  !! deepest in it (violet)
    real(wp),parameter :: no_data_color(3) = [0.6_wp, 0.6_wp, 0.6_wp]  !! (no eclipse data)

    !> colormap for the distance from the Moon: near (red) to far (blue)
    real(wp),parameter :: cmap(3,5) = reshape([ 1.00_wp, 0.25_wp, 0.20_wp, &
                                                1.00_wp, 0.75_wp, 0.15_wp, &
                                                0.55_wp, 1.00_wp, 0.35_wp, &
                                                0.20_wp, 0.85_wp, 1.00_wp, &
                                                0.35_wp, 0.45_wp, 1.00_wp ], [3,5])

    real(wp),parameter :: line_color(3) = [0.2_wp, 0.9_wp, 1.0_wp]  !! the trajectory with "Color by" None (cyan)
    real(wp),parameter :: start_color(3) = [0.2_wp, 1.0_wp, 0.3_wp]  !! segment start markers (green)
    real(wp),parameter :: end_color(3) = [1.0_wp, 0.2_wp, 0.2_wp]  !! segment end markers (red)
    !> the plots' colors (dark): [r, g, b], 0 to 255
    integer,parameter :: plot_background(3) = [30, 31, 34]
    integer,parameter :: plot_grid(3)       = [62, 64, 70]
    integer,parameter :: plot_frame(3)      = [110, 112, 120]
    integer,parameter :: plot_text(3)       = [190, 192, 198]  !! tick and axis labels, legend
    integer,parameter :: plot_dim_text(3)   = [140, 142, 150]  !! messages (no data)
    integer,parameter :: plot_title(3)      = [235, 236, 240]
    integer,parameter :: plot_dr(3)         = [90, 170, 255]   !! |dr| (blue)
    integer,parameter :: plot_dv(3)         = [255, 110, 110]  !! |dv| (red)
    integer,parameter :: plot_phi(3)        = [175, 145, 255]  !! eclipses' phi (violet)
    integer,parameter :: plot_now(3)        = [255, 180, 50]   !! the current time (orange)

    !> the window's colors (dark, like the plots): [r, g, b], 0 to 255
    integer,parameter :: ui_window(3)    = [30, 31, 34]     !! panels' background
    integer,parameter :: ui_base(3)      = [43, 45, 49]     !! fields' background
    integer,parameter :: ui_button(3)    = [52, 54, 60]
    integer,parameter :: ui_text(3)      = [220, 221, 225]
    integer,parameter :: ui_disabled(3)  = [100, 103, 110]
    integer,parameter :: ui_accent(3)    = [61, 139, 253]   !! (blue) the selection, focus, the main button
    character(len=*),parameter :: ui_dim = '#8b8e96'        !! secondary text (names, hints, status)
    !> the controls' look (Fusion style, with the colors above)
    character(len=*),parameter :: ui_style_sheet = &
        'QMainWindow, QDockWidget > QWidget { background: #1e1f22; }'// &
        'QDockWidget { color: #b4b6bc; }'// &
        'QDockWidget::title { background: #18191b; padding: 6px 8px; }'// &
        'QScrollArea { border: none; }'// &
        'QLabel#section { color: '//ui_dim//'; font-weight: bold; padding-top: 14px; }'// &
        'QLabel#dim { color: '//ui_dim//'; }'// &
        'QLineEdit, QComboBox, QDoubleSpinBox { background: #2b2d31; border: 1px solid #3c3f45;'// &
        ' border-radius: 6px; padding: 4px 8px; min-height: 18px; }'// &
        'QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus { border-color: #3d8bfd; }'// &
        'QLineEdit:disabled, QComboBox:disabled { color: #64676e; }'// &
        'QComboBox QAbstractItemView { background: #2b2d31; border: 1px solid #3c3f45;'// &
        ' selection-background-color: #3d8bfd; outline: none; }'// &
        'QDoubleSpinBox::up-button, QDoubleSpinBox::down-button { width: 0; border: none; }'// &
        'QPushButton, QToolButton { background: #34363c; border: 1px solid #3c3f45; border-radius: 6px;'// &
        ' padding: 5px 12px; min-height: 18px; }'// &
        'QPushButton:hover, QToolButton:hover { background: #3e4148; }'// &
        'QPushButton:pressed, QToolButton:pressed { background: #2b2d31; }'// &
        'QPushButton:disabled, QToolButton:disabled { background: #26282c; color: #5c5f66;'// &
        ' border-color: #303237; }'// &
        'QPushButton#primary { background: #3d8bfd; border-color: #3d8bfd; color: white; font-weight: bold; }'// &
        'QPushButton#primary:hover { background: #5a9cff; }'// &
        'QPushButton#primary:disabled { background: #2a3d5c; border-color: #2a3d5c; color: #8b8e96; }'// &
        'QToolButton { padding: 5px 8px; min-width: 18px; }'// &
        'QSlider::groove:horizontal { height: 4px; background: #3c3f45; border-radius: 2px; }'// &
        'QSlider::sub-page:horizontal { background: #3d8bfd; border-radius: 2px; }'// &
        'QSlider::handle:horizontal { background: #dcdde1; width: 14px; height: 14px; margin: -5px 0;'// &
        ' border-radius: 7px; }'// &
        'QCheckBox { spacing: 8px; }'// &
        'QScrollBar:vertical { background: transparent; width: 8px; margin: 2px; }'// &
        'QScrollBar::handle:vertical { background: #3c3f45; border-radius: 3px; min-height: 30px; }'// &
        'QScrollBar::handle:vertical:hover { background: #4a4e55; }'// &
        'QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }'// &
        'QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical { background: none; }'// &
        'QTabWidget::pane { border: none; }'// &
        'QTabBar::tab { background: transparent; color: '//ui_dim//'; padding: 6px 14px;'// &
        ' border: none; border-bottom: 2px solid transparent; }'// &
        'QTabBar::tab:selected { color: #dcdde1; border-bottom: 2px solid #3d8bfd; }'// &
        'QTabBar::tab:hover { color: #dcdde1; }'// &
        'QPlainTextEdit { background: #18191b; border: none; color: #c8cad0; }'// &
        'QToolTip { background: #2b2d31; color: #dcdde1; border: 1px solid #3c3f45; }'
    character(len=*),parameter :: play_symbol = '▶'   !! (on the play/pause button)
    character(len=*),parameter :: pause_symbol = '❚❚'
    character(len=*),parameter :: reset_symbol = '↺'

    real(wp),parameter :: segment_marker_size = 8.0_wp  !! (window units: the largest they get, with the
                                                        !! camera within the markers' extent; smaller
                                                        !! farther away)
    real(wp),parameter :: state_marker_size = 14.0_wp  !! the current state's marker

    !> The trajectory in one frame, and what draws it.
    type :: frame_t
        type(Trajectory) :: traj
        type(Trajectory) :: starts, ends  !! each segment's start and end points
        type(DrawableTrajectory) :: drawable
        type(CurveArtist) :: curve
        type(MarkerArtist) :: start_markers, end_markers
        real(wp) :: radius = 1.0e5_wp  !! the trajectory's farthest distance from the frame's center (km)
    end type frame_t

    !> One iteration's trajectory, made by the worker thread and handed to the
    !> GUI thread (which deletes it) through `solver_job%iteration`.
    type :: iteration_points
        integer :: iter = 0
        real(wp),dimension(:),allocatable :: t  !! (n) ET (sec)
        real(wp),dimension(:,:),allocatable :: rot, inr, se  !! (3,n) km
        real(wp),dimension(:,:),allocatable :: seg_t  !! (2,nseg) each segment's start and end times
        real(wp),dimension(:,:,:),allocatable :: seg_rot, seg_inr, seg_se  !! (3,2,nseg) and positions
        real(wp),dimension(:,:),allocatable :: defects  !! (2,nseg/2) |dr| (km), |dv| (km/s)
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
        type(QComboBox) :: color_mode  !! how the trajectory is colored (COLOR_*)
        type(QCheckBox) :: show_segment_markers  !! the segments' start and end markers
        type(QDoubleSpinBox) :: time_scale
        type(QToolButton) :: pause_time  !! play/pause (checked: paused)
        type(QLabel) :: time_label
        type(QSlider) :: time_slider
        logical :: slider_held = .false.
        logical :: slider_setting = .false.
        type(QPlainTextEdit) :: log
        type(QTabWidget) :: bottom_tabs  !! the solver output, and the plots
        type(QtfWidget) :: defects_plot, eclipse_plot
        ! the plots' data
        real(wp),dimension(:,:),allocatable :: defects  !! (2,npairs) |dr| (km), |dv| (km/s) at the interfaces
        character(len=:),allocatable :: defects_title
        logical :: has_eclipses = .false.  !! eclipse data for the trajectory shown
        real(wp),dimension(:),allocatable :: ecl_t, ecl_phi  !! its times (ET sec, increasing) and phi
        real(wp) :: ecl_phi_min = 0.0_wp  !! the most negative phi (0: no eclipses)
        real(wp) :: rmin = 0.0_wp, rmax = 0.0_wp  !! the trajectory's distances from the Moon (km)
        ! a plot being painted: its area (pixels), its data ranges, its legend's entries
        real(wp) :: plot_x0 = 0.0_wp, plot_x1 = 1.0_wp, plot_y0 = 0.0_wp, plot_y1 = 1.0_wp
        real(wp) :: plot_xmin = 0.0_wp, plot_xmax = 1.0_wp, plot_ymin = 0.0_wp, plot_ymax = 1.0_wp
        integer :: plot_legend_row = 0
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
        type(ReferenceFrame) :: root  !! the frame's center (the Moon, or the Earth)
        type(Sphere) :: moon, earth
        type(Polyhedron) :: umbra, penumbra  !! the Earth's shadow (Sun-Earth rotating frame)
        type(ReferenceFrame) :: body  !! follows the trajectory (the current state)
        type(Trajectory) :: state_point  !! one point, at the body's origin
        type(DrawableTrajectory) :: state_drawable
        type(MarkerArtist) :: state_marker
        type(frame_t) :: frames(3)  !! rotating, inertial, Sun-Earth rotating
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
        integer :: frame_at_start = ROTATING  !! --frame n
        integer :: color_at_start = COLOR_ALTITUDE  !! --color altitude|eclipse|none
        integer :: tab_at_start = 0  !! --tab n (the bottom tabs: 0 output, 1 defects, 2 eclipses)
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

    app%root = ReferenceFrame('Center')
    call app%root%show_decorations(.false.)

    app%moon = Sphere('Moon', radius=moon_radius, color=[0.6_wp, 0.6_wp, 0.6_wp])
    call app%moon%set_texture_map('lroc_color_2k.jpg')  ! (vkframes' data/, on VKFRAMES_DATA_PATH)
    call app%moon%show_axes(NO_AXES)
    call app%moon%show_axes_labels(NO_AXES)
    if (.not. app%root%add_child(app%moon)) write(*,'(A)') 'halo_gui: add_child failed'

    ! the Earth and its shadow (Sun-Earth rotating frame: the Sun along -x)
    app%earth = Sphere('Earth', radius=earth_radius, color=[0.3_wp, 0.5_wp, 0.9_wp])
    call app%earth%set_texture_map('2k_earth_daymap.jpg')
    call app%earth%show_axes(NO_AXES)
    call app%earth%show_axes_labels(NO_AXES)
    if (.not. app%root%add_child(app%earth)) write(*,'(A)') 'halo_gui: add_child failed'
    call shadow_surface(app%umbra, 'Umbra', &
                        earth_radius - shadow_length * (sun_radius - earth_radius) / au, &
                        [0.55_wp, 0.6_wp, 0.8_wp], 0.35_wp)
    call shadow_surface(app%penumbra, 'Penumbra', &
                        earth_radius + shadow_length * (sun_radius + earth_radius) / au, &
                        [0.55_wp, 0.6_wp, 0.8_wp], 0.12_wp)
    if (.not. app%root%add_child(app%umbra)) write(*,'(A)') 'halo_gui: add_child failed'
    if (.not. app%root%add_child(app%penumbra)) write(*,'(A)') 'halo_gui: add_child failed'

    do i = 1, size(app%frames)
        associate (f => app%frames(i))
        f%traj = Trajectory()
        f%curve = CurveArtist(f%traj, color=[1.0_wp, 1.0_wp, 1.0_wp], width=2.0_wp)
        call f%curve%set_point_colors(.true.)  ! (the colors: see trajectory_colors)
        ! the segments' starts and ends: a marker at each point
        f%starts = Trajectory()
        f%ends = Trajectory()
        call segment_markers(f%start_markers, f%starts, start_color)
        call segment_markers(f%end_markers, f%ends, end_color)
        f%drawable = DrawableTrajectory(trim(frame_names(i)))
        call f%drawable%show_decorations(.false.)
        call f%drawable%add_artist(f%curve)
        call f%drawable%add_artist(f%start_markers)
        call f%drawable%add_artist(f%end_markers)
        if (.not. app%root%add_child(f%drawable)) write(*,'(A)') 'halo_gui: add_child failed'
        end associate
    end do

    ! the current state: a marker at the origin of a frame following the trajectory
    app%body = ReferenceFrame('HALO')
    call app%body%show_axes(NO_AXES)
    call app%body%show_axes_labels(NO_AXES)
    app%state_point = Trajectory()
    call app%state_point%add([0.0_wp], pos=reshape([0.0_wp, 0.0_wp, 0.0_wp], [3,1]))
    app%state_marker = MarkerArtist(app%state_point, markers=MARKER_START, size=state_marker_size, &
                                    color=[1.0_wp, 1.0_wp, 1.0_wp])
    call app%state_marker%set_marker_color(MARKER_START, [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp])
    call app%state_marker%set_marker_shape(MARKER_TARGET)
    call app%state_marker%set_auto_attenuate(.false.)  ! (one point: no size to attenuate by)
    app%state_drawable = DrawableTrajectory('state')
    call app%state_drawable%show_decorations(.false.)
    call app%state_drawable%add_artist(app%state_marker)
    if (.not. app%body%add_child(app%state_drawable)) write(*,'(A)') 'halo_gui: add_child failed'
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
!  A MarkerArtist drawing a marker of `color` at each point of `traj`.

    subroutine segment_markers(artist, traj, color)

    type(MarkerArtist),intent(out) :: artist
    type(Trajectory),intent(in),target :: traj
    real(wp),dimension(3),intent(in) :: color

    artist = MarkerArtist(traj, markers=ior(ior(MARKER_START, MARKER_INTERMEDIATE), MARKER_END), &
                          size=segment_marker_size, color=color)
    call artist%set_intermediate_type(INTERMEDIATE_DATA)
    call artist%set_intermediate_spacing(1.0_wp)  ! (every point)
    call artist%set_auto_attenuate(.true.)  ! (smaller as the camera moves away: size * min(1, sqrt(R/d)),
                                            !  R the markers' extent, d the distance)
    call artist%set_marker_color(MARKER_START, [color, 1.0_wp])
    call artist%set_marker_color(MARKER_INTERMEDIATE, [color, 1.0_wp])
    call artist%set_marker_color(MARKER_END, [color, 1.0_wp])

    end subroutine segment_markers
!*****************************************************************************************

!*****************************************************************************************
!>
!  A part of the Earth's shadow, along +x (away from the Sun) in the Sun-Earth
!  rotating frame: a translucent surface from the Earth (its radius) to
!  `shadow_length` (radius `r_end`): the umbra narrows (to its apex, about 1.38
!  million km away), the penumbra widens. (Cones tangent to the Sun and the
!  Earth, with them 1 au apart.)

    subroutine shadow_surface(poly, name, r_end, color, opacity)

    type(Polyhedron),intent(out) :: poly
    character(len=*),intent(in) :: name
    real(wp),intent(in) :: r_end  !! the radius at `shadow_length` (km)
    real(wp),dimension(3),intent(in) :: color
    real(wp),intent(in) :: opacity

    integer,parameter :: nu = 72  !! points around
    real(wp),dimension(:,:,:),allocatable :: points
    real(wp) :: angle
    integer :: i

    allocate(points(3,nu,2))
    do i = 1, nu
        angle = 2.0_wp * acos(-1.0_wp) * real(i - 1, wp) / nu
        points(:,i,1) = [0.0_wp,        earth_radius * cos(angle), earth_radius * sin(angle)]
        points(:,i,2) = [shadow_length, r_end * cos(angle),        r_end * sin(angle)]
    end do
    poly = Polyhedron(name, color=color)
    call poly%show_decorations(.false.)
    call poly%set_grid(points, closed_u=.true., lines=[0, 0])
    call poly%set_draw_mode(POLY_FACES)
    call poly%set_lit(.false.)
    call poly%set_face_color(color)
    call poly%set_opacity(opacity)

    end subroutine shadow_surface
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
    real(wp),dimension(:),allocatable :: et, xr, yr, zr, xi, yi, zi, xs, ys, zs, vx, vy, vz
    real(wp),dimension(:),allocatable :: t_all
    real(wp),dimension(:,:),allocatable :: rot_all, inr_all, se_all
    real(wp),dimension(:,:),allocatable :: seg_t  !! (2,nseg) each segment's start and end times
    real(wp),dimension(:,:,:),allocatable :: seg_rot, seg_inr, seg_se  !! (3,2,nseg) and positions
    real(wp),dimension(:,:),allocatable :: seg_v_end  !! (3,nseg) each segment's final velocity (inertial)
    real(wp),dimension(:,:),allocatable :: defects  !! (2,nseg/2)
    integer :: nsegs, iseg, n, ntot, k, kstart, kend, kstep, nseg_ok
    logical :: found, status_ok, has_se, has_v
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
    allocate(t_all(0), rot_all(3,0), inr_all(3,0), se_all(3,0))
    allocate(seg_t(2,nsegs), seg_rot(3,2,nsegs), seg_inr(3,2,nsegs), seg_se(3,2,nsegs), seg_v_end(3,nsegs))
    nseg_ok = 0
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
        ! the Sun-Earth rotating frame (zero if not in the file):
        has_se = .true.
        call json%get(seg, 'x_se_rotating', xs, found);  has_se = has_se .and. found
        call json%get(seg, 'y_se_rotating', ys, found);  has_se = has_se .and. found
        call json%get(seg, 'z_se_rotating', zs, found);  has_se = has_se .and. found
        if (has_se) has_se = min(size(xs), size(ys), size(zs)) >= n
        if (.not. has_se) then
            if (allocated(xs)) deallocate(xs, ys, zs)
            allocate(xs(n), ys(n), zs(n))
            xs = 0.0_wp;  ys = 0.0_wp;  zs = 0.0_wp
        end if
        ! the final velocity (inertial: for the defects)
        has_v = .true.
        call json%get(seg, 'vx_inertial', vx, found);  has_v = has_v .and. found
        call json%get(seg, 'vy_inertial', vy, found);  has_v = has_v .and. found
        call json%get(seg, 'vz_inertial', vz, found);  has_v = has_v .and. found
        if (has_v) has_v = min(size(vx), size(vy), size(vz)) >= n
        if (n == 0) cycle
        ! the segment's start and end (as propagated: the file's order)
        nseg_ok = nseg_ok + 1
        seg_t(1,nseg_ok) = et(1)
        seg_t(2,nseg_ok) = et(n)
        seg_rot(1,1,nseg_ok) = xr(1);  seg_rot(2,1,nseg_ok) = yr(1);  seg_rot(3,1,nseg_ok) = zr(1)
        seg_rot(1,2,nseg_ok) = xr(n);  seg_rot(2,2,nseg_ok) = yr(n);  seg_rot(3,2,nseg_ok) = zr(n)
        seg_inr(1,1,nseg_ok) = xi(1);  seg_inr(2,1,nseg_ok) = yi(1);  seg_inr(3,1,nseg_ok) = zi(1)
        seg_inr(1,2,nseg_ok) = xi(n);  seg_inr(2,2,nseg_ok) = yi(n);  seg_inr(3,2,nseg_ok) = zi(n)
        seg_se(1,1,nseg_ok) = xs(1);  seg_se(2,1,nseg_ok) = ys(1);  seg_se(3,1,nseg_ok) = zs(1)
        seg_se(1,2,nseg_ok) = xs(n);  seg_se(2,2,nseg_ok) = ys(n);  seg_se(3,2,nseg_ok) = zs(n)
        if (has_v) then
            seg_v_end(:,nseg_ok) = [vx(n), vy(n), vz(n)]
        else
            seg_v_end(:,nseg_ok) = 0.0_wp
        end if
        ! the joined trajectory:
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
            se_all(1,ntot) = xs(k);  se_all(2,ntot) = ys(k);  se_all(3,ntot) = zs(k)
        end do
    end do
    call json%destroy(root)

    if (ntot < 2) then
        call error_message('No trajectory points in '//file)
        return
    end if

    ! the defects: segments 2k-1 and 2k end at the same point (when solved).
    ! (|dv| from the inertial velocities: the rotating frame's differs by w x dr, negligible)
    allocate(defects(2, nseg_ok/2))
    do k = 1, nseg_ok/2
        defects(1,k) = norm2(seg_inr(:,2,2*k-1) - seg_inr(:,2,2*k))
        defects(2,k) = norm2(seg_v_end(:,2*k-1) - seg_v_end(:,2*k))
    end do

    ! the eclipses, from the solution's eclipse file (if there is one):
    call load_eclipses(eclipse_file_for(file))

    call plot_points(t_all(1:ntot), rot_all(:,1:ntot), inr_all(:,1:ntot), se_all(:,1:ntot), &
                     seg_t(:,1:nseg_ok), seg_rot(:,:,1:nseg_ok), seg_inr(:,:,1:nseg_ok), seg_se(:,:,1:nseg_ok), &
                     defects, basename(file), int_text(nsegs)//' segments', keep_view)
    call app%main%setWindowTitle('HALO - '//basename(file))
    ok = .true.

    contains

        subroutine grow(nmax)  !! make room for nmax points (copies on the heap)
            integer,intent(in) :: nmax
            integer :: m
            real(wp),dimension(:),allocatable :: t_new
            real(wp),dimension(:,:),allocatable :: rot_new, inr_new, se_new
            if (size(t_all) >= nmax) return
            m = max(nmax, 2 * size(t_all))
            allocate(t_new(m), rot_new(3,m), inr_new(3,m), se_new(3,m))
            if (ntot > 0) then
                t_new(1:ntot) = t_all(1:ntot)
                rot_new(:,1:ntot) = rot_all(:,1:ntot)
                inr_new(:,1:ntot) = inr_all(:,1:ntot)
                se_new(:,1:ntot) = se_all(:,1:ntot)
            end if
            call move_alloc(t_new, t_all)
            call move_alloc(rot_new, rot_all)
            call move_alloc(inr_new, inr_all)
            call move_alloc(se_new, se_all)
        end subroutine grow

    end function load_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  The eclipse file written with a trajectory file: `eclipse_*.json` for
!  `traj_*.json`, in the same folder ('' if the name isn't like that).

    function eclipse_file_for(traj_file) result(file)

    character(len=*),intent(in) :: traj_file
    character(len=:),allocatable :: file

    integer :: i

    i = index(traj_file, '/', back=.true.)
    file = ''
    if (len(traj_file) < i + 5) return
    if (traj_file(i+1:i+5) /= 'traj_') return
    file = traj_file(1:i)//'eclipse_'//traj_file(i+6:)

    end function eclipse_file_for
!*****************************************************************************************

!*****************************************************************************************
!>
!  Read the eclipse data (`et_se` and `phi_se`: < 0 in the Earth's shadow) from
!  an `eclipse_*.json` file, for the coloring and the plot. None if the file
!  isn't there (or isn't one).

    subroutine load_eclipses(file)

    character(len=*),intent(in) :: file

    type(json_core) :: json
    type(json_value),pointer :: root
    real(wp),dimension(:),allocatable :: et, phi
    logical :: exists, found1, found2
    integer :: i, n, m

    call clear_eclipses()
    if (len(file) == 0) return
    inquire(file=file, exist=exists)
    if (.not. exists) return
    call json%initialize()
    call json%load(file=file, p=root)
    if (json%failed()) then
        call json%destroy()
        return
    end if
    call json%get(root, 'et_se', et, found1)
    call json%get(root, 'phi_se', phi, found2)
    call json%destroy(root)
    if (.not. (found1 .and. found2)) return
    n = min(size(et), size(phi))
    if (n < 2) return

    ! increasing times (the segments' shared ends only once):
    allocate(app%ecl_t(n), app%ecl_phi(n))
    m = 0
    do i = 1, n
        if (m > 0) then
            if (et(i) <= app%ecl_t(m)) cycle
        end if
        m = m + 1
        app%ecl_t(m) = et(i)
        app%ecl_phi(m) = phi(i)
    end do
    if (m < 2) then
        call clear_eclipses()
        return
    end if
    app%ecl_t = app%ecl_t(1:m)
    app%ecl_phi = app%ecl_phi(1:m)
    app%ecl_phi_min = min(0.0_wp, minval(app%ecl_phi))
    app%has_eclipses = .true.

    end subroutine load_eclipses
!*****************************************************************************************

!*****************************************************************************************
    subroutine clear_eclipses()  !! no eclipse data
    if (allocated(app%ecl_t)) deallocate(app%ecl_t)
    if (allocated(app%ecl_phi)) deallocate(app%ecl_phi)
    app%has_eclipses = .false.
    app%ecl_phi_min = 0.0_wp
    end subroutine clear_eclipses
!*****************************************************************************************

!*****************************************************************************************
!>
!  The eclipse data's phi at time `t` (linear between its times; held outside them).

    real(wp) function eclipse_phi(t) result(phi)

    real(wp),intent(in) :: t  !! ET (sec)

    integer :: lo, hi, mid
    real(wp) :: w

    associate (et => app%ecl_t, p => app%ecl_phi)
        if (t <= et(1)) then
            phi = p(1)
        else if (t >= et(size(et))) then
            phi = p(size(p))
        else
            lo = 1
            hi = size(et)
            do while (hi - lo > 1)  ! et(lo) <= t < et(hi)
                mid = (lo + hi) / 2
                if (et(mid) <= t) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            w = (t - et(lo)) / (et(hi) - et(lo))
            phi = (1.0_wp - w) * p(lo) + w * p(hi)
        end if
    end associate

    end function eclipse_phi
!*****************************************************************************************

!*****************************************************************************************
!>
!  Plot a trajectory (in both frames): its points put into the vkframes
!  trajectories, colored by the distance from the Moon.

    subroutine plot_points(t, rot, inr, se, seg_t, seg_rot, seg_inr, seg_se, defects, title, details, keep_view)

    real(wp),dimension(:),intent(in) :: t  !! (n) ET (sec)
    real(wp),dimension(:,:),intent(in) :: rot  !! (3,n) Moon-centered rotating frame (km)
    real(wp),dimension(:,:),intent(in) :: inr  !! (3,n) Moon-centered inertial frame (km)
    real(wp),dimension(:,:),intent(in) :: se  !! (3,n) Earth-centered Sun-Earth rotating frame (km)
    real(wp),dimension(:,:),intent(in) :: seg_t  !! (2,nseg) each segment's start and end times (sec)
    real(wp),dimension(:,:,:),intent(in) :: seg_rot  !! (3,2,nseg) their positions, rotating frame (km)
    real(wp),dimension(:,:,:),intent(in) :: seg_inr  !! (3,2,nseg) their positions, inertial frame (km)
    real(wp),dimension(:,:,:),intent(in) :: seg_se  !! (3,2,nseg) their positions, Sun-Earth rotating frame (km)
    real(wp),dimension(:,:),intent(in) :: defects  !! (2,nseg/2) |dr| (km), |dv| (km/s) at the interfaces
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
    app%rmin = rmin
    app%rmax = rmax
    call trajectory_colors(t, r, colors)

    call app%frames(ROTATING)%traj%clear()
    call app%frames(ROTATING)%traj%add(t, pos=rot, color=colors)
    call app%frames(INERTIAL)%traj%clear()
    call app%frames(INERTIAL)%traj%add(t, pos=inr, color=colors)
    call app%frames(SUN_EARTH)%traj%clear()
    call app%frames(SUN_EARTH)%traj%add(t, pos=se, color=colors)
    call set_bounds(app%frames(ROTATING), rot, moon_radius)
    call set_bounds(app%frames(INERTIAL), inr, moon_radius)
    call set_bounds(app%frames(SUN_EARTH), se, earth_radius)
    call set_segment_points(app%frames(ROTATING), seg_t, seg_rot)
    call set_segment_points(app%frames(INERTIAL), seg_t, seg_inr)
    call set_segment_points(app%frames(SUN_EARTH), seg_t, seg_se)

    ! the plots:
    app%defects = defects
    app%defects_title = title
    call app%defects_plot%update()
    call app%eclipse_plot%update()

    app%t0 = t(1)
    app%tf = t(n)
    app%loaded = .true.

    call app%hud_title%set_text(title)
    call show_legend()
    call app%info_label%setText('<b>'//breakable(title)//'</b>'// &
                                '<table cellspacing="0" cellpadding="1">'// &
                                info_row('Points', int_text(n)//', '//details)// &
                                info_row('Duration', real_text((app%tf - app%t0)/day, '(F0.2)')//' days')// &
                                info_row('|r|', real_text(rmin, '(F0.0)')//' to '// &
                                                real_text(rmax, '(F0.0)')//' km')// &
                                max_defect_text()//'</table>')
    call show_frame(app%current_frame)

    if (keep_view) then
        ! the same time (within the new range), the same camera:
        if (app%proxy%is_open()) then
            time = max(app%t0, min(app%tf, app%proxy%time()))
            call app%proxy%set_time_range(app%t0, app%tf)
            call app%proxy%set_time(time)
        end if
    else
        call update_time_range()  ! (at the start; the time scale and pause as they are)
        call reset_view()
    end if

    contains

        subroutine set_segment_points(f, st, sp)  !! the segments' start and end markers' points
            type(frame_t),intent(inout) :: f
            real(wp),dimension(:,:),intent(in) :: st  !! (2,nseg)
            real(wp),dimension(:,:,:),intent(in) :: sp  !! (3,2,nseg)
            real(wp),dimension(:,:),allocatable :: pos
            integer :: k, m
            m = size(st,2)
            call f%starts%clear()
            call f%ends%clear()
            if (m == 0) return
            allocate(pos(3,m))
            do k = 1, m
                pos(:,k) = sp(:,1,k)
            end do
            call f%starts%add(st(1,:), pos=pos)
            do k = 1, m
                pos(:,k) = sp(:,2,k)
            end do
            call f%ends%add(st(2,:), pos=pos)
        end subroutine set_segment_points

        subroutine set_bounds(f, pos, body_radius)  !! the trajectory's farthest distance from the center
            type(frame_t),intent(inout) :: f
            real(wp),dimension(:,:),intent(in) :: pos
            real(wp),intent(in) :: body_radius  !! the central body's (km)
            integer :: i
            f%radius = body_radius
            do i = 1, size(pos,2)
                f%radius = max(f%radius, norm2(pos(:,i)))
            end do
        end subroutine set_bounds

    end subroutine plot_points
!*****************************************************************************************

!*****************************************************************************************
!>
!  The largest defects, for the panel ('' if none).

    function max_defect_text() result(text)

    character(len=:),allocatable :: text

    text = ''
    if (.not. allocated(app%defects)) return
    if (size(app%defects,2) == 0) return
    text = info_row('Max defect', real_text(maxval(app%defects(1,:)), '(ES9.2)')//' km, '// &
                                  real_text(maxval(app%defects(2,:)), '(ES9.2)')//' km/s')

    end function max_defect_text
!*****************************************************************************************

!*****************************************************************************************
!>
!  A row of the trajectory's summary table: a dim name, and its value.

    function info_row(name, value) result(text)

    character(len=*),intent(in) :: name, value
    character(len=:),allocatable :: text

    text = '<tr><td style="color:'//ui_dim//'; padding-right:12px;">'//name//'</td><td>'//value//'</td></tr>'

    end function info_row
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
!  The trajectory's colors, as "Color by" says: by the distance from the Moon
!  `r` (along the colormap, nearest to farthest); by eclipses at times `t`
!  (yellow in sunlight; in the Earth's shadow cyan to violet, by how deep; gray
!  without eclipse data); or all cyan.

    subroutine trajectory_colors(t, r, colors)

    real(wp),dimension(:),intent(in) :: t  !! (n) ET (sec)
    real(wp),dimension(:),intent(in) :: r  !! (n) km
    real(wp),dimension(:,:),intent(out) :: colors  !! (3,n)

    real(wp) :: rmin, rmax, phi, w
    integer :: j

    select case (app%color_mode%currentIndex())
    case (COLOR_ALTITUDE)
        rmin = minval(r)
        rmax = maxval(r)
        do j = 1, size(r)
            colors(:,j) = colormap((r(j) - rmin) / max(rmax - rmin, 1.0e-9_wp))
        end do
    case (COLOR_ECLIPSE)
        do j = 1, size(t)
            if (.not. app%has_eclipses) then
                colors(:,j) = no_data_color
                cycle
            end if
            phi = eclipse_phi(t(j))
            if (phi >= 0.0_wp) then
                colors(:,j) = lit_color
            else
                w = min(1.0_wp, phi / min(app%ecl_phi_min, -tiny(1.0_wp)))
                colors(:,j) = (1.0_wp - w) * shadow_color_edge + w * shadow_color_deep
            end if
        end do
    case default
        do j = 1, size(r)
            colors(:,j) = line_color
        end do
    end select

    end subroutine trajectory_colors
!*****************************************************************************************

!*****************************************************************************************
!>
!  The legend for "Color by" (bottom right): the nearest and farthest distances,
!  or the eclipses' colors; none for a single color.

    subroutine show_legend()

    logical :: shown

    shown = app%loaded
    select case (app%color_mode%currentIndex())
    case (COLOR_ALTITUDE)
        call app%hud_near%set_color(cmap(:,1))
        call app%hud_far%set_color(cmap(:,size(cmap,2)))
        call app%hud_near%set_text('|r| = '//real_text(app%rmin, '(F0.0)')//' km')
        call app%hud_far%set_text('|r| = '//real_text(app%rmax, '(F0.0)')//' km')
    case (COLOR_ECLIPSE)
        if (app%has_eclipses) then
            call app%hud_near%set_color(lit_color)
            call app%hud_near%set_text('in sunlight')
            call app%hud_far%set_color(shadow_color_edge)
            if (app%ecl_phi_min < 0.0_wp) then
                call app%hud_far%set_text('in the Earth''s shadow')
            else
                call app%hud_far%set_text('(no eclipses)')
            end if
        else
            call app%hud_near%set_text('')
            call app%hud_far%set_color(no_data_color)
            call app%hud_far%set_text('no eclipse data (eclipse_*.json)')
        end if
    case default
        shown = .false.
    end select
    call app%hud_near%show(shown)
    call app%hud_far%show(shown)

    end subroutine show_legend
!*****************************************************************************************

!*****************************************************************************************
!>
!  "Color by" changed: the trajectory's points added again (from the vkframes
!  trajectories, in each frame) with the new colors.

    subroutine on_color_mode(index, data)

    integer,intent(in) :: index
    class(*),intent(inout),optional :: data

    call show_legend()
    call recolor()

    end subroutine on_color_mode
!*****************************************************************************************

!*****************************************************************************************
!>
!  The segments' start and end markers shown or not (in each frame).

    subroutine on_show_segment_markers(checked, data)

    logical,intent(in) :: checked
    class(*),intent(inout),optional :: data

    integer :: i, markers

    markers = 0
    if (checked) markers = ior(ior(MARKER_START, MARKER_INTERMEDIATE), MARKER_END)
    do i = 1, size(app%frames)
        call app%frames(i)%start_markers%set_markers(markers)
        call app%frames(i)%end_markers%set_markers(markers)
    end do

    end subroutine on_show_segment_markers
!*****************************************************************************************

!*****************************************************************************************
    subroutine recolor()  !! the trajectory colored again (in each frame)

    real(wp),dimension(:),allocatable :: t, r
    real(wp),dimension(:,:),allocatable :: pos, colors
    integer :: i, k, n

    if (.not. app%loaded) return
    n = app%frames(ROTATING)%traj%num_times()
    if (n == 0) return
    allocate(t(n), r(n), pos(3,n), colors(3,n))
    do k = 1, n
        t(k) = app%frames(ROTATING)%traj%time(k)
        r(k) = norm2(app%frames(ROTATING)%traj%position(k))  ! (from the Moon)
    end do
    call trajectory_colors(t, r, colors)
    do i = 1, size(app%frames)
        associate (traj => app%frames(i)%traj)
            if (traj%num_times() /= n) cycle
            do k = 1, n
                pos(:,k) = traj%position(k)
            end do
            call traj%clear()
            call traj%add(t, pos=pos, color=colors)
        end associate
    end do

    end subroutine recolor
!*****************************************************************************************

!*****************************************************************************************
!>
!  Show the trajectory in frame `iframe` (ROTATING, INERTIAL or SUN_EARTH), with the
!  current state following it, and the frame's central body.

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
    ! the central body: the Moon, or the Earth with its shadow
    call app%moon%show_contents(iframe /= SUN_EARTH)
    call app%moon%show_name_label(iframe /= SUN_EARTH)
    call app%earth%show_contents(iframe == SUN_EARTH)
    call app%earth%show_name_label(iframe == SUN_EARTH)
    call app%umbra%show_contents(iframe == SUN_EARTH)
    call app%penumbra%show_contents(iframe == SUN_EARTH)
    if (app%loaded) then
        call app%hud_frame%set_text(trim(frame_names(iframe))//', '//trim(frame_centers(iframe)))
    else
        call app%hud_frame%set_text('')
    end if

    end subroutine show_frame
!*****************************************************************************************

!*****************************************************************************************
!>
!  The camera put back to see the whole trajectory in the frame shown: looking
!  at the frame's center (the Moon or the Earth: the view rotates about it) from
!  an isometric direction (+x, +y, +z), far enough for the trajectory's farthest
!  point, z up.

    subroutine reset_view()

    type(View),pointer :: v
    real(wp),dimension(3) :: eye
    real(wp),dimension(3),parameter :: center = 0.0_wp  ! the Moon or the Earth
    real(wp),dimension(3),parameter :: direction = [1.0_wp, 1.0_wp, 1.0_wp] / sqrt(3.0_wp)
    real(wp),dimension(3),parameter :: up = [0.0_wp, 0.0_wp, 1.0_wp]

    if (.not. app%loaded) return
    if (app%current_frame == SUN_EARTH) then  ! (the trajectory goes around the Earth)
        eye = 3.0_wp * app%frames(app%current_frame)%radius * direction
    else  ! (the trajectory is on one side of the Moon)
        eye = 2.0_wp * app%frames(app%current_frame)%radius * direction
    end if
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
    call app%proxy%pause_time(app%pause_time%isChecked())
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
    if (app%loaded .and. app%bottom_tabs%currentIndex() == 2) call app%eclipse_plot%update()  ! (the time marker)
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
    call clear_eclipses()  ! (the iterations have none; the solution's are read when it's done)
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

    subroutine send_iteration(iter, et, r_rotating, r_inertial, r_se_rotating, &
                              seg_et, seg_r_rotating, seg_r_inertial, seg_r_se_rotating, defects)

    integer,intent(in) :: iter
    real(hwp),dimension(:),intent(in) :: et
    real(hwp),dimension(:,:),intent(in) :: r_rotating, r_inertial, r_se_rotating
    real(hwp),dimension(:,:),intent(in) :: seg_et
    real(hwp),dimension(:,:,:),intent(in) :: seg_r_rotating, seg_r_inertial, seg_r_se_rotating
    real(hwp),dimension(:,:),intent(in) :: defects

    type(iteration_points),pointer :: p
    integer(c_intptr_t) :: address
    integer(int32),dimension(2) :: halves
    integer :: i, n, nseg, npairs

    flush(output_unit)  ! (the iteration's output shown now)
    n = size(et)
    if (n < 2) return
    allocate(p)
    allocate(p%t(n), p%rot(3,n), p%inr(3,n), p%se(3,n))
    p%iter = iter
    do i = 1, n
        p%t(i)     = real(et(i), wp)
        p%rot(:,i) = real(r_rotating(:,i), wp)
        p%inr(:,i) = real(r_inertial(:,i), wp)
        p%se(:,i)  = real(r_se_rotating(:,i), wp)
    end do
    nseg = size(seg_et,2)
    allocate(p%seg_t(2,nseg), p%seg_rot(3,2,nseg), p%seg_inr(3,2,nseg), p%seg_se(3,2,nseg))
    do i = 1, nseg
        p%seg_t(:,i)     = real(seg_et(:,i), wp)
        p%seg_rot(:,:,i) = real(seg_r_rotating(:,:,i), wp)
        p%seg_inr(:,:,i) = real(seg_r_inertial(:,:,i), wp)
        p%seg_se(:,:,i)  = real(seg_r_se_rotating(:,:,i), wp)
    end do
    npairs = size(defects,2)
    allocate(p%defects(2,npairs))
    do i = 1, npairs
        p%defects(:,i) = real(defects(:,i), wp)
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

    call plot_points(p%t, p%rot, p%inr, p%se, p%seg_t, p%seg_rot, p%seg_inr, p%seg_se, p%defects, &
                     'Iteration '//int_text(p%iter), &
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
    call show_pause_state(checked)
    end subroutine on_pause_time

    subroutine show_pause_state(paused)  !! the play/pause button's symbol: what a click does
    logical,intent(in) :: paused
    if (paused) then
        call app%pause_time%setText(play_symbol)
        call app%pause_time%setToolTip('Play')
    else
        call app%pause_time%setText(pause_symbol)
        call app%pause_time%setToolTip('Pause')
    end if
    end subroutine show_pause_state

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
    type(QGridLayout) :: grid
    type(QHBoxLayout) :: hbox
    type(QPushButton) :: reset_view_button
    type(QToolButton) :: reset
    type(QLabel) :: hint
    integer :: r  !! the grid's next row
    type(QMenuBar) :: bar
    type(QMenu) :: menu
    type(QAction) :: action

    call apply_theme()
    app%main = QMainWindow()
    call app%main%setStyleSheet(ui_style_sheet)
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

    ! the controls: one grid (names in column 0, the fields in 1), in sections
    panel = QWidget()
    grid = QGridLayout(panel)
    call grid%setContentsMargins(14, 4, 14, 14)
    call grid%setSpacing(8)
    call grid%setColumnStretch(1, 1)
    r = 0

    ! the solver
    call add_section('SOLVER')
    row = QWidget()
    hbox = QHBoxLayout(row)
    call hbox%setContentsMargins(0, 0, 0, 0)
    call hbox%setSpacing(6)
    app%config_edit = QLineEdit(app%config_file)
    call app%config_edit%setPlaceholderText('config file (.json)')
    call app%config_edit%setToolTip('The solver''s config file')
    call hbox%addWidget(app%config_edit, 1)
    app%browse_button = QPushButton('Browse')
    call app%browse_button%setToolTip('Choose a config file')
    call app%browse_button%clicked%connect(on_browse)
    call hbox%addWidget(app%browse_button)
    call add_wide(row)
    app%plot_iterations = QCheckBox('Plot each iteration')
    call app%plot_iterations%setChecked(app%iterations_at_start)
    call app%plot_iterations%setToolTip('Show the trajectory of each solver iteration as it solves '// &
                                        '(propagates all the segments at each iteration).')
    call add_wide(app%plot_iterations)
    row = QWidget()
    hbox = QHBoxLayout(row)
    call hbox%setContentsMargins(0, 0, 0, 0)
    call hbox%setSpacing(6)
    app%run_button = QPushButton('Run solver')
    call app%run_button%setObjectName('primary')
    call app%run_button%clicked%connect(on_run)
    call hbox%addWidget(app%run_button, 1)
    app%stop_button = QPushButton('Stop')
    call app%stop_button%setEnabled(.false.)
    call app%stop_button%clicked%connect(on_stop)
    call hbox%addWidget(app%stop_button)
    call add_wide(row)
    app%status_label = QLabel('Choose a config file and run the solver.')
    call app%status_label%setObjectName('dim')
    call app%status_label%setWordWrap(.true.)
    call add_wide(app%status_label)

    ! the trajectory
    call add_section('TRAJECTORY')
    app%info_label = QLabel('No trajectory loaded')
    call app%info_label%setWordWrap(.true.)
    call add_wide(app%info_label)
    app%frame_choice = QComboBox()
    call app%frame_choice%addItems(frame_names)
    call app%frame_choice%setCurrentIndex(app%frame_at_start - 1)
    call show_frame(app%frame_at_start)
    call app%frame_choice%currentIndexChanged%connect(on_frame_choice)
    call add_field('Frame', app%frame_choice)
    app%color_mode = QComboBox()
    call app%color_mode%addItems(color_names)
    call app%color_mode%setCurrentIndex(app%color_at_start)
    call app%color_mode%setToolTip('Color the trajectory by the distance from the Moon (red: nearest, '// &
                                   'blue: farthest), by eclipses (yellow: in sunlight; cyan to violet: in '// &
                                   'the Earth''s shadow), or cyan.')
    call app%color_mode%currentIndexChanged%connect(on_color_mode)
    call add_field('Color', app%color_mode)
    app%show_segment_markers = QCheckBox('Segment start/end markers')
    call app%show_segment_markers%setChecked(.true.)
    call app%show_segment_markers%toggled%connect(on_show_segment_markers)
    call add_wide(app%show_segment_markers)
    reset_view_button = QPushButton('Reset view')
    call reset_view_button%clicked%connect(on_reset_view)
    call add_wide(reset_view_button)

    ! time: play/pause, back to the start, the slider; the readout; the speed
    call add_section('TIME')
    row = QWidget()
    hbox = QHBoxLayout(row)
    call hbox%setContentsMargins(0, 0, 0, 0)
    call hbox%setSpacing(6)
    app%pause_time = QToolButton()
    call app%pause_time%setCheckable(.true.)
    call app%pause_time%setChecked(.true.)  ! (paused at first)
    call show_pause_state(.true.)
    call app%pause_time%toggled%connect(on_pause_time)
    call hbox%addWidget(app%pause_time)
    reset = QToolButton()
    call reset%setText(reset_symbol)
    call reset%setToolTip('Back to the start')
    call reset%clicked%connect(on_reset_time)
    call hbox%addWidget(reset)
    app%time_slider = QSlider(Qt_Horizontal)
    call app%time_slider%setRange(0, slider_steps)
    call app%time_slider%setSingleStep(slider_steps / 100)
    call app%time_slider%setPageStep(slider_steps / 10)
    call app%time_slider%valueChanged%connect(on_time_slider)
    call app%time_slider%sliderPressed%connect(on_slider_pressed)
    call app%time_slider%sliderReleased%connect(on_slider_released)
    call hbox%addWidget(app%time_slider, 1)
    call add_wide(row)
    app%time_label = QLabel('')
    call app%time_label%setObjectName('dim')
    call add_wide(app%time_label)
    app%time_scale = QDoubleSpinBox()
    call app%time_scale%setRange(0.0_c_double, 1.0e9_c_double)
    call app%time_scale%setValue(1.0_c_double)  ! (real time)
    call app%time_scale%setSuffix(' x')
    call app%time_scale%setToolTip('Simulation seconds per second (type a value, or use the wheel)')
    call app%time_scale%valueChanged%connect(on_time_scale)
    call add_field('Speed', app%time_scale)

    ! the mouse
    call add_section('VIEW')
    hint = QLabel('Left drag: rotate<br>Right drag or wheel: zoom<br>Shift + left drag, or middle drag: pan')
    call hint%setObjectName('dim')
    call add_wide(hint)
    call grid%setRowStretch(r, 1)

    scroll = QScrollArea()
    call scroll%setWidgetResizable(.true.)
    call scroll%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
    call scroll%setWidget(panel)
    call scroll%setMinimumWidth(340)  ! (the panel, and the scroll bar beside it)
    dock = QDockWidget('Controls', app%main)
    call dock%setFeatures(ior(QDockWidget_DockWidgetMovable, QDockWidget_DockWidgetFloatable))
    call dock%setWidget(scroll)
    call app%main%addDockWidget(Qt_LeftDockWidgetArea, dock)

    ! the solver's output
    app%log = QPlainTextEdit()
    call app%log%setReadOnly(.true.)
    call app%log%setMaximumBlockCount(20000)
    call app%log%setStyleSheet('font-family: Menlo, Monaco, Consolas, monospace; font-size: 11pt;')
    ! the plots
    app%defects_plot = QtfWidget()
    call app%defects_plot%setMinimumSize(300, 150)
    call app%defects_plot%paintEvent%connect(on_paint_defects)
    app%eclipse_plot = QtfWidget()
    call app%eclipse_plot%setMinimumSize(300, 150)
    call app%eclipse_plot%paintEvent%connect(on_paint_eclipses)
    app%bottom_tabs = QTabWidget()
    call app%bottom_tabs%addTab(app%log, 'Solver output')
    call app%bottom_tabs%addTab(app%defects_plot, 'Defects')
    call app%bottom_tabs%addTab(app%eclipse_plot, 'Eclipses')
    dock = QDockWidget('Output and plots', app%main)
    call dock%setFeatures(ior(QDockWidget_DockWidgetMovable, QDockWidget_DockWidgetFloatable))
    call dock%setWidget(app%bottom_tabs)
    call app%main%addDockWidget(Qt_BottomDockWidgetArea, dock)
    call app%bottom_tabs%setCurrentIndex(app%tab_at_start)

    ! the iterations' trajectories (from the solver's thread), and its output (from the reader's)
    app%job%iteration = QSignalInt2(app%main)
    call app%job%iteration%connect(on_iteration)
    app%rjob%text = QSignalString(app%main)
    call app%rjob%text%connect(on_output)

    contains

        subroutine add_section(title)  !! a section's heading (across the grid)
        character(len=*),intent(in) :: title
        type(QLabel) :: label
        label = QLabel(title)
        call label%setObjectName('section')
        call grid%addWidget(label, r, 0, 1, 2)
        r = r + 1
        end subroutine add_section

        subroutine add_field(name, widget)  !! a field and its name
        character(len=*),intent(in) :: name
        class(QWidget),intent(in) :: widget
        type(QLabel) :: label
        label = QLabel(name)
        call grid%addWidget(label, r, 0)
        call grid%addWidget(widget, r, 1)
        r = r + 1
        end subroutine add_field

        subroutine add_wide(widget)  !! a widget across the grid
        class(QWidget),intent(in) :: widget
        call grid%addWidget(widget, r, 0, 1, 2)
        r = r + 1
        end subroutine add_wide

    end subroutine build_ui
!*****************************************************************************************

!*****************************************************************************************
!>
!  The application's look: the Fusion style (the same on all platforms), with a
!  dark palette (the controls' details are in `ui_style_sheet`).

    subroutine apply_theme()

    type(QPalette) :: palette

    if (.not. QApplication_setStyle('Fusion')) return  ! (keep the platform's look)
    palette = QPalette()
    call palette%setColor(QPalette_Window, plot_color(ui_window))
    call palette%setColor(QPalette_WindowText, plot_color(ui_text))
    call palette%setColor(QPalette_Base, plot_color(ui_base))
    call palette%setColor(QPalette_AlternateBase, plot_color(ui_button))
    call palette%setColor(QPalette_Text, plot_color(ui_text))
    call palette%setColor(QPalette_Button, plot_color(ui_button))
    call palette%setColor(QPalette_ButtonText, plot_color(ui_text))
    call palette%setColor(QPalette_Highlight, plot_color(ui_accent))
    call palette%setColor(QPalette_HighlightedText, QColor(255, 255, 255))
    call palette%setColor(QPalette_ToolTipBase, plot_color(ui_base))
    call palette%setColor(QPalette_ToolTipText, plot_color(ui_text))
    call palette%setColor(QPalette_PlaceholderText, plot_color(ui_disabled))
    call palette%setColor(QPalette_Text, plot_color(ui_disabled), QPalette_Disabled)
    call palette%setColor(QPalette_WindowText, plot_color(ui_disabled), QPalette_Disabled)
    call palette%setColor(QPalette_ButtonText, plot_color(ui_disabled), QPalette_Disabled)
    call QApplication_setPalette(palette)
    call palette%delete()

    end subroutine apply_theme
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
!>
!  The defects plot: |dr| (km) and |dv| (km/s) at each interface, log scale.

    subroutine on_paint_defects(painter, data)

    type(QPainter),intent(inout) :: painter
    class(*),intent(inout),optional :: data

    real(wp),parameter :: floor = 1.0e-16_wp  !! (zero on the log scale)
    real(wp),dimension(:),allocatable :: x, y
    real(wp) :: ymin, ymax
    integer :: i, n

    call begin_plot(painter)
    if (.not. allocated(app%defects)) then
        call plot_message(painter, 'No defects: load or solve a trajectory.')
        return
    end if
    n = size(app%defects,2)
    if (n == 0) then
        call plot_message(painter, 'No defects: load or solve a trajectory.')
        return
    end if
    allocate(x(n), y(n))
    do i = 1, n
        x(i) = real(i, wp)
    end do
    ymin = log10(max(floor, minval(app%defects)))
    ymax = log10(max(floor, maxval(app%defects)))
    ymin = floor_to(ymin, 1.0_wp)
    ymax = max(ymin + 1.0_wp, ceiling_to(ymax, 1.0_wp))
    call plot_axes(painter, 1.0_wp, max(2.0_wp, real(n, wp)), ymin, ymax, .true., &
                   'Interface (segments 2k-1, 2k)', 'Defect', &
                   'Constraint defects: '//app%defects_title)
    do i = 1, n
        y(i) = log10(max(floor, app%defects(1,i)))
    end do
    call plot_series(painter, x, y, plot_color(plot_dr), 1, '|dr| (km)')
    do i = 1, n
        y(i) = log10(max(floor, app%defects(2,i)))
    end do
    call plot_series(painter, x, y, plot_color(plot_dv), 2, '|dv| (km/s)')

    end subroutine on_paint_defects
!*****************************************************************************************

!*****************************************************************************************
!>
!  The eclipse plot: min(phi, 0) over time (< 0: in the Earth's shadow), at most
!  one point per pixel column (the deepest), with the current time.

    subroutine on_paint_eclipses(painter, data)

    type(QPainter),intent(inout) :: painter
    class(*),intent(inout),optional :: data

    integer,parameter :: nbins = 2000
    real(wp),dimension(:),allocatable :: x, y
    real(wp) :: xmax, ymin, t, xpix
    integer :: i, k, n

    call begin_plot(painter)
    if (.not. app%has_eclipses) then
        call plot_message(painter, 'No eclipse data: the solution''s eclipse_*.json file '// &
                                   '(generate_eclipse_files, eclipse_filetype 2).')
        return
    end if
    n = size(app%ecl_t)
    xmax = max((app%ecl_t(n) - app%ecl_t(1)) / day, 1.0e-6_wp)
    ymin = min(app%ecl_phi_min, -1.0e-6_wp)
    ymin = ymin * 1.1_wp
    if (app%ecl_phi_min < 0.0_wp) then
        call plot_axes(painter, 0.0_wp, xmax, ymin, 0.0_wp - 0.05_wp * ymin, .false., &
                       'Time (days)', 'phi', 'Eclipses: phi < 0 in the Earth''s shadow')
    else
        call plot_axes(painter, 0.0_wp, xmax, -1.0_wp, 1.0_wp, .false., &
                       'Time (days)', 'phi', 'Eclipses: none (phi > 0 throughout)')
    end if

    ! the deepest value in each of nbins time bins
    allocate(x(nbins), y(nbins))
    x = 0.0_wp
    y = huge(1.0_wp)
    do i = 1, n
        k = 1 + min(nbins - 1, int(real(nbins, wp) * (app%ecl_t(i) - app%ecl_t(1)) / day / xmax))
        x(k) = (app%ecl_t(i) - app%ecl_t(1)) / day
        y(k) = min(y(k), min(app%ecl_phi(i), 0.0_wp))
    end do
    k = 0
    do i = 1, nbins
        if (y(i) == huge(1.0_wp)) cycle
        k = k + 1
        x(k) = x(i)
        y(k) = y(i)
    end do
    call plot_series(painter, x(1:k), y(1:k), plot_color(plot_phi), 0, '')

    ! the current time
    if (app%loaded .and. app%proxy%is_open()) then
        t = (app%proxy%time() - app%ecl_t(1)) / day
        if (t >= 0.0_wp .and. t <= xmax) then
            xpix = app%plot_x0 + (app%plot_x1 - app%plot_x0) * t / xmax
            call painter%setPen(QPen(plot_color(plot_now), 1.5_wp))
            call painter%drawLine(xpix, app%plot_y0, xpix, app%plot_y1)
        end if
    end if

    end subroutine on_paint_eclipses
!*****************************************************************************************

!*****************************************************************************************
!  plotting (QPainter): the plot area and its data ranges are kept in `app` while painting

    subroutine begin_plot(painter)  !! the background, and the plot area
    type(QPainter),intent(inout) :: painter
    integer :: w, h
    w = painter%deviceWidth()
    h = painter%deviceHeight()
    call painter%setRenderHint(QPainter_Antialiasing)
    call painter%fillRect(0, 0, w, h, plot_color(plot_background))
    app%plot_x0 = 70.0_wp
    app%plot_x1 = max(app%plot_x0 + 10.0_wp, real(w, wp) - 20.0_wp)
    app%plot_y0 = 30.0_wp
    app%plot_y1 = max(app%plot_y0 + 10.0_wp, real(h, wp) - 40.0_wp)
    end subroutine begin_plot

    subroutine plot_message(painter, text)  !! a message in place of a plot
    type(QPainter),intent(inout) :: painter
    character(len=*),intent(in) :: text
    call painter%setPen(plot_color(plot_dim_text))
    call painter%setFont(QFont(pointSize=11.0_wp))
    call painter%drawText(0.0_wp, 0.0_wp, real(painter%deviceWidth(), wp), real(painter%deviceHeight(), wp), &
                          Qt_AlignCenter, text)
    end subroutine plot_message

    !> The axes (box, grid, tick labels, axis labels, title) for data ranges
    !> [xmin, xmax] by [ymin, ymax]; with `logy`, the y values are log10 (ticks
    !> at each decade, labeled 1e-n).
    subroutine plot_axes(painter, xmin, xmax, ymin, ymax, logy, xlabel, ylabel, title)
    type(QPainter),intent(inout) :: painter
    real(wp),intent(in) :: xmin, xmax, ymin, ymax
    logical,intent(in) :: logy
    character(len=*),intent(in) :: xlabel, ylabel, title
    real(wp) :: step, v, px, py
    character(len=32) :: label
    integer :: i, k

    app%plot_xmin = xmin
    app%plot_xmax = xmax
    app%plot_ymin = ymin
    app%plot_ymax = ymax
    call painter%setFont(QFont(pointSize=9.0_wp))

    ! y grid and labels
    if (logy) then
        step = max(1.0_wp, real(ceiling((ymax - ymin) / 8.0_wp), wp))
    else
        step = nice_step((ymax - ymin) / 5.0_wp)
    end if
    k = 0
    v = ceiling_to(ymin, step)
    do while (v <= ymax + 1.0e-9_wp * step .and. k < 50)
        py = plot_y(v)
        call painter%setPen(QPen(plot_color(plot_grid), 1.0_wp, Qt_DotLine))
        call painter%drawLine(app%plot_x0, py, app%plot_x1, py)
        if (logy) then
            write(label, '(A,I0)') '1e', nint(v)
        else
            label = number_label(v, step)
        end if
        call painter%setPen(plot_color(plot_text))
        call painter%drawText(0.0_wp, py - 8.0_wp, app%plot_x0 - 6.0_wp, 16.0_wp, &
                              ior(Qt_AlignRight, Qt_AlignVCenter), trim(label))
        v = v + step
        k = k + 1
    end do

    ! x grid and labels
    step = nice_step((xmax - xmin) / 8.0_wp)
    k = 0
    v = ceiling_to(xmin, step)
    do while (v <= xmax + 1.0e-9_wp * step .and. k < 50)
        px = plot_x(v)
        call painter%setPen(QPen(plot_color(plot_grid), 1.0_wp, Qt_DotLine))
        call painter%drawLine(px, app%plot_y0, px, app%plot_y1)
        label = number_label(v, step)
        call painter%setPen(plot_color(plot_text))
        call painter%drawText(px - 40.0_wp, app%plot_y1 + 2.0_wp, 80.0_wp, 16.0_wp, Qt_AlignCenter, trim(label))
        v = v + step
        k = k + 1
    end do

    ! box, labels, title
    call painter%setPen(QPen(plot_color(plot_frame), 1.0_wp))
    call painter%setBrush(QBrush(style=Qt_NoBrush))
    call painter%drawRect(app%plot_x0, app%plot_y0, app%plot_x1 - app%plot_x0, app%plot_y1 - app%plot_y0)
    call painter%drawText(app%plot_x0, app%plot_y1 + 18.0_wp, app%plot_x1 - app%plot_x0, 16.0_wp, &
                          Qt_AlignCenter, xlabel)
    call painter%drawText(4.0_wp, app%plot_y0 - 22.0_wp, ylabel)
    call painter%setPen(plot_color(plot_title))
    call painter%setFont(QFont(pointSize=11.0_wp, bold=.true.))
    call painter%drawText(app%plot_x0, 4.0_wp, app%plot_x1 - app%plot_x0, 20.0_wp, Qt_AlignCenter, title)
    call painter%setFont(QFont(pointSize=9.0_wp))
    app%plot_legend_row = 0
    i = 0
    end subroutine plot_axes

    !> A line through the points (data coordinates), with markers if `markers` > 0
    !> (a size in pixels) and an entry in the legend (top right) if `name` isn't ''.
    subroutine plot_series(painter, x, y, color, markers, name)
    type(QPainter),intent(inout) :: painter
    real(wp),dimension(:),intent(in) :: x, y
    type(QColor),intent(in) :: color
    integer,intent(in) :: markers
    character(len=*),intent(in) :: name
    real(wp),dimension(:),allocatable :: px, py
    real(wp) :: ly
    integer :: i, n
    n = min(size(x), size(y))
    if (n == 0) return
    allocate(px(n), py(n))
    do i = 1, n
        px(i) = plot_x(x(i))
        py(i) = plot_y(y(i))
    end do
    call painter%setPen(QPen(color, 1.5_wp))
    if (n > 1) call painter%drawPolyline(px, py)
    if (markers > 0 .and. n <= 2000) then
        call painter%setBrush(color)
        do i = 1, n
            call painter%drawEllipse(px(i) - markers, py(i) - markers, 2.0_wp * markers, 2.0_wp * markers)
        end do
        call painter%setBrush(QBrush(style=Qt_NoBrush))
    end if
    if (len(name) > 0) then
        app%plot_legend_row = app%plot_legend_row + 1
        ly = app%plot_y0 + 6.0_wp + 16.0_wp * (app%plot_legend_row - 1)
        call painter%fillRect(app%plot_x1 - 114.0_wp, ly - 1.0_wp, 112.0_wp, 16.0_wp, plot_color(plot_background))
        call painter%setPen(QPen(color, 1.5_wp))
        call painter%drawLine(app%plot_x1 - 110.0_wp, ly + 7.0_wp, app%plot_x1 - 90.0_wp, ly + 7.0_wp)
        call painter%setPen(plot_color(plot_text))
        call painter%drawText(app%plot_x1 - 86.0_wp, ly, 84.0_wp, 14.0_wp, &
                              ior(Qt_AlignLeft, Qt_AlignVCenter), name)
    end if
    end subroutine plot_series

    function plot_color(rgb) result(c)  !! a plot palette color as a QColor
    integer,dimension(3),intent(in) :: rgb
    type(QColor) :: c
    c = QColor(rgb(1), rgb(2), rgb(3))
    end function plot_color

    real(wp) function plot_x(v)  !! data x to pixels
    real(wp),intent(in) :: v
    plot_x = app%plot_x0 + (app%plot_x1 - app%plot_x0) * (v - app%plot_xmin) / &
             max(app%plot_xmax - app%plot_xmin, tiny(1.0_wp))
    end function plot_x

    real(wp) function plot_y(v)  !! data y to pixels (up)
    real(wp),intent(in) :: v
    plot_y = app%plot_y1 - (app%plot_y1 - app%plot_y0) * (v - app%plot_ymin) / &
             max(app%plot_ymax - app%plot_ymin, tiny(1.0_wp))
    end function plot_y

    real(wp) function nice_step(raw)  !! 1, 2 or 5 times a power of 10, >= raw
    real(wp),intent(in) :: raw
    real(wp) :: p, f
    if (raw <= 0.0_wp) then
        nice_step = 1.0_wp
        return
    end if
    p = 10.0_wp**floor(log10(raw))
    f = raw / p
    if (f <= 1.0_wp) then
        nice_step = p
    else if (f <= 2.0_wp) then
        nice_step = 2.0_wp * p
    else if (f <= 5.0_wp) then
        nice_step = 5.0_wp * p
    else
        nice_step = 10.0_wp * p
    end if
    end function nice_step

    real(wp) function floor_to(v, step)  !! v rounded down to a multiple of step
    real(wp),intent(in) :: v, step
    floor_to = step * floor(v / step)
    end function floor_to

    real(wp) function ceiling_to(v, step)  !! v rounded up to a multiple of step
    real(wp),intent(in) :: v, step
    ceiling_to = step * ceiling(v / step - 1.0e-9_wp)
    end function ceiling_to

    function number_label(v, step) result(label)  !! a tick label (decimals for the step)
    real(wp),intent(in) :: v, step
    character(len=32) :: label
    character(len=16) :: fmt
    integer :: decimals
    if (abs(v) >= 1.0e6_wp .or. (abs(v) < 1.0e-3_wp .and. v /= 0.0_wp .and. step < 1.0e-3_wp)) then
        write(label, '(ES10.2)') v
    else
        decimals = max(0, -floor(log10(step) + 1.0e-9_wp))
        write(fmt, '(A,I0,A)') '(F0.', decimals, ')'
        write(label, fmt) v
        if (label(1:1) == '.') label = '0'//label
        if (label(1:2) == '-.') label = '-0'//label(2:)
        if (len_trim(label) > 1 .and. label(len_trim(label):len_trim(label)) == '.') &
            label = label(1:len_trim(label)-1)  ! ("50." -> "50")
        if (abs(v) < 0.5_wp * step) label = '0'
    end if
    label = adjustl(label)
    end function number_label
!*****************************************************************************************

!*****************************************************************************************
!  small utilities

    subroutine set_status(text)
    character(len=*),intent(in) :: text
    call app%status_label%setText(breakable(text))
    end subroutine set_status

    function breakable(text) result(b)  !! the text, with line breaks allowed after _ and / (for file names)
    character(len=*),intent(in) :: text
    character(len=:),allocatable :: b
    character(len=*),parameter :: zero_width_space = '​'
    integer :: i
    b = ''
    do i = 1, len(text)
        b = b//text(i:i)
        if (text(i:i) == '_' .or. text(i:i) == '/') b = b//zero_width_space
    end do
    end function breakable

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
!  Arguments: `[--dir DIR] [config.json] [--run] [--iterations] [--traj traj_file.json]
!  [--frame 1|2|3] [--color altitude|eclipse|none] [--tab 0|1|2]`.
!  `--dir` sets the current directory (the halo root, which the config files'
!  paths are relative to; give it first: the other paths are relative to it).
!  The config file defaults to `examples/example_sparse.json` (if it exists);
!  `--run` runs the solver on it at startup; `--iterations` checks "Plot each
!  iteration"; `--traj` plots a trajectory file at startup; `--frame` the frame shown
!  (1: Earth-Moon rotating, 2: inertial, 3: Sun-Earth rotating); `--color` the
!  coloring; `--tab` the bottom tab (0: solver output, 1: defects, 2: eclipses).
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
        case ('--frame')  ! 1: Earth-Moon rotating, 2: inertial, 3: Sun-Earth rotating
            call get_command_argument(i + 1, arg)
            read(arg, *) app%frame_at_start
            app%frame_at_start = max(1, min(3, app%frame_at_start))
            i = i + 1
        case ('--color')
            call get_command_argument(i + 1, arg)
            select case (arg)
            case ('eclipse', 'eclipses')
                app%color_at_start = COLOR_ECLIPSE
            case ('none', 'cyan', 'red')
                app%color_at_start = COLOR_NONE
            case default
                app%color_at_start = COLOR_ALTITUDE
            end select
            i = i + 1
        case ('--tab')  ! 0: solver output, 1: defects, 2: eclipses
            call get_command_argument(i + 1, arg)
            read(arg, *) app%tab_at_start
            i = i + 1
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
