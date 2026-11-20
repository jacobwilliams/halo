!*****************************************************************************************
!>
!  A GUI for the Halo solver: choose a config file (e.g. `examples/example_sparse.json`),
!  run the solver, and plot the resulting trajectory in 3D.
!
!  The window is made with [qt-fortran](https://github.com/jacobwilliams/qt-fortran),
!  and the 3D view with [vkframes](https://github.com/jacobwilliams/vkframes).
!
!  * The solver is run as a separate process (`halo_solver config.json --gui`, the
!    executable next to this one), in the current directory (the config files' paths
!    are relative to it). Its output is shown in the "Solver output" panel, and it can
!    be stopped. With `--gui`, it always writes the JSON trajectory file (`traj_*.json`)
!    and names it in its output; this is plotted when the solver is done.
!  * File > Open Trajectory plots an existing `traj_*.json` file without solving.
!  * The trajectory is drawn in the Moon-centered Earth-Moon rotating frame or the
!    Moon-centered inertial (J2000) frame (as in `python/plot_utilities.py`), colored by
!    the distance from the Moon, with a body following it in time.
!
!### Usage
!
!```
!    > pixi run gui                                    (build and run)
!    > python tools/fpm_gui.py run halo_gui -- [config.json] [--run] [--traj traj_file.json]
!```
!
!### Author
!  * Jacob Williams : Nov. 2026

    module halo_gui_module

    use qt
    use vkframes
    use vkframes_qt,  only: connect_view, ViewLink, DemoOptions
    use vkf_text,     only: get_font
    use json_module,  only: json_core, json_value
    use iso_c_binding

    implicit none

    private

    public :: app, build_scene, build_ui, start_up

    integer,parameter :: slider_steps = 1000  !! the time slider's range
    real(wp),parameter :: moon_radius = 1737.4_wp  !! km
    real(wp),parameter :: day = 86400.0_wp  !! sec
    character(len=*),parameter :: traj_marker = ' * Trajectory file: '  !! the solver's line naming the file

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

    !> The GUI, the scene and the solver process.
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
        type(QPushButton) :: run_button, stop_button
        type(QLabel) :: status_label, info_label
        type(QComboBox) :: frame_choice
        type(QDoubleSpinBox) :: time_scale
        type(QCheckBox) :: pause_time
        type(QLabel) :: time_label
        type(QSlider) :: time_slider
        logical :: slider_held = .false.
        logical :: slider_setting = .false.
        type(QPlainTextEdit) :: log
        ! solver process
        type(QProcess) :: process
        character(len=:),allocatable :: out_buffer  !! output not yet shown (no newline yet)
        character(len=:),allocatable :: traj_file  !! the trajectory file the solver wrote (if any)
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
        real(wp) :: extent = 1.0e5_wp  !! the farthest point from the Moon (km)
        ! startup options
        character(len=:),allocatable :: config_file, initial_traj
        logical :: run_at_start = .false.  !! --run: run the solver at startup
        type(DemoOptions) :: options  !! --capture, --quit-after (for testing)
    end type app_t

    type(app_t),target,save :: app

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
    call app%view%set_default_view_distance(3.0_wp * app%extent)

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
!  `segs` of segments, each with `et` and the positions in each frame), and put
!  it in the scene. The segments are joined in time order (backward-propagated
!  ones reversed), as in `python/plot_utilities.py`.

    logical function load_trajectory(file) result(ok)

    character(len=*),intent(in) :: file  !! the file to read

    type(json_core) :: json
    type(json_value),pointer :: root, segs, seg
    real(wp),dimension(:),allocatable :: et, xr, yr, zr, xi, yi, zi
    real(wp),dimension(:),allocatable :: t_all, r_all
    real(wp),dimension(:,:),allocatable :: rot_all, inr_all, colors
    integer :: nsegs, iseg, n, ntot, k, j
    logical :: found, status_ok
    character(len=:),allocatable :: msg
    real(wp) :: rmin, rmax

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
        n = minval([size(et), size(xr), size(yr), size(zr), size(xi), size(yi), size(zi)])
        if (n == 0) cycle
        if (n > 1) then
            if (et(1) > et(n)) then  ! backward-propagated segment
                et(1:n) = et(n:1:-1)
                xr(1:n) = xr(n:1:-1);  yr(1:n) = yr(n:1:-1);  zr(1:n) = zr(n:1:-1)
                xi(1:n) = xi(n:1:-1);  yi(1:n) = yi(n:1:-1);  zi(1:n) = zi(n:1:-1)
            end if
        end if
        call grow(ntot + n)
        do k = 1, n
            ! skip points not after the last one (the segments' shared ends)
            if (ntot > 0) then
                if (et(k) <= t_all(ntot)) cycle
            end if
            ntot = ntot + 1
            t_all(ntot) = et(k)
            rot_all(:,ntot) = [xr(k), yr(k), zr(k)]
            inr_all(:,ntot) = [xi(k), yi(k), zi(k)]
        end do
    end do
    call json%destroy(root)

    if (ntot < 2) then
        call error_message('No trajectory points in '//file)
        return
    end if

    ! colors: distance from the Moon (the same in both frames)
    r_all = norm2(rot_all(:,1:ntot), dim=1)
    rmin = minval(r_all)
    rmax = maxval(r_all)
    allocate(colors(3,ntot))
    do j = 1, ntot
        colors(:,j) = colormap((r_all(j) - rmin) / max(rmax - rmin, 1.0e-9_wp))
    end do

    ! into the scene:
    call app%frames(ROTATING)%traj%clear()
    call app%frames(ROTATING)%traj%add(t_all(1:ntot), pos=rot_all(:,1:ntot), color=colors)
    call app%frames(INERTIAL)%traj%clear()
    call app%frames(INERTIAL)%traj%add(t_all(1:ntot), pos=inr_all(:,1:ntot), color=colors)

    call set_bounds(app%frames(ROTATING), rot_all(:,1:ntot))
    call set_bounds(app%frames(INERTIAL), inr_all(:,1:ntot))
    app%t0 = t_all(1)
    app%tf = t_all(ntot)
    app%extent = max(maxval(norm2(inr_all(:,1:ntot), dim=1)), rmax, 2.0_wp * moon_radius)
    call app%body%set_radius(0.01_wp * app%extent)
    app%loaded = .true.

    call app%hud_title%set_text(basename(file))
    call app%hud_near%set_text('|r| = '//real_text(rmin, '(F0.0)')//' km')
    call app%hud_far%set_text('|r| = '//real_text(rmax, '(F0.0)')//' km')
    call app%info_label%setText(basename(file)//'<br>'// &
                                int_text(ntot)//' points, '//int_text(nsegs)//' segments<br>'// &
                                real_text((app%tf - app%t0)/day, '(F0.2)')//' days<br>'// &
                                '|r|: '//real_text(rmin, '(F0.0)')//' to '//real_text(rmax, '(F0.0)')//' km')
    call app%main%setWindowTitle('HALO - '//basename(file))

    ! time and view:
    call app%time_scale%setValue(real(max(1.0_wp, (app%tf - app%t0) / 60.0_wp), c_double))  ! whole thing in ~1 min
    call update_time_range()
    call show_frame(app%current_frame)
    call reset_view()
    ok = .true.

    contains

        subroutine set_bounds(f, pos)  !! the bounding box of the trajectory (including the Moon)
            type(frame_t),intent(inout) :: f
            real(wp),dimension(:,:),intent(in) :: pos
            real(wp),dimension(3) :: lo, hi
            lo = min(minval(pos, dim=2), -moon_radius)
            hi = max(maxval(pos, dim=2),  moon_radius)
            f%center = 0.5_wp * (lo + hi)
            f%radius = 0.5_wp * norm2(hi - lo)
        end subroutine set_bounds

        subroutine grow(nmax)  !! make room for nmax points
            integer,intent(in) :: nmax
            integer :: m
            if (size(t_all) >= nmax) return
            m = max(nmax, 2 * size(t_all))
            t_all = [t_all, [(0.0_wp, j = 1, m - size(t_all))]]
            rot_all = reshape(rot_all, [3, m], pad=[0.0_wp])
            inr_all = reshape(inr_all, [3, m], pad=[0.0_wp])
        end subroutine grow

    end function load_trajectory
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
!  trajectory given on the command line, if any.

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
        if (load_trajectory(app%initial_traj)) call set_status('Loaded '//basename(app%initial_traj))
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

    file = QFileDialog_getOpenFileName('Open HALO trajectory file', QDir_currentPath(), &
                                       'Trajectory files (traj_*.json);;JSON files (*.json);;All files (*)', &
                                       app%main)
    if (len(file) == 0) return
    if (load_trajectory(file)) call set_status('Loaded '//basename(file))

    end subroutine on_open_trajectory
!*****************************************************************************************

!*****************************************************************************************
!>
!  Run the solver on the config file: `halo_solver <config> --gui` (next to this
!  program), in the current directory, its output shown as it comes.

    subroutine on_run(data)

    class(*),intent(inout),optional :: data

    character(len=:),allocatable :: config, solver
    logical :: exists
    type(QProcessEnvironment) :: env

    if (app%process%state() /= QProcess_NotRunning) return
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
    solver = QCoreApplication_applicationDirPath()//'/halo_solver'
    inquire(file=solver, exist=exists)
    if (.not. exists) then
        call error_message('The solver was not found: '//solver)
        return
    end if

    call app%log%clear()
    app%out_buffer = ''
    if (allocated(app%traj_file)) deallocate(app%traj_file)
    call app%log%appendPlainText('$ '//solver//' '//config//' --gui')

    env = QProcessEnvironment_systemEnvironment()
    call env%insert('GFORTRAN_UNBUFFERED_PRECONNECTED', 'y')  ! (output shown as it is written)
    call app%process%setProcessEnvironment(env)
    call env%delete()
    call app%process%setWorkingDirectory(QDir_currentPath())
    call app%process%start(solver, [character(len=max(len(config),5)) :: config, '--gui'])

    call app%run_button%setEnabled(.false.)
    call app%stop_button%setEnabled(.true.)
    call app%config_edit%setEnabled(.false.)
    call set_status('Solving '//basename(config)//'...')

    end subroutine on_run
!*****************************************************************************************

!*****************************************************************************************
    subroutine on_stop(data)
    class(*),intent(inout),optional :: data
    if (app%process%state() == QProcess_NotRunning) return
    call app%log%appendPlainText('*** stopped by the user')
    call app%process%kill()
    end subroutine on_stop
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver wrote output: show the complete lines (and look for the trajectory file).

    subroutine on_solver_output(data)

    class(*),intent(inout),optional :: data

    app%out_buffer = app%out_buffer//app%process%readAllStandardOutput()
    call flush_lines(.false.)

    end subroutine on_solver_output
!*****************************************************************************************

!*****************************************************************************************
!>
!  Show the output's complete lines in the log (with `all`: what's left too).

    subroutine flush_lines(all)

    logical,intent(in) :: all

    integer :: i
    character(len=:),allocatable :: line

    do
        i = index(app%out_buffer, new_line('a'))
        if (i == 0) exit
        line = app%out_buffer(1:i-1)
        app%out_buffer = app%out_buffer(i+1:)
        call add_line(line)
    end do
    if (all .and. len(app%out_buffer) > 0) then
        call add_line(app%out_buffer)
        app%out_buffer = ''
    end if

    contains

        subroutine add_line(line)
            character(len=*),intent(in) :: line
            character(len=:),allocatable :: l
            l = line
            if (len(l) > 0) then
                if (l(len(l):len(l)) == achar(13)) l = l(1:len(l)-1)
            end if
            call app%log%appendPlainText(l)
            if (index(l, traj_marker) == 1) app%traj_file = trim(adjustl(l(len(traj_marker)+1:)))
        end subroutine add_line

    end subroutine flush_lines
!*****************************************************************************************

!*****************************************************************************************
!>
!  The solver finished: plot its trajectory (if it wrote one).

    subroutine on_solver_finished(exit_code, exit_status, data)

    integer,intent(in) :: exit_code, exit_status
    class(*),intent(inout),optional :: data

    character(len=:),allocatable :: file

    app%out_buffer = app%out_buffer//app%process%readAllStandardOutput()
    call flush_lines(.true.)
    call app%run_button%setEnabled(.true.)
    call app%stop_button%setEnabled(.false.)
    call app%config_edit%setEnabled(.true.)

    if (exit_status /= QProcess_NormalExit) then
        call set_status('The solver was stopped or crashed.')
        return
    end if
    if (exit_code /= 0) then
        call set_status('The solver failed (exit code '//int_text(exit_code)//'). See its output.')
        return
    end if
    if (.not. allocated(app%traj_file)) then
        call set_status('Solver done, but it wrote no trajectory file.')
        return
    end if
    file = app%traj_file
    if (file(1:1) /= '/') file = app%process%workingDirectory()//'/'//file
    if (load_trajectory(file)) call set_status('Done: plotted '//app%traj_file)

    end subroutine on_solver_finished
!*****************************************************************************************

!*****************************************************************************************
    subroutine on_solver_error(error, data)
    integer,intent(in) :: error
    class(*),intent(inout),optional :: data
    if (error == QProcess_FailedToStart) then
        call app%run_button%setEnabled(.true.)
        call app%stop_button%setEnabled(.false.)
        call app%config_edit%setEnabled(.true.)
        call error_message('The solver could not be started: '//app%process%errorString())
    end if
    end subroutine on_solver_error
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
!  Closing: stop a running solver (after asking), and close the view.

    subroutine on_close(event, data)

    type(QCloseEvent),intent(in) :: event
    class(*),intent(inout),optional :: data

    integer :: answer

    if (app%process%state() /= QProcess_NotRunning) then
        answer = QMessageBox_question('Solver running', 'The solver is still running. Stop it and quit?', &
                                      ior(QMessageBox_Yes, QMessageBox_No), QMessageBox_No, app%main)
        if (answer /= QMessageBox_Yes) then
            call event%ignore()
            return
        end if
        call app%process%kill()
        if (app%process%waitForFinished(3000)) continue
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
    type(QPushButton) :: browse, reset, reset_view_button
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
    browse = QPushButton('Browse...')
    call browse%clicked%connect(on_browse)
    call hbox%addWidget(browse)
    call form%addRow('Config file', row)
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

    ! the solver process
    app%process = QProcess(app%main)
    call app%process%setProcessChannelMode(QProcess_MergedChannels)
    call app%process%readyReadStandardOutput%connect(on_solver_output)
    call app%process%finished%connect(on_solver_finished)
    call app%process%errorOccurred%connect(on_solver_error)

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
!  Arguments: `[config.json] [--run] [--traj traj_file.json]`. The config file
!  defaults to `examples/example_sparse.json` (if it exists); `--run` runs the
!  solver on it at startup; `--traj` plots a trajectory file at startup.
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
        case ('--run')
            app%run_at_start = .true.
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
