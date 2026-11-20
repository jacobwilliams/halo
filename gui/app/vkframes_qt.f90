!> Copied from vkframes (https://github.com/jacobwilliams/vkframes, rev a2f551a):
!> demos/qt_common/src/vkframes_qt.f90 [MIT license]. It is not part of the
!> vkframes library itself, and fpm cannot depend on a subfolder of a git repo.
!> Changed here: Shift + left drag pans (sent as the middle button, VSG's
!> trackball's pan button), for mice and trackpads without a middle button.
!>
!> What the vkframes demos in qt-fortran GUIs share (vkframes itself doesn't
!> depend on Qt): Qt's mouse, wheel and key events turned into VSG's for a
!> WindowProxy rendering into a QWindow (`forward_mouse`, `forward_wheel`,
!> `forward_key`; the conversions `vsg_button`, `vsg_mask`, `vsg_key`), and
!> the window's native surface going away (`watch_surface`); and all of
!> those connected at once for a QWindow showing a WindowProxy
!> (`connect_view`, with its ViewLink); and the options the demos share, for
!> screenshots and tests (DemoOptions).
!>
!> Keys: VSG's key symbols are X11 keysyms. Letters are lower case without
!> Shift (the modified symbol upper case with it); other printable keys are
!> their character (the modified symbol the character typed, e.g. '!' for
!> Shift+1); the keypad's keys (Qt_KeypadModifier) are VSG's KEY_KP_*; the
!> editing, function (F1 to F20), modifier and lock keys their keysyms. On
!> macOS Qt reports Command as Control and Control as Meta, as VSG's own
!> macOS window does.
module vkframes_qt

    use qt
    use vkframes, only: wp, WindowProxy, MOUSE_MOVE, MOUSE_PRESS, MOUSE_RELEASE

    implicit none

    private

    public :: vsg_button, vsg_mask, vsg_key
    public :: forward_mouse, forward_wheel, forward_key
    public :: watch_surface, connect_view

    !> The options the demos share (for screenshots and tests), taken from the
    !> command line among each demo's own (`take`): `--capture file` saves
    !> frame `--capture-frame n` (default: the demo's, set before parsing) to
    !> the file; `--quit-after seconds` closes the main window then
    !> (`start_quit_timer`); `--recreate-surface n` destroys the window's
    !> native surface at frame n and creates it again (a test: the view
    !> closes, and opens again on the next expose). `after_frame` does the
    !> capture and the recreation.
    type, public :: DemoOptions
        character(len=:), allocatable :: capture            !! the file a frame is saved to (unallocated: none)
        integer :: capture_frame = 120                      !! the frame saved
        real(wp) :: quit_after = -1.0_wp                    !! seconds before quitting (< 0: never)
        integer :: recreate_frame = -1                      !! the frame the surface is made again at (< 0: never)
    contains
        procedure :: take => options_take
        procedure :: after_frame => options_after_frame
        procedure :: start_quit_timer => options_start_quit_timer
    end type DemoOptions

    !> A QWindow showing a WindowProxy (`connect_view`): what its connections
    !> forward to. It must stay alive while the window exists (e.g. a
    !> component of the program's saved state).
    type, public :: ViewLink
        type(WindowProxy), pointer :: proxy => null()       !! the proxy rendering into the window (not owned)
        type(QWindow), pointer :: window => null()          !! the window (not owned)
    end type ViewLink

    ! VSG's modifier bits (vsg::KeyModifier)
    integer, parameter :: MODKEY_SHIFT = 1, MODKEY_CONTROL = 4, MODKEY_ALT = 8, MODKEY_META = 128

    ! VSG's button masks: left (the trackball turns), middle (it pans)
    integer, parameter :: MASK_LEFT = 256, MASK_MIDDLE = 512

    !> A Shift + left drag is in progress: its left button is sent as the
    !> middle one (the trackball pans) until it is released.
    logical, save :: shift_panning = .false.

contains

    !> VSG's button number (1 left, 2 middle, 3 right; 0 other) for a Qt button.
    pure integer function vsg_button(qt_button)
        integer, intent(in) :: qt_button                    !! a Qt button (Qt_LeftButton, ...)
        select case (qt_button)
        case (Qt_LeftButton)
            vsg_button = 1
        case (Qt_MiddleButton)
            vsg_button = 2
        case (Qt_RightButton)
            vsg_button = 3
        case default
            vsg_button = 0
        end select
    end function vsg_button

    !> VSG's mask of the buttons held (256 left, 512 middle, 1024 right) for Qt's.
    pure integer function vsg_mask(qt_buttons)
        integer, intent(in) :: qt_buttons                   !! Qt's buttons held (an ior of Qt_*Button)
        vsg_mask = 0
        if (iand(qt_buttons, Qt_LeftButton) /= 0) vsg_mask = ior(vsg_mask, 256)
        if (iand(qt_buttons, Qt_MiddleButton) /= 0) vsg_mask = ior(vsg_mask, 512)
        if (iand(qt_buttons, Qt_RightButton) /= 0) vsg_mask = ior(vsg_mask, 1024)
    end function vsg_mask

    !> VSG's key symbol (without modifiers), modified symbol and modifier bits
    !> for a Qt key event; symbol 0 for a key VSG has no symbol for.
    subroutine vsg_key(event, symbol, modified_symbol, modifier)
        type(QKeyEvent), intent(in) :: event                !! the key event
        integer, intent(out) :: symbol                      !! VSG's key symbol without the modifiers (0: none)
        integer, intent(out) :: modified_symbol             !! with them (e.g. upper case with Shift)
        integer, intent(out) :: modifier                    !! VSG's modifier bits held
        character(len=:), allocatable :: text
        integer :: key, mods
        key = event%key()
        mods = event%modifiers()
        modifier = 0
        if (iand(mods, Qt_ShiftModifier) /= 0) modifier = ior(modifier, MODKEY_SHIFT)
        if (iand(mods, Qt_ControlModifier) /= 0) modifier = ior(modifier, MODKEY_CONTROL)
        if (iand(mods, Qt_AltModifier) /= 0) modifier = ior(modifier, MODKEY_ALT)
        if (iand(mods, Qt_MetaModifier) /= 0) modifier = ior(modifier, MODKEY_META)
        if (iand(mods, Qt_KeypadModifier) /= 0) then
            symbol = keypad_symbol(key)
            if (symbol /= 0) then
                modified_symbol = symbol
                return
            end if
        end if
        symbol = special_symbol(key)
        if (symbol /= 0) then
            modified_symbol = symbol
            return
        end if
        if (key >= 32 .and. key < 127) then                 ! a printable character (Qt's key: its code, letters upper case)
            symbol = key
            if (key >= Qt_Key_A .and. key <= Qt_Key_Z) symbol = key - Qt_Key_A + iachar('a')
            modified_symbol = symbol
            text = event%text()
            if (len(text) == 1) then
                if (iachar(text) >= 32 .and. iachar(text) < 127) modified_symbol = iachar(text)
            end if
            if (key >= Qt_Key_A .and. key <= Qt_Key_Z .and. iand(modifier, MODKEY_SHIFT) /= 0) modified_symbol = key
            return
        end if
        symbol = 0
        modified_symbol = 0
    end subroutine vsg_key

    !> VSG's keysym of a non-printing Qt key (0 if none).
    pure integer function special_symbol(key) result(symbol)
        integer, intent(in) :: key                          !! Qt's key (Qt_Key_*)
        select case (key)
        case (Qt_Key_Escape)
            symbol = int(z'FF1B')
        case (Qt_Key_Tab, Qt_Key_Backtab)                  ! (Backtab: Shift+Tab)
            symbol = int(z'FF09')
        case (Qt_Key_Backspace)
            symbol = int(z'FF08')
        case (Qt_Key_Return)
            symbol = int(z'FF0D')
        case (Qt_Key_Enter)
            symbol = int(z'FF8D')
        case (Qt_Key_Insert)
            symbol = int(z'FF63')
        case (Qt_Key_Delete)
            symbol = int(z'FFFF')
        case (Qt_Key_Pause)
            symbol = int(z'FF13')
        case (Qt_Key_Print)
            symbol = int(z'FF61')
        case (Qt_Key_SysReq)
            symbol = int(z'FF15')
        case (Qt_Key_Clear)
            symbol = int(z'FF0B')
        case (Qt_Key_Home)
            symbol = int(z'FF50')
        case (Qt_Key_End)
            symbol = int(z'FF57')
        case (Qt_Key_Left)
            symbol = int(z'FF51')
        case (Qt_Key_Up)
            symbol = int(z'FF52')
        case (Qt_Key_Right)
            symbol = int(z'FF53')
        case (Qt_Key_Down)
            symbol = int(z'FF54')
        case (Qt_Key_PageUp)
            symbol = int(z'FF55')
        case (Qt_Key_PageDown)
            symbol = int(z'FF56')
        case (Qt_Key_Shift)
            symbol = int(z'FFE1')
        case (Qt_Key_Control)
            symbol = int(z'FFE3')
        case (Qt_Key_Meta)
            symbol = int(z'FFE7')
        case (Qt_Key_Alt)
            symbol = int(z'FFE9')
        case (Qt_Key_CapsLock)
            symbol = int(z'FFE5')
        case (Qt_Key_NumLock)
            symbol = int(z'FF7F')
        case (Qt_Key_ScrollLock)
            symbol = int(z'FF14')
        case (Qt_Key_Menu)
            symbol = int(z'FF67')
        case (Qt_Key_Help)
            symbol = int(z'FF6A')
        case (Qt_Key_F1:Qt_Key_F20)                         ! (consecutive in both)
            symbol = int(z'FFBE') + (key - Qt_Key_F1)
        case default
            symbol = 0
        end select
    end function special_symbol

    !> VSG's keysym of a keypad key (0 if none: e.g. its arrows, as the others).
    pure integer function keypad_symbol(key) result(symbol)
        integer, intent(in) :: key                          !! Qt's key (Qt_Key_*), from the keypad
        select case (key)
        case (Qt_Key_0:Qt_Key_9)
            symbol = int(z'FFB0') + (key - Qt_Key_0)
        case (Qt_Key_Asterisk)
            symbol = int(z'FFAA')
        case (Qt_Key_Plus)
            symbol = int(z'FFAB')
        case (Qt_Key_Minus)
            symbol = int(z'FFAD')
        case (Qt_Key_Period)
            symbol = int(z'FFAE')
        case (Qt_Key_Slash)
            symbol = int(z'FFAF')
        case (Qt_Key_Equal)
            symbol = int(z'FFBD')
        case (Qt_Key_Enter)
            symbol = int(z'FF8D')
        case default
            symbol = 0
        end select
    end function keypad_symbol

    !> A Qt mouse event of `window` to `proxy` (kind MOUSE_MOVE, MOUSE_PRESS
    !> or MOUSE_RELEASE), in framebuffer pixels; accepted. A left drag started
    !> with Shift held pans: its left button is sent as the middle one (VSG's
    !> trackball: left turns, middle pans, right zooms) until it is released.
    subroutine forward_mouse(proxy, window, kind, event)
        type(WindowProxy), intent(inout) :: proxy           !! the view the event goes to
        type(QWindow), intent(in) :: window                 !! its window (its pixel ratio)
        integer, intent(in) :: kind                         !! MOUSE_MOVE, MOUSE_PRESS or MOUSE_RELEASE
        type(QMouseEvent), intent(in) :: event              !! the event (logical pixels)
        real(wp) :: dpr
        integer :: mask, button
        logical :: release_pan
        dpr = window%devicePixelRatio()
        mask = vsg_mask(event%buttons())
        button = vsg_button(event%button())
        release_pan = .false.
        if (kind == MOUSE_PRESS .and. button == 1 .and. iand(event%modifiers(), Qt_ShiftModifier) /= 0) &
            shift_panning = .true.
        if (shift_panning) then
            if (iand(mask, MASK_LEFT) /= 0) mask = ior(iand(mask, not(MASK_LEFT)), MASK_MIDDLE)
            if (button == 1) then
                button = 2
                release_pan = kind == MOUSE_RELEASE
            end if
        end if
        call proxy%mouse(kind, nint(event%x() * dpr), nint(event%y() * dpr), mask, button)
        if (release_pan) shift_panning = .false.
        call event%accept()
    end subroutine forward_mouse

    !> A Qt wheel event to `proxy`, in wheel notches (a trackpad's are
    !> fractions); accepted (ignored if it has no vertical part).
    subroutine forward_wheel(proxy, event)
        type(WindowProxy), intent(inout) :: proxy           !! the view the event goes to
        type(QWheelEvent), intent(in) :: event              !! the event (angle deltas: 120 a notch)
        if (event%angleDeltaY() == 0) return
        call proxy%scroll(event%angleDeltaY() / 120.0_wp)
        call event%accept()
    end subroutine forward_wheel

    !> A Qt key press (`pressed`) or release to `proxy`, as VSG's key symbols;
    !> accepted, unless VSG has no symbol for the key.
    subroutine forward_key(proxy, pressed, event)
        type(WindowProxy), intent(inout) :: proxy           !! the view the event goes to
        logical, intent(in) :: pressed                      !! pressed (else released)
        type(QKeyEvent), intent(in) :: event                !! the event
        integer :: symbol, modified_symbol, modifier
        call vsg_key(event, symbol, modified_symbol, modifier)
        if (symbol == 0) return
        call proxy%key(pressed, symbol, modified_symbol, modifier)
        call event%accept()
    end subroutine forward_key

    !> Connect a QWindow's events (`events`, its QtfEventFilter) to the
    !> WindowProxy rendering into it, through `link` (set here; it must stay
    !> alive while the window exists): the mouse (buttons, moves and the
    !> wheel) and, unless `keys` is false (a program handling the keys
    !> itself), the keys forwarded (forward_mouse, forward_wheel,
    !> forward_key); a resize given to the proxy (at the window's pixel
    !> ratio), with a frame requested; and, unless `surface` is false (a
    !> program doing more when it goes), the native surface watched
    !> (watch_surface). The window's expose and update requests are the
    !> program's (when to open the view and draw).
    subroutine connect_view(link, events, proxy, window, keys, surface)
        type(ViewLink), intent(inout), target :: link       !! the link the connections use
        type(QtfEventFilter), intent(in) :: events          !! the window's event filter
        type(WindowProxy), intent(inout), target :: proxy   !! the proxy rendering into the window
        type(QWindow), intent(in), target :: window         !! the window
        logical, intent(in), optional :: keys               !! the keys forwarded (default true)
        logical, intent(in), optional :: surface            !! the surface watched (default true)
        logical :: with_keys, with_surface
        with_keys = .true.
        if (present(keys)) with_keys = keys
        with_surface = .true.
        if (present(surface)) with_surface = surface
        link%proxy => proxy
        link%window => window
        call events%mousePressEvent%connect(on_view_press, link)
        call events%mouseReleaseEvent%connect(on_view_release, link)
        call events%mouseMoveEvent%connect(on_view_move, link)
        call events%wheelEvent%connect(on_view_wheel, link)
        if (with_keys) then
            call events%keyPressEvent%connect(on_view_key_press, link)
            call events%keyReleaseEvent%connect(on_view_key_release, link)
        end if
        call events%resizeEvent%connect(on_view_resize, link)
        if (with_surface) call watch_surface(events, proxy)
    end subroutine connect_view

    !> A mouse button pressed over a linked window: to its proxy.
    subroutine on_view_press(event, data)
        type(QMouseEvent), intent(in) :: event              !! the button and the pointer's position
        class(*), intent(inout), optional :: data           !! the ViewLink
        call forward_linked_mouse(MOUSE_PRESS, event, data)
    end subroutine on_view_press

    !> A mouse button released over a linked window: to its proxy.
    subroutine on_view_release(event, data)
        type(QMouseEvent), intent(in) :: event              !! the button and the pointer's position
        class(*), intent(inout), optional :: data           !! the ViewLink
        call forward_linked_mouse(MOUSE_RELEASE, event, data)
    end subroutine on_view_release

    !> The mouse moved over a linked window: to its proxy.
    subroutine on_view_move(event, data)
        type(QMouseEvent), intent(in) :: event              !! the pointer's position and buttons
        class(*), intent(inout), optional :: data           !! the ViewLink
        call forward_linked_mouse(MOUSE_MOVE, event, data)
    end subroutine on_view_move

    !> A mouse event of a linked window (`data`) to its proxy.
    subroutine forward_linked_mouse(kind, event, data)
        integer, intent(in) :: kind                         !! MOUSE_MOVE, MOUSE_PRESS or MOUSE_RELEASE
        type(QMouseEvent), intent(in) :: event              !! the event
        class(*), intent(inout), optional :: data           !! the ViewLink
        if (.not. present(data)) return
        select type (data)
        type is (ViewLink)
            call forward_mouse(data%proxy, data%window, kind, event)
        end select
    end subroutine forward_linked_mouse

    !> The wheel over a linked window: to its proxy.
    subroutine on_view_wheel(event, data)
        type(QWheelEvent), intent(in) :: event              !! the wheel's steps
        class(*), intent(inout), optional :: data           !! the ViewLink
        if (.not. present(data)) return
        select type (data)
        type is (ViewLink)
            call forward_wheel(data%proxy, event)
        end select
    end subroutine on_view_wheel

    !> A key pressed in a linked window: to its proxy.
    subroutine on_view_key_press(event, data)
        type(QKeyEvent), intent(in) :: event                !! the key and its modifiers
        class(*), intent(inout), optional :: data           !! the ViewLink
        if (.not. present(data)) return
        select type (data)
        type is (ViewLink)
            call forward_key(data%proxy, .true., event)
        end select
    end subroutine on_view_key_press

    !> A key released in a linked window: to its proxy.
    subroutine on_view_key_release(event, data)
        type(QKeyEvent), intent(in) :: event                !! the key and its modifiers
        class(*), intent(inout), optional :: data           !! the ViewLink
        if (.not. present(data)) return
        select type (data)
        type is (ViewLink)
            call forward_key(data%proxy, .false., event)
        end select
    end subroutine on_view_key_release

    !> A linked window resized: its proxy too (framebuffer pixels, at the
    !> window's pixel ratio), and a frame requested (drawn at the new size).
    subroutine on_view_resize(event, data)
        type(QResizeEvent), intent(in) :: event             !! the new size (logical pixels)
        class(*), intent(inout), optional :: data           !! the ViewLink
        real(wp) :: dpr
        if (.not. present(data)) return
        select type (data)
        type is (ViewLink)
            if (.not. data%proxy%is_open()) return
            dpr = data%window%devicePixelRatio()
            call data%proxy%resize(nint(event%width() * dpr), nint(event%height() * dpr), pixel_ratio=dpr)
            call data%window%requestUpdate()
        end select
    end subroutine on_view_resize

    !> Take command-line argument i if it is one of the shared options (see
    !> DemoOptions): its value read, and i moved past them both. False (i
    !> unchanged) for any other argument: the demo's own.
    logical function options_take(self, i) result(taken)
        class(DemoOptions), intent(inout) :: self
        integer, intent(inout) :: i                         !! the argument (the next one's after, if taken)
        character(len=1024) :: arg, value
        call get_command_argument(i, arg)
        call get_command_argument(i + 1, value)
        taken = .true.
        select case (arg)
        case ('--capture')
            self%capture = trim(value)
        case ('--capture-frame')
            read(value, *) self%capture_frame
        case ('--quit-after')
            read(value, *) self%quit_after
        case ('--recreate-surface')
            read(value, *) self%recreate_frame
        case default
            taken = .false.
            return
        end select
        i = i + 2
    end function options_take

    !> After frame `frame` of `proxy` in `window`: the frame saved if it is
    !> --capture-frame (window `index` > 1 of several: to the file name with
    !> _index before its extension); the native surface destroyed and created
    !> again if it is --recreate-surface. True if the surface was: the caller
    !> asks for no next frame (the expose that follows opens the view again).
    logical function options_after_frame(self, proxy, window, frame, index) result(recreated)
        class(DemoOptions), intent(in) :: self
        type(WindowProxy), intent(inout) :: proxy           !! the view drawn
        type(QWindow), intent(inout) :: window              !! its window
        integer, intent(in) :: frame                        !! the frame just drawn (from 1)
        integer, intent(in), optional :: index              !! the window's number (default 1)
        character(len=:), allocatable :: file
        character(len=16) :: suffix
        integer :: dot
        recreated = .false.
        if (allocated(self%capture) .and. frame == self%capture_frame) then
            file = self%capture
            if (present(index)) then
                if (index > 1) then                         ! file_2.png, ...
                    dot = scan(file, '.', back=.true.)
                    if (dot == 0) dot = len(file) + 1
                    write(suffix, '(a,i0)') '_', index
                    file = file(1:dot-1)//trim(suffix)//file(dot:)
                end if
            end if
            if (proxy%capture(file)) write(*, '(a)') 'captured '//file
        end if
        if (frame == self%recreate_frame) then             ! (a test: Qt destroys and remakes the native window)
            write(*, '(a,i0)') 'native window destroyed and created again at frame ', frame
            call window%destroy()
            call window%create()
            call window%show()
            recreated = .true.
        end if
    end function options_after_frame

    !> With --quit-after: `main` closed after that many seconds (it must stay
    !> alive until then).
    subroutine options_start_quit_timer(self, main)
        class(DemoOptions), intent(in) :: self
        type(QMainWindow), intent(inout), target :: main    !! the window closed
        if (self%quit_after > 0.0_wp) call QTimer_singleShot(nint(1000.0_wp * self%quit_after), on_quit_timer, main)
    end subroutine options_start_quit_timer

    !> The --quit-after time is up: the main window (`data`) closed.
    subroutine on_quit_timer(data)
        class(*), intent(inout), optional :: data           !! the QMainWindow
        if (.not. present(data)) return
        select type (data)
        type is (QMainWindow)
            call data%close()
        end select
    end subroutine on_quit_timer

    !> Close `proxy` when the window `events` filters is about to lose its
    !> native surface (Qt destroying it: re-parenting, `destroy`, the window
    !> deleted), while the surface still exists. The host opens it again on
    !> the window's next expose (as at first), with the proxy's views, HUD
    !> texts and scene as before.
    subroutine watch_surface(events, proxy)
        type(QtfEventFilter), intent(in) :: events          !! the window's event filter
        type(WindowProxy), intent(inout), target :: proxy   !! the view closed then (it must stay alive)
        call events%surfaceAboutToBeDestroyed%connect(on_surface_destroyed, proxy)
    end subroutine watch_surface

    !> The surface is going: stop rendering into it (waiting for the GPU).
    subroutine on_surface_destroyed(data)
        class(*), intent(inout), optional :: data           !! the WindowProxy
        select type (data)
        type is (WindowProxy)
            call data%close()
        end select
    end subroutine on_surface_destroyed

end module vkframes_qt
