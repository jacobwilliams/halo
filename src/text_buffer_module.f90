!*****************************************************************************************
!>
!  A growable text buffer, with fast number formatting, for writing large
!  JSON files directly (see [[export_trajectory_json_file]]).
!
!  The numbers are formatted with integer arithmetic instead of formatted
!  `write` statements. This is much faster, and it can run in parallel
!  (gfortran's formatted I/O is serialized by a lock, even for internal
!  writes). The numbers read back exactly (see [[add_real]]).

    module text_buffer_module

    use halo_kinds_module, only: wp
    use iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    private

    type,public :: text_buffer
        !! a text buffer. text is appended with the `add` methods, and the
        !! result is `buf%s(1:buf%n)`. the storage is on the heap, and grows
        !! (doubles) as needed.
        character(len=:),allocatable :: s  !! the storage
        integer :: n = 0  !! number of characters used
    contains
        procedure,public :: reserve   !! make room for this many more characters
        procedure,public :: clear     !! empty it (keeping the storage)
        procedure,public :: add_str   !! append a string
        procedure,public :: add_int   !! append an integer
        procedure,public :: add_real  !! append a real number (JSON syntax)
        procedure,public :: add_real_array  !! append a JSON array of real numbers
    end type text_buffer

    real(wp),dimension(0:22),parameter :: p10 = [1.0e0_wp,  1.0e1_wp,  1.0e2_wp,  1.0e3_wp,  1.0e4_wp,  &
                                                 1.0e5_wp,  1.0e6_wp,  1.0e7_wp,  1.0e8_wp,  1.0e9_wp,  &
                                                 1.0e10_wp, 1.0e11_wp, 1.0e12_wp, 1.0e13_wp, 1.0e14_wp, &
                                                 1.0e15_wp, 1.0e16_wp, 1.0e17_wp, 1.0e18_wp, 1.0e19_wp, &
                                                 1.0e20_wp, 1.0e21_wp, 1.0e22_wp]
        !! the powers of ten that are exactly representable in double precision

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Make sure the buffer has room for `nchars` more characters.

    pure subroutine reserve(me, nchars)

    class(text_buffer),intent(inout) :: me
    integer,intent(in) :: nchars

    character(len=:),allocatable :: tmp
    integer :: new_size

    if (.not. allocated(me%s)) then
        allocate(character(len=max(nchars, 1024)) :: me%s)
        me%n = 0
    else if (me%n + nchars > len(me%s)) then
        new_size = max(2*len(me%s), me%n + nchars)
        allocate(character(len=new_size) :: tmp)
        tmp(1:me%n) = me%s(1:me%n)
        call move_alloc(tmp, me%s)
    end if

    end subroutine reserve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Empty the buffer (the storage is kept, to be reused).

    pure subroutine clear(me)
    class(text_buffer),intent(inout) :: me
    me%n = 0
    end subroutine clear
!*****************************************************************************************

!*****************************************************************************************
!>
!  Append a string.

    pure subroutine add_str(me, str)

    class(text_buffer),intent(inout) :: me
    character(len=*),intent(in) :: str

    call me%reserve(len(str))
    me%s(me%n+1:me%n+len(str)) = str
    me%n = me%n + len(str)

    end subroutine add_str
!*****************************************************************************************

!*****************************************************************************************
!>
!  Append an integer.

    pure subroutine add_int(me, i)

    class(text_buffer),intent(inout) :: me
    integer,intent(in) :: i

    call add_digits(me, int(i,int64), negative=i<0)

    end subroutine add_int
!*****************************************************************************************

!*****************************************************************************************
!>
!  Append a real number, in JSON syntax, with 17 significant digits
!  (e.g., `8.2031911203980994E+8`; trailing zeros are removed, e.g. `1.5E+0`).
!
!  The 17 digits are correctly rounded, so the number reads back as exactly the
!  same double precision value. |x| is scaled to a 17-digit integer by an exact
!  power of ten (\( 10^k, |k| \le 22 \)), and the rounding error of that product
!  (or quotient) is computed exactly (Dekker's two-product), so the integer is the
!  correctly rounded one. Numbers too small or too large for this
!  (|x| < 1e-6 or |x| >= 1e39) are formatted with a `write` statement instead.
!  Non-finite numbers are written as `null`, since JSON has no NaN or Infinity.

    pure subroutine add_real(me, x)

    class(text_buffer),intent(inout) :: me
    real(wp),intent(in) :: x

    integer(int64),parameter :: ten16 = 10000000000000000_int64  !! 10^16
    integer(int64),parameter :: ten17 = 10_int64*ten16           !! 10^17

    real(wp) :: ax
    integer :: e10, i, last
    integer(int64) :: im
    logical :: ok
    character(len=17) :: digits
    character(len=30) :: str

    if (.not. ieee_is_finite(x)) then
        call me%add_str('null')
        return
    end if
    if (x == 0.0_wp) then
        call me%add_str('0.0')
        return
    end if

    ! im = |x| * 10^(16-e10), rounded, should have 17 digits:
    ax  = abs(x)
    e10 = floor(log10(ax))
    call scaled(ax, 16-e10, im, ok)
    ! (log10 may be off by one near powers of ten)
    if (ok .and. im >= ten17) then
        e10 = e10 + 1
        call scaled(ax, 16-e10, im, ok)
    else if (ok .and. im < ten16) then
        e10 = e10 - 1
        call scaled(ax, 16-e10, im, ok)
    end if
    if (ok .and. im == ten17) then
        ! rounded up to the next power of ten (e.g., 9.99999999999999999)
        im  = ten16
        e10 = e10 + 1
    end if
    if (.not. ok) then
        ! out of the range of the exact powers of ten:
        write(str,'(ES24.16E3)') x
        call me%add_str(trim(adjustl(str)))
        return
    end if

    ! the 17 digits:
    do i = 17, 1, -1
        digits(i:i) = achar(iachar('0') + int(mod(im, 10_int64)))
        im = im / 10_int64
    end do
    ! remove the trailing zeros (keeping at least one digit after the decimal point):
    last = 17
    do while (last > 2 .and. digits(last:last) == '0')
        last = last - 1
    end do

    call me%reserve(30)
    if (x < 0.0_wp) call me%add_str('-')
    call me%add_str(digits(1:1)//'.'//digits(2:last)//'E')
    if (e10 >= 0) call me%add_str('+')
    call add_digits(me, int(e10,int64), negative=e10<0)

    end subroutine add_real
!*****************************************************************************************

!*****************************************************************************************
!>
!  \( a \cdot 10^k \) (\( a > 0 \)), correctly rounded to an integer, when the
!  result is about \( 10^{16} \) to \( 10^{17} \). `ok` is false if
!  \( |k| > 22 \) (no exact power of ten).
!
!  Since doubles of that size are integers (their spacing is 2 or more), the
!  rounded product or quotient is an integer, and only its exact rounding error
!  has to be added (rounded).

    pure subroutine scaled(a, k, im, ok)

    real(wp),intent(in) :: a
    integer,intent(in) :: k
    integer(int64),intent(out) :: im
    logical,intent(out) :: ok

    real(wp) :: p, q, hi, lo, r

    im = 0
    ok = abs(k) <= 22
    if (.not. ok) return

    if (k >= 0) then
        ! a*p = hi + lo exactly:
        p = p10(k)
        call two_prod(a, p, hi, lo)
        im = int(hi, int64) + nint(lo, int64)
    else
        ! a/p = q + r/p, with the remainder r = a - q*p computed exactly:
        p = p10(-k)
        q = a / p
        call two_prod(q, p, hi, lo)
        r = (a - hi) - lo
        im = int(q, int64) + nint(r/p, int64)
    end if

    end subroutine scaled
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dekker's exact product: \( a b = p + e \) exactly (barring overflow), without FMA.

    pure subroutine two_prod(a, b, p, e)

    real(wp),intent(in) :: a, b
    real(wp),intent(out) :: p, e

    real(wp) :: ah, al, bh, bl

    p = a * b
    call split(a, ah, al)
    call split(b, bh, bl)
    e = ((ah*bh - p) + ah*bl + al*bh) + al*bl

    end subroutine two_prod
!*****************************************************************************************

!*****************************************************************************************
!>
!  Veltkamp's split of a double into two halves of 26 bits: \( a = h + l \).

    pure subroutine split(a, h, l)

    real(wp),intent(in) :: a
    real(wp),intent(out) :: h, l

    real(wp),parameter :: factor = 134217729.0_wp  !! 2^27 + 1
    real(wp) :: c

    c = factor * a
    h = c - (c - a)
    l = a - h

    end subroutine split
!*****************************************************************************************

!*****************************************************************************************
!>
!  Append a JSON array of real numbers (`[x1, x2, ...]`).

    pure subroutine add_real_array(me, x)

    class(text_buffer),intent(inout) :: me
    real(wp),dimension(:),intent(in) :: x

    integer :: i

    call me%reserve(2 + 26*size(x))  ! (at most 24 characters per number, plus ", ")
    call me%add_str('[')
    do i = 1, size(x)
        if (i > 1) call me%add_str(', ')
        call me%add_real(x(i))
    end do
    call me%add_str(']')

    end subroutine add_real_array
!*****************************************************************************************

!*****************************************************************************************
!>
!  Append the decimal digits of an integer.

    pure subroutine add_digits(me, i, negative)

    class(text_buffer),intent(inout) :: me
    integer(int64),intent(in) :: i
    logical,intent(in) :: negative

    character(len=20) :: str
    integer(int64) :: j
    integer :: k

    j = abs(i)
    k = len(str) + 1
    do
        k = k - 1
        str(k:k) = achar(iachar('0') + int(mod(j, 10_int64)))
        j = j / 10_int64
        if (j == 0) exit
    end do
    if (negative) call me%add_str('-')
    call me%add_str(str(k:))

    end subroutine add_digits
!*****************************************************************************************

    end module text_buffer_module
!*****************************************************************************************
