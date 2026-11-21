!*****************************************************************************************
!>
!  Test of `text_buffer_module`'s number formatting: numbers over a wide range of
!  magnitudes are formatted, read back, and compared to the original values.
!
!  `fpm test text_buffer_test --profile release`

    program text_buffer_test

    use halo_kinds_module, only: wp
    use text_buffer_module

    implicit none

    integer,parameter :: n = 1000000
    type(text_buffer) :: buf
    real(wp) :: x, y, u, err, worst
    real(wp),dimension(*),parameter :: special = [0.0_wp, 1.0_wp, -1.0_wp, 0.1_wp, 1.5_wp, &
                                                  9.9999999999999999e5_wp, 1.0e16_wp, 1.0e17_wp, &
                                                  820319112.03980994_wp, -1738.0_wp, 1.0e-300_wp, &
                                                  -2.5e300_wp, tiny(1.0_wp), huge(1.0_wp), &
                                                  epsilon(1.0_wp), 123456789012345678.0_wp]
    integer :: i, n_exact, n_ok
    logical :: all_ok

    all_ok = .true.
    worst = 0.0_wp
    n_exact = 0
    n_ok = 0
    do i = 1, size(special) + n
        if (i <= size(special)) then
            x = special(i)
        else
            ! random numbers, with random signs and magnitudes 1e-12 .. 1e12:
            call random_number(u); x = 1.0_wp + u
            call random_number(u); x = x * 10.0_wp**int(24*u - 12)
            call random_number(u); if (u < 0.5_wp) x = -x
        end if
        call buf%clear()
        call buf%add_real(x)
        read(buf%s(1:buf%n),*) y
        if (x == 0.0_wp) then
            err = abs(y)
        else
            err = abs(y - x) / abs(x)
        end if
        worst = max(worst, err)
        if (y == x) n_exact = n_exact + 1
        if (i <= size(special)) write(*,'(A30,A,A)') buf%s(1:buf%n), '   ', merge('exact', '     ', y==x)
    end do
    write(*,'(A,ES10.2)') 'max relative error:            ', worst
    write(*,'(A,F8.4,A)') 'read back as the same number:  ', 100.0_wp*n_exact/(size(special)+n), ' %'
    if (worst > 0.0_wp .or. n_exact /= size(special)+n) then
        write(*,'(A)') 'FAILED: not all numbers read back exactly'
        all_ok = .false.
    end if

    ! JSON syntax checks:
    call buf%clear(); call buf%add_real_array([1.0_wp, -2.5_wp]); call buf%add_str(' ')
    call buf%add_int(-42); call buf%add_str(' '); call buf%add_int(0)
    write(*,'(A)') buf%s(1:buf%n)
    if (buf%s(1:buf%n) /= '[1.0E+0, -2.5E+0] -42 0') then
        write(*,'(A)') 'FAILED: unexpected array/integer text'
        all_ok = .false.
    end if

    if (.not. all_ok) error stop 'text_buffer_test FAILED'
    write(*,'(A)') 'text_buffer_test passed'

    end program text_buffer_test
!*****************************************************************************************
