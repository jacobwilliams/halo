!*****************************************************************************************
!>
!  Test of `splined_moon_frame_module` (the MOON_PA frame with shared spline
!  coefficients): its rotation matrices must be identical to those of the
!  `moonframe` library's `moon_frame_interpolater` (a dev dependency, used only
!  here), for the same table and time span.
!
!  Checked: many epochs over a 20-year span (random, so the interval caches are
!  exercised in a scrambled order, plus the ends of the span), and copies of the
!  object (one per segment, as the solver makes them) evaluated in parallel.
!
!  Run from the halo root directory (it reads `data/moon_pa_2000_2100.csv`):
!
!    fpm test moon_frame_test --profile release --flag "-fopenmp"

    program moon_frame_test

    use halo_kinds_module, only: wp
    use splined_moon_frame_module, only: moon_pa_frame
    use moon_frame_module, only: moon_frame_interpolater

    implicit none

    character(len=*),parameter :: file = 'data/moon_pa_2000_2100.csv'
    real(wp),parameter :: day = 86400.0_wp
    real(wp),parameter :: et0 = 820303200.0_wp - 2.0_wp*6.5624_wp*day     !! the 20-year baseline's span
    real(wp),parameter :: etf = 820303200.0_wp + 1122.0_wp*6.5624_wp*day
    integer,parameter :: n = 200000   !! number of test epochs
    integer,parameter :: n_copies = 16  !! number of copies (parallel)

    type(moon_frame_interpolater) :: ref
    type(moon_pa_frame) :: new
    type(moon_pa_frame),dimension(:),allocatable :: copies
    real(wp),dimension(:),allocatable :: et
    real(wp),dimension(:,:,:),allocatable :: r_ref, r_par  !! (3,3,n)
    real(wp) :: u, worst, worst_par, worst_ortho
    integer :: i, j
    logical :: all_ok

    call ref%initialize(file, et0=et0, etf=etf)
    call new%initialize(file, et0=et0, etf=etf)

    ! test epochs: random, plus the ends of the span
    allocate(et(n))
    do i = 1, n
        call random_number(u)
        et(i) = et0 + u*(etf - et0)
    end do
    et(1:4) = [et0, etf, et0 + 1.0_wp, etf - 1.0_wp]

    ! serial: the new module vs. moonframe
    allocate(r_ref(3,3,n), r_par(3,3,n))
    worst = 0.0_wp
    worst_ortho = 0.0_wp
    do i = 1, n
        r_ref(:,:,i) = ref%j2000_to_frame(et(i))
        r_par(:,:,i) = new%j2000_to_frame(et(i))
        worst = max(worst, maxval(abs(r_par(:,:,i) - r_ref(:,:,i))))
        worst_ortho = max(worst_ortho, maxval(abs(matmul(r_par(:,:,i), transpose(r_par(:,:,i))) - identity())))
    end do

    ! parallel: copies of the object (as each segment has), each used for every
    ! n_copies-th epoch in its own thread
    allocate(copies(n_copies), source=new)
    r_par = 0.0_wp
    !$OMP PARALLEL DO PRIVATE(i)
    do j = 1, n_copies
        do i = j, n, n_copies
            r_par(:,:,i) = copies(j)%j2000_to_frame(et(i))
        end do
    end do
    !$OMP END PARALLEL DO
    worst_par = maxval(abs(r_par - r_ref))

    write(*,'(A,I0,A)') ' MOON_PA frame: ', n, ' epochs over the 20-year span'
    write(*,'(A,ES10.2)') '   max |R_new - R_moonframe|:             ', worst
    write(*,'(A,ES10.2)') '   max |R_new - R_moonframe| (parallel):  ', worst_par
    write(*,'(A,ES10.2)') '   max |R R^T - I|:                       ', worst_ortho

    all_ok = worst == 0.0_wp .and. worst_par == 0.0_wp .and. worst_ortho < 1.0e-14_wp
    if (.not. all_ok) error stop 'moon_frame_test FAILED'
    write(*,'(A)') ' moon_frame_test passed (identical to moonframe)'

    contains

        pure function identity() result(x)
        real(wp),dimension(3,3) :: x
        x = 0.0_wp
        x(1,1) = 1.0_wp; x(2,2) = 1.0_wp; x(3,3) = 1.0_wp
        end function identity

    end program moon_frame_test
!*****************************************************************************************
