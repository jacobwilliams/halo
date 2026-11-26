!*****************************************************************************************
!>
!  Experiment: the basis choice of sqpopt's sparse QP solver (`choose_basis`, a rank-revealing
!  LUSOL factorization) on a Jacobian written by `jacobian_dump_test`.
!
!  Usage:
!
!    fpm test basis_choice_test --profile release -- jac_file
!
!  For several variants, the matrix that `choose_basis` gives LUSOL is built (the rows
!  normalized by the norms of their variable columns, each variable column weighted 1, and,
!  for the equality rows, a slack column with entry `-w/|row|`), factorized by `lu1fac`
!  as `independent_columns` does, and the rank found, the number of columns picked, and the
!  number of rows without a pivot (which a basis repair would cover with their slacks) are
!  printed.

    program basis_choice_test

    use lusol,           only: lu1fac
    use lusol_precision, only: ip, rp
    use iso_fortran_env, only: int32, real64, output_unit

    implicit none

    integer(int32) :: m, n, nnz
    integer(int32),dimension(:),allocatable :: irow, icol, keep
    real(real64),dimension(:),allocatable :: val, rn
    character(len=1024) :: jac_file
    integer :: iunit, k

    if (command_argument_count() < 1) then
        ! (a tool, run by hand: under a plain `fpm test` there are no arguments, and it is skipped)
        write(output_unit,'(A)') 'usage: basis_choice_test jac_file (no arguments: skipped)'
        stop
    end if
    call get_command_argument(1, jac_file)
    open(newunit=iunit, file=trim(jac_file), access='stream', form='unformatted', status='old')
    read(iunit) m, n, nnz
    allocate(irow(nnz), icol(nnz), val(nnz), keep(n))
    read(iunit) irow, icol
    read(iunit) val
    read(iunit) keep
    close(iunit)
    write(output_unit,'(A,I0,A,I0,A,I0)') 'J: ', m, ' x ', n, ', nnz ', nnz

    allocate(rn(m))
    rn = 0.0_real64
    do k = 1, nnz
        rn(irow(k)) = rn(irow(k)) + val(k)**2
    end do
    rn = sqrt(rn)

    write(output_unit,'(A)') 'variant                                    pivot   Utol2     rank  picked  rows w/o pivot  time(s)'
    ! as choose_basis: TCP, rel 1e-8, equality slacks weighted 1e-6
    call try('choose_basis (TCP, 1e-8, slack w 1e-6)',  2, 1.0e-8_real64,  .true., 1.0e-6_real64, .true.)
    call try('  same, slack w 1e-3',                     2, 1.0e-8_real64,  .true., 1.0e-3_real64, .true.)
    call try('  same, slack w 1',                        2, 1.0e-8_real64,  .true., 1.0_real64,    .true.)
    call try('  slack w 1e-6, Utol2 1e-12',              2, 1.0e-12_real64, .true., 1.0e-6_real64, .true.)
    call try('variables only (TCP, 1e-8)',               2, 1.0e-8_real64,  .false., 0.0_real64,   .true.)
    call try('variables only, rows not normalized',      2, 1.0e-8_real64,  .false., 0.0_real64,   .false.)
    call try('variables only, TRP',                      1, 1.0e-8_real64,  .false., 0.0_real64,   .true.)
    call try('variables only, TPP',                      0, 1.0e-8_real64,  .false., 0.0_real64,   .true.)
    call try('variables only, TCP, Utol2 1e-11',         2, 1.0e-11_real64, .false., 0.0_real64,   .true.)
    call try('slack w 1e-6, TRP',                        1, 1.0e-8_real64,  .true., 1.0e-6_real64, .true.)

    contains

    subroutine try(name, pivot, utol2, slacks, w, normalize)
    !! factorize one variant and print what LUSOL found
    character(len=*),intent(in) :: name   !! a description
    integer,intent(in)          :: pivot  !! LUSOL's luparm(6): 0 TPP, 1 TRP, 2 TCP
    real(real64),intent(in)     :: utol2  !! LUSOL's relative singularity tolerance
    logical,intent(in)          :: slacks !! include a slack column for each row
    real(real64),intent(in)     :: w      !! the slacks' weight
    logical,intent(in)          :: normalize !! normalize the rows by their variable columns' norms

    integer(ip) :: mm, nn, nelem, lena, inform, luparm(30), nrank, j
    real(rp) :: parmlu(30)
    real(rp),dimension(:),allocatable :: a, wk
    integer(ip),dimension(:),allocatable :: indc, indr, p, q, lenc, lenr, locc, locr, &
                                            iploc, iqloc, ipinv, iqinv
    logical,dimension(:),allocatable :: covered
    integer :: c0, c1, rate, attempt

    call system_clock(c0, rate)
    mm = m
    nn = n
    if (slacks) nn = n + m
    nelem = nnz
    if (slacks) nelem = nnz + m
    lena = 1 + max(5*nelem, 10*mm, 10*nn, 10000_ip)
    allocate(p(mm), q(nn), lenc(nn), lenr(mm), locc(nn), locr(mm), iploc(nn), iqloc(mm), &
             ipinv(mm), iqinv(nn), wk(nn))
    do attempt = 1, 3
        if (allocated(a)) deallocate(a, indc, indr)
        allocate(a(lena), indc(lena), indr(lena))
        do k = 1, nnz
            a(k) = val(k)
            if (normalize) a(k) = a(k) / rn(irow(k))
            indc(k) = irow(k)
            indr(k) = icol(k)
        end do
        if (slacks) then
            do k = 1, m
                a(nnz+k) = -w
                if (normalize) a(nnz+k) = -w / rn(k)
                indc(nnz+k) = k
                indr(nnz+k) = n + k
            end do
        end if
        luparm = 0
        luparm(1) = 6
        luparm(2) = -1
        luparm(3) = 5
        luparm(6) = pivot
        luparm(8) = 1
        parmlu = 0.0_rp
        parmlu(1) = 5.0_rp
        parmlu(2) = 5.0_rp
        parmlu(3) = epsilon(1.0_rp)**0.8_rp
        parmlu(4) = 1.0e-6_rp*epsilon(1.0_rp)**0.67_rp
        parmlu(5) = utol2
        parmlu(6) = 3.0_rp
        parmlu(7) = 0.3_rp
        parmlu(8) = 0.5_rp
        call lu1fac(mm, nn, nelem, lena, luparm, parmlu, a, indc, indr, p, q, &
                    lenc, lenr, locc, locr, iploc, iqloc, ipinv, iqinv, wk, inform)
        if (inform /= 7) exit
        lena = max(2*lena, luparm(13) + 1)
    end do
    nrank = luparm(16)
    allocate(covered(mm))
    covered = .false.
    do j = 1, nrank
        covered(p(j)) = .true.
    end do
    call system_clock(c1)
    block
        logical,dimension(:),allocatable :: pivcol
        allocate(pivcol(nn)); pivcol = .false.
        pivcol(q(1:nrank)) = .true.
        write(output_unit,'(A42,I4,ES10.1,I8,I8,I12,F10.2,A,I0,A,I0,A,I0)') name, pivot, utol2, nrank, count(wk > 0.0_rp), &
            count(.not. covered), real(c1-c0)/real(rate), '   inform ', inform, &
            '  | pivot cols with w<=0: ', count(pivcol .and. wk <= 0.0_rp), ', w>0 non-pivot: ', count(.not. pivcol .and. wk > 0.0_rp)
    end block

    end subroutine try

    end program basis_choice_test
!*****************************************************************************************
