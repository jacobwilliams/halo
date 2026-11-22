!*****************************************************************************************
!>
!  The sensitivities of the constraints to each kind of variable (the technical
!  note's sensitivity table): the norm of each column of the scaled analytic
!  Jacobian at the initial guess, per unit change of the (unscaled) variable,
!  with the median over the nodes for each kind of variable (epoch per hour,
!  position per km, velocity per m/s) and each node position in the rev
!  (periapsis, quarter, apoapsis, three-quarter).
!
!  Usage (from the halo root directory):
!
!    fpm test sensitivity_test --profile release --flag "-fopenmp" -- [config.json]
!
!  The default config file is `docs/halo-paper/runs/baseline_stm.json`. The problem is set up
!  with `jacobian_mode=2` and the splined ephemeris.

    program sensitivity_test

    use halo_module
    use halo_kinds_module, only: wp
    use iso_fortran_env, only: output_unit

    implicit none

    character(len=*),dimension(3),parameter :: type_name = ['epoch (per hour)  ', &
                                                            'position (per km) ', &
                                                            'velocity (per m/s)']
    character(len=*),dimension(0:3),parameter :: node_name = ['periapsis    ', 'quarter      ', &
                                                              'apoapsis     ', 'three-quarter']

    type(my_solver_type) :: solver
    real(wp),dimension(:),allocatable :: x, jac, colnorm
    integer,dimension(:),allocatable :: irow, icol
    real(wp),dimension(:,:,:),allocatable :: vals  !! (n, variable type 1..3, node position 0..3)
    integer,dimension(3,0:3) :: cnt
    character(len=1024) :: arg
    character(len=:),allocatable :: config
    character(len=:),allocatable :: label
    integer :: k, n, j, iseg, ityp, nodepos

    config = 'docs/halo-paper/runs/baseline_stm.json'
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        config = trim(arg)
    end if

    call solver%read_config_file(config)
    solver%mission%jacobian_mode = 2
    solver%mission%use_splined_ephemeris = .true.
    call solver%mission%init(x)
    call solver%mission%compute_jacobian_stm(x, jac)
    call solver%mission%get_sparsity_pattern(irow, icol)

    ! column norms, per unit change of the unscaled variable:
    n = size(x)
    allocate(colnorm(n))
    colnorm = 0.0_wp
    do k = 1, size(jac)
        colnorm(icol(k)) = colnorm(icol(k)) + jac(k)**2
    end do
    colnorm = sqrt(colnorm) / solver%mission%xscale

    allocate(vals(n,3,0:3))
    cnt = 0
    do j = 1, n
        label = trim(solver%mission%xname(j))  ! e.g. 'SEG8 T0 (day)'
        read(label(4:index(label,' ')-1),*) iseg
        ! the node of this variable: segment 1 starts at node 1, segment 2k at node k+1 (4 nodes per rev)
        if (iseg == 1) then
            nodepos = 0
        else
            nodepos = mod(iseg/2, 4)
        end if
        if (index(label,'T0') > 0) then
            ityp = 1;  colnorm(j) = colnorm(j) / 24.0_wp     ! per hour
        else if (index(label,' R') > 0) then
            ityp = 2                                          ! per km
        else
            ityp = 3;  colnorm(j) = colnorm(j) / 1000.0_wp   ! per m/s
        end if
        cnt(ityp,nodepos) = cnt(ityp,nodepos) + 1
        vals(cnt(ityp,nodepos),ityp,nodepos) = colnorm(j)
    end do

    write(output_unit,'(A)') ''
    write(output_unit,'(A)') ' median column norm of the scaled Jacobian (per unit change), by node position:'
    write(output_unit,'(A20,4A15)') '', node_name
    do ityp = 1, 3
        write(output_unit,'(A20,4ES15.3)') type_name(ityp), &
            (median(vals(1:cnt(ityp,nodepos),ityp,nodepos)), nodepos=0,3)
    end do

    contains

        real(wp) function median(a)
        real(wp),dimension(:),intent(in) :: a
        real(wp),dimension(:),allocatable :: b
        real(wp) :: t
        integer :: i, m
        b = a
        do i = 2, size(b)  ! (insertion sort)
            t = b(i)
            m = i - 1
            do while (m >= 1)
                if (b(m) <= t) exit
                b(m+1) = b(m)
                m = m - 1
            end do
            b(m+1) = t
        end do
        median = b((size(b)+1)/2)
        end function median

    end program sensitivity_test
!*****************************************************************************************
