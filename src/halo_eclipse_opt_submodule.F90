!*****************************************************************************************
!>
!  Eclipse optimization: moves a converged solution out of the Earth's shadow.
!
!  Starting from the converged solution `x0`, it finds the closest continuous trajectory
!  (in halo's scaled variables) that stays out of the shadow near the eclipses:
!
!  * minimize \( \frac{1}{2} \lVert \mathbf{x} - \mathbf{x}_0 \rVert^2 \)
!  * subject to the continuity constraints \( \mathbf{c}(\mathbf{x}) = \mathbf{0} \) (halo's)
!  * and, for each segment that comes near the shadow, \( \mathrm{KS}(\phi_k) \ge \phi_{min} \)
!
!  with the SQP solver [sqpopt](https://github.com/jacobwilliams/sqpopt). \( \phi \) is the
!  cubic shadow model (Williams et al., 2023), which is smooth inside and outside the
!  shadow. \( \phi_k \) are its values at samples every `eclipse_opt_sample_step` seconds
!  from the segment's start, across the part of the segment near the shadow, and KS is
!  their Kreisselmeier-Steinhauser aggregate
!  \( \mathrm{KS} = \phi_{lo} - \frac{1}{\rho} \ln \sum_k e^{-\rho(\phi_k - \phi_{lo})} \),
!  \( \phi_{lo} = \min_k \phi_k \): a smooth lower bound of the smallest sample (within
!  \( \ln(N)/\rho \) of it), so the constraint keeps every sample at or above
!  \( \phi_{min} \). One row per segment, not one per sample: the samples of a segment all
!  depend on the same few variables (its start node), and a row for each would make the
!  QP's working sets dependent (sqpopt's direct method then cycles).
!
!  Following the 2023 paper's process, the constraints are added only where the trajectory
!  comes near the shadow, and the trajectory is checked again after each optimization (with
!  new rows where needed). Each pass ends with a Newton solve of the continuity constraints
!  and the eclipse rows at their bound, as equalities (`eclipse_opt_newton_finish`), to
!  bring the defects back to the usual level without undoing the eclipse fix.
!
!  The continuity rows of the Jacobian are halo's analytic ones. The eclipse rows (7
!  nonzeros each: the epoch and state of the segment's start node; the samples are at fixed
!  offsets from the start) are computed by central differences, propagating only that
!  segment.
!
!  See `docs/eclipse_optimization_plan.md`.

    submodule (halo_module) halo_eclipse_opt_submodule

    use sqpopt_module,           only: sqpopt_type
    use sqpopt_problem_module,   only: sqpopt_problem_type
    use sqpopt_options_module,   only: sqpopt_options_type
    use sqpopt_qp_solver_module, only: sqpopt_qp_solver_type
    use sqpopt_types_module,     only: sqpopt_infinity, sqpopt_results_type
    use iso_fortran_env,         only: int64

    implicit none

    ! for each of the 8 segments of a rev: its start node, and the node whose epoch sets its final
    ! epoch (as in [[compute_jacobian_stm]])
    integer,dimension(8),parameter :: start_node = [1,2,2,3,3,4,4,5]
    integer,dimension(8),parameter :: tf_node    = [1,1,3,3,3,3,5,5]

    type :: eclipse_row
        !! an eclipse constraint: the KS aggregate of the shadow function at the samples
        !! `k_lo` to `k_hi` of a segment (sample `k` is `k*eclipse_opt_sample_step` seconds
        !! from the segment's start, along its direction of propagation)
        integer :: iseg = 0
        integer :: k_lo = 0
        integer :: k_hi = 0
        logical :: tf_with_t0 = .false. !! if the segment's final epoch moves with its start node's epoch
        integer :: ncol = 0             !! number of variables it depends on
        integer,dimension(7) :: col = 0  !! their indices in `x`
        integer,dimension(7) :: kind = 0 !! 0: the start node's epoch, 1-6: its state
    end type eclipse_row

    type :: eclipse_opt_data
        !! the data the sqpopt functions need
        type(mission_type),pointer :: mission => null()
        real(wp),dimension(:),allocatable :: x0  !! the solution before the optimization (scaled)
        real(wp),dimension(:),allocatable :: r   !! the optimization's variables are `x/r` (see `eclipse_opt_epoch_scale`)
        type(eclipse_row),dimension(:),allocatable :: rows
        integer :: m_halo = 0    !! number of halo's constraints
        integer,dimension(:),allocatable :: icol  !! the columns of halo's Jacobian nonzeros
        integer,dimension(:),allocatable :: irow  !! the rows of halo's Jacobian nonzeros
        real(wp),dimension(:),allocatable :: w    !! weights of halo's constraint rows in the optimization
                                                  !! (see `eclipse_opt_equilibrate`; 1 otherwise)
        integer :: nnz_halo = 0  !! number of nonzeros of halo's Jacobian
        logical :: single = .false. !! one eclipse row for all the segments: the KS aggregate of the
                                    !! segments' rows (see `eclipse_opt_single_row`)
        integer,dimension(:),allocatable :: ucol  !! with `single`: the row's columns (the union of the segments')
        integer,dimension(:),allocatable :: upos  !! with `single`: the index in `ucol` of each segment row's column
    end type eclipse_opt_data

    type,extends(nlesolver_type) :: eclipse_newton_type
        !! the Newton finish (see [[newton_finish]]): halo's continuity constraints, then
        !! the active eclipse rows as equalities
        type(eclipse_opt_data),pointer :: d => null()
        integer,dimension(:),allocatable :: active  !! the indices in `d%rows` of the rows that are equalities
        logical,dimension(:),allocatable :: mask    !! the same, as a mask of `d%rows`
    end type eclipse_newton_type

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Move a converged solution out of the Earth's shadow (see the module documentation).

    module subroutine eclipse_optimize(me,x)

    class(my_solver_type),intent(inout),target :: me
    real(wp),dimension(:),intent(inout) :: x !! the converged solution (scaled); replaced by the result

    type(eclipse_opt_data),target :: d
    type(sqpopt_problem_type) :: problem
    type(sqpopt_options_type) :: options
    type(sqpopt_qp_solver_type) :: qp_solver
    type(sqpopt_type) :: sqp
    type(sqpopt_results_type) :: res
    type(eclipse_row),dimension(:),allocatable :: new
    real(wp),dimension(:),allocatable :: lambda, x_lb, x_ub, c_lb, c_ub, ks, limits
    integer,dimension(:),allocatable :: irow, icol, irow_all, icol_all
    integer :: n, m, mt, iouter, istat, k, j, nnz, nr
    real(wp) :: phimin, shadow_min
    character(len=:),allocatable :: message
    integer(int64) :: count0, count1, rate

    associate (mission => me%mission)

    if (mission%jacobian_mode /= 2) error stop 'eclipse optimization requires jacobian_mode=2'

    call system_clock(count0, rate)
    d%mission => me%mission
    call mission%define_problem_size(n=n, m=m)
    d%x0 = x
    ! the optimization's own scaling of the variables: y = x/r
    allocate(d%r(n))
    d%r = one
    if (mission%eclipse_opt_epoch_scale > zero) then
        where (mission%x_is_epoch) d%r = mission%eclipse_opt_epoch_scale / mission%xscale
    end if
    d%m_halo = m
    d%single = mission%eclipse_opt_single_row
    call mission%get_sparsity_pattern(irow, icol)
    d%nnz_halo = size(irow)
    d%icol = icol
    d%irow = irow
    allocate(d%w(m))
    d%w = one
    if (mission%eclipse_opt_equilibrate) call equilibrate(d, x)
    allocate(d%rows(0))
    mission%segs%soft_fail = .true.  ! (a trial point can be far off: an integrator failure is reported to sqpopt)

    write(*,'(A)') ''
    write(*,'(A)') ' * Eclipse optimization'
    if (mission%sqpopt_solve) then
        ! TEMPORARY experiment: solve halo's problem (continuity only) with sqpopt
        allocate(new(0))
        write(*,'(A)') '   (sqpopt_solve: the continuity constraints only, from the initial guess)'
    else
        call find_rows(mission, x, new, phimin, shadow_min)
        write(*,'(A,F8.4,A,F8.1,A,I0,A)') '   converged solution: min phi = ', phimin, ', ', shadow_min, &
                                          ' min in the shadow, ', size(new), ' segments near the shadow'
    end if

    do iouter = 1, mission%eclipse_opt_max_outer

        call merge_rows(d%rows, new)
        nr = size(d%rows)
        if (nr==0 .and. .not. mission%sqpopt_solve) then
            write(*,'(A)') '   nothing near the shadow: no optimization needed'
            exit
        end if
        if (d%single) then
            mt = m + 1
            call union_columns(d)
            write(*,'(A,I0,A,I0,A,I0,A)') '   pass ', iouter, ': 1 eclipse constraint (', nr, ' segments, ', &
                                          sum(d%rows%k_hi - d%rows%k_lo + 1), ' samples)'
        else
            mt = m + nr
            write(*,'(A,I0,A,I0,A,I0,A)') '   pass ', iouter, ': ', nr, ' eclipse constraints (', &
                                          sum(d%rows%k_hi - d%rows%k_lo + 1), ' samples)'
        end if
        if (nr > 0) then
            allocate(ks(nr))
            call mission%put_x_in_segments(x)
            call eval_rows(d, ks)
            write(*,'(A,F8.4,A,F8.4)') '   smallest eclipse constraint (KS) at the start: ', minval(ks), &
                                       ', their aggregate: ', ks_aggregate(ks, mission%eclipse_opt_ks_rho)
            deallocate(ks)
        end if

        ! the problem: halo's continuity rows (equalities), then the eclipse rows
        ! (inequalities, KS >= phi_min)
        if (d%single) then
            nnz = d%nnz_halo + size(d%ucol)
        else
            nnz = d%nnz_halo + sum(d%rows%ncol)
        end if
        allocate(irow_all(nnz), icol_all(nnz))
        irow_all(1:d%nnz_halo) = irow
        icol_all(1:d%nnz_halo) = icol
        if (d%single) then
            irow_all(d%nnz_halo+1:) = m + 1
            icol_all(d%nnz_halo+1:) = d%ucol
        else
            j = d%nnz_halo
            do k = 1, nr
                irow_all(j+1:j+d%rows(k)%ncol) = m + k
                icol_all(j+1:j+d%rows(k)%ncol) = d%rows(k)%col(1:d%rows(k)%ncol)
                j = j + d%rows(k)%ncol
            end do
        end if
        allocate(x_lb(n), x_ub(n), c_lb(mt), c_ub(mt), lambda(mt))
        x_lb = -sqpopt_infinity
        x_ub =  sqpopt_infinity
        c_lb(1:m) = zero
        c_ub(1:m) = zero
        c_lb(m+1:) = mission%eclipse_opt_phi_min
        c_ub(m+1:) = sqpopt_infinity

        call problem%set_problem_size(n=n, m=mt)
        call problem%set_bounds(x_lb, x_ub, c_lb, c_ub)
        call problem%set_jacobian_sparsity(nnz, irow_all, icol_all)
        call problem%set_functions(fc=eclipse_fc, gjac=eclipse_gjac, data=d)
        if (mission%eclipse_opt_max_step_epoch > zero .or. mission%eclipse_opt_max_step_position > zero .or. &
            mission%eclipse_opt_max_step_velocity > zero) then
            limits = step_limits(d)
            call problem%set_max_step(limits)
        end if

        options%max_iter             = mission%eclipse_opt_max_iter
        options%print_level          = mission%eclipse_opt_print_level
        options%hessian_mode         = mission%eclipse_opt_hessian_mode
        options%direct_qp            = mission%eclipse_opt_direct_qp
        options%direct_least_squares = mission%eclipse_opt_direct_least_squares
        options%qp_solver_mode       = mission%eclipse_opt_qp_solver_mode
        options%linesearch_mode      = mission%eclipse_opt_linesearch_mode
        options%lbfgs_memory         = mission%eclipse_opt_lbfgs_memory
        options%scaling              = mission%eclipse_opt_scaling
        options%ctol                 = mission%eclipse_opt_ctol
        options%ktol                 = mission%eclipse_opt_ktol
        options%elastic_multiplier_limit = mission%eclipse_opt_elastic_multiplier_limit
        options%scaling_max_gradient = mission%eclipse_opt_scaling_max_gradient
        options%hessian_scale0       = mission%eclipse_opt_hessian_scale0
        options%inertia_control      = mission%eclipse_opt_inertia_control
        options%merit_mode           = mission%eclipse_opt_merit_mode
        options%restoration_mode     = mission%eclipse_opt_restoration_mode
        options%max_time             = mission%eclipse_opt_max_time
        options%linear_solver        = mission%eclipse_opt_linear_solver
        options%factorization_threads = mission%eclipse_opt_factorization_threads
        qp_solver%max_step           = mission%eclipse_opt_qp_max_step
        qp_solver%direct_max_changes = mission%eclipse_opt_direct_max_changes
        qp_solver%direct_tol         = mission%eclipse_opt_direct_tol
        qp_solver%sparse_qp%basis_rel_tol         = mission%eclipse_opt_basis_rel_tol
        qp_solver%sparse_qp%basis_repair          = mission%eclipse_opt_basis_repair
        qp_solver%sparse_qp%basis_eq_slack_weight = mission%eclipse_opt_basis_eq_slack_weight
        qp_solver%direct_relative_tol   = mission%eclipse_opt_direct_relative_tol
        qp_solver%direct_accept_inexact = mission%eclipse_opt_direct_accept_inexact

        call sqp%initialize(problem=problem, options=options, qp_solver=qp_solver)
        call sqp%solve(x/d%r, istat)
        call sqp%get_solution(x, lambda)
        x = x * d%r
        call sqp%get_results(res)
        write(*,'(A)') '   sqpopt: '//sqp%status_message()
        write(*,'(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3,A,F9.2,A)') '   sqpopt: ', res%iterations, ' iterations, ', &
            res%n_eval_fc, ' fc, ', res%n_eval_gjac, ' gjac, f = ', res%f, ', infeasibility = ', &
            res%feasibility_error, ', ', res%time, ' s'
        call sqp%destroy()
        deallocate(irow_all, icol_all, x_lb, x_ub, c_lb, c_ub, lambda)

        if (mission%sqpopt_solve) then
            ! TEMPORARY experiment: how well sqpopt alone converged halo's problem
            block
                real(wp),dimension(:),allocatable :: c
                allocate(c(m))
                call mission%constraint_violations(x, c)  ! (propagates the segments at x)
                write(*,'(A,ES12.4,A,ES12.4,A,ES10.2)') '   sqpopt_solve: halo''s constraint norm ', norm2(c), &
                    ', largest |c| ', maxval(abs(c)), ' (nlesolver_tol ', mission%nlesolver_tol
            end block
            exit
        end if

        ! converge the continuity, with the eclipse rows at their bound held there:
        if (mission%eclipse_opt_newton_finish) call newton_finish(d, x)

        ! re-converge the continuity (only):
        if (mission%eclipse_opt_polish) then
            write(*,'(A)') '   polish (Newton, continuity only):'
            call me%solve(x)
            call me%status(istat=istat, message=message)
            write(*,'(A)') '   polish: '//message
        end if

        ! check the whole trajectory:
        call find_rows(mission, x, new, phimin, shadow_min)
        write(*,'(A,I0,A,F8.4,A,F8.1,A)') '   after pass ', iouter, ': min phi = ', phimin, ', ', &
                                          shadow_min, ' min in the shadow'
        if (phimin >= zero) then
            write(*,'(A)') '   eclipse free'
            exit
        end if

    end do

    mission%segs%soft_fail = .false.
    call report_change(mission, d%x0, x)
    call system_clock(count1)
    write(*,'(A,F10.2,A)') '   eclipse optimization wall time: ', real(count1-count0,wp)/real(rate,wp), ' s'
    write(*,'(A)') ''

    end associate

    end subroutine eclipse_optimize
!*****************************************************************************************

!*****************************************************************************************
!>
!  The objective and the constraints for sqpopt: halo's continuity constraints, then the
!  eclipse rows.

    subroutine eclipse_fc(x, f, c, status, data)

    real(wp),dimension(:),intent(in)    :: x
    real(wp),intent(out)                :: f
    real(wp),dimension(:),intent(out)   :: c
    integer,intent(inout)               :: status
    class(*),optional,intent(inout)     :: data

    real(wp),dimension(:),allocatable :: ks
    logical :: ok

    if (.not. present(data)) error stop 'eclipse_fc: no data'
    select type (data)
    type is (eclipse_opt_data)
        ! (x here is the optimization's variable y = x/r)
        if (data%mission%eclipse_opt_objective == 0) then
            f = zero  ! (a feasibility problem: as SNOPT's feasible point mode)
        else
            f = 0.5_wp * sum((x - data%x0/data%r)**2)
        end if
        if (.not. epochs_in_range(data, x*data%r)) then
            c = zero
            status = 1  ! (outside the splined ephemeris: can't evaluate here)
            return
        end if
        call data%mission%constraint_violations(x*data%r, c(1:data%m_halo))  ! (also puts x in the segments)
        c(1:data%m_halo) = c(1:data%m_halo) * data%w
        if (data%single) then
            allocate(ks(size(data%rows)))
            call eval_rows(data, ks, ok=ok)
            c(data%m_halo+1) = ks_aggregate(ks, data%mission%eclipse_opt_ks_rho)
        else
            call eval_rows(data, c(data%m_halo+1:), ok=ok)
        end if
        if (.not. ok .or. any(data%mission%segs%failed)) status = 1  ! (can't evaluate here: sqpopt shortens the step)
    class default
        error stop 'eclipse_fc: invalid data'
    end select

    end subroutine eclipse_fc
!*****************************************************************************************

!*****************************************************************************************
!>
!  The gradient of the objective and the Jacobian values for sqpopt: halo's analytic
!  Jacobian, then the eclipse rows (central differences).

    subroutine eclipse_gjac(x, g, jac_val, accuracy, status, data)

    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(out)   :: g
    real(wp),dimension(:),intent(out)   :: jac_val
    integer,intent(in)                  :: accuracy
    integer,intent(inout)               :: status
    class(*),optional,intent(inout)     :: data

    real(wp),dimension(:),allocatable :: jac, ks, jr, w
    real(wp) :: rho, kmin
    integer :: k, j, o
    logical :: ok

    if (.not. present(data)) error stop 'eclipse_gjac: no data'
    select type (data)
    type is (eclipse_opt_data)
        ! (x here is the optimization's variable y = x/r: the columns are multiplied by r)
        if (data%mission%eclipse_opt_objective == 0) then
            g = zero
        else
            g = x - data%x0/data%r
        end if
        if (.not. epochs_in_range(data, x*data%r)) then
            jac_val = zero
            status = 1
            return
        end if
        call data%mission%compute_jacobian_stm(x*data%r, jac)  ! (also puts x in the segments)
        jac_val(1:data%nnz_halo) = jac * data%r(data%icol) * data%w(data%irow)
        allocate(ks(size(data%rows)), jr(sum(data%rows%ncol)))
        call eval_rows(data, ks, jr, ok=ok)
        if (.not. ok .or. any(data%mission%segs%failed)) status = 1
        call scale_row_columns(data, jr)
        if (data%single) then
            ! the aggregate's partials: the segment rows' partials, weighted by d(KS)/d(ks_k)
            rho  = data%mission%eclipse_opt_ks_rho
            kmin = minval(ks)
            w = exp(-rho*(ks - kmin))
            w = w / sum(w)
            jac_val(data%nnz_halo+1:) = zero
            o = 0
            do k = 1, size(data%rows)
                do j = 1, data%rows(k)%ncol
                    jac_val(data%nnz_halo + data%upos(o+j)) = jac_val(data%nnz_halo + data%upos(o+j)) + w(k)*jr(o+j)
                end do
                o = o + data%rows(k)%ncol
            end do
        else
            jac_val(data%nnz_halo+1:) = jr
        end if
    class default
        error stop 'eclipse_gjac: invalid data'
    end select

    end subroutine eclipse_gjac
!*****************************************************************************************

!*****************************************************************************************
!>
!  The eclipse rows (and, optionally, their partials, in the order of the rows' columns),
!  for the `x` already in the segments. The rows (segments) are done in parallel.

    subroutine eval_rows(d, ks, jac, ok, only)

    type(eclipse_opt_data),intent(inout) :: d
    real(wp),dimension(:),intent(out) :: ks
    real(wp),dimension(:),intent(out),optional :: jac
    logical,intent(out),optional :: ok  !! false if a row couldn't be evaluated (see [[row_ks]])
    logical,dimension(:),intent(in),optional :: only !! evaluate only these rows (the others' `ks` are 0,
                                                     !! and their partials are not set)

    integer :: k, j, offset
    integer,dimension(:),allocatable :: first_nz  !! index in `jac` of each row's first nonzero, minus 1
    real(wp) :: t0, tf, dt, h, fp, fm
    real(wp),dimension(6) :: x0r, xp
    real(wp) :: tp0, tpf
    logical :: ok_k, all_ok

    allocate(first_nz(size(d%rows)))
    offset = 0
    do k = 1, size(d%rows)
        first_nz(k) = offset
        offset = offset + d%rows(k)%ncol
    end do
    h = d%mission%eclipse_opt_fd_step
    all_ok = .true.

!$OMP PARALLEL DO PRIVATE(k,j,t0,tf,dt,x0r,xp,tp0,tpf,fp,fm,ok_k) SCHEDULE(dynamic)
    do k = 1, size(d%rows)
        if (present(only)) then
            if (.not. only(k)) then
                ks(k) = zero
                cycle
            end if
        end if
        associate (row => d%rows(k), seg => d%mission%segs(d%rows(k)%iseg))
            t0  = seg%data%t0
            tf  = seg%data%tf
            x0r = seg%data%x0_rotating
            ks(k) = row_ks(d%mission, seg, row, t0, tf, x0r, ok_k)
            if (.not. ok_k) all_ok = .false.
            if (present(jac)) then
                do j = 1, row%ncol
                    dt = h * d%mission%xscale(row%col(j))  ! the unscaled change [days, km or km/s]
                    call perturbed(row, j, t0, tf, x0r, +dt, tp0, tpf, xp)
                    fp = row_ks(d%mission, seg, row, tp0, tpf, xp, ok_k)
                    if (.not. ok_k) all_ok = .false.
                    call perturbed(row, j, t0, tf, x0r, -dt, tp0, tpf, xp)
                    fm = row_ks(d%mission, seg, row, tp0, tpf, xp, ok_k)
                    if (.not. ok_k) all_ok = .false.
                    jac(first_nz(k)+j) = (fp - fm) / (2.0_wp*h)
                end do
            end if
        end associate
    end do
!$OMP END PARALLEL DO
    if (present(ok)) ok = all_ok

    contains

        pure subroutine perturbed(row, j, t0, tf, x0r, delta, tp0, tpf, xp)
        !! the segment's inputs with variable `j` of the row changed by `delta`
        !! (everything is passed in: an internal procedure would see the host's
        !! variables, not the OpenMP thread's private copies)
        type(eclipse_row),intent(in) :: row
        integer,intent(in) :: j
        real(wp),intent(in) :: t0, tf, delta
        real(wp),dimension(6),intent(in) :: x0r
        real(wp),intent(out) :: tp0, tpf
        real(wp),dimension(6),intent(out) :: xp
        tp0 = t0; tpf = tf; xp = x0r
        select case (row%kind(j))
        case (0)
            tp0 = t0 + delta
            if (row%tf_with_t0) tpf = tf + delta
        case (1:6)
            xp(row%kind(j)) = x0r(row%kind(j)) + delta
        end select
        end subroutine perturbed

    end subroutine eval_rows
!*****************************************************************************************

!*****************************************************************************************
!>
!  The Newton finish (`eclipse_opt_newton_finish`): from the optimization's point, solve
!  halo's continuity constraints and, as equalities \( \mathrm{KS}_k = \phi_{min} \), the
!  eclipse rows that are at or near their bound (\( \mathrm{KS}_k < \phi_{min} \) plus
!  `eclipse_opt_newton_active_tol`), with halo's minimum-norm Newton solver (nlesolver, with
!  halo's sparse linear solver). The optimization's point is near feasible but its
!  continuity isn't converged; a solve of the continuity alone (`eclipse_opt_polish`) moves
!  the trajectory back into the shadow, but with the active rows held at their bound it
!  can't. A row left free that falls below its bound is added to the equalities, and the
!  solve repeated (from where it ended).

    subroutine newton_finish(d, x)

    type(eclipse_opt_data),intent(inout),target :: d
    real(wp),dimension(:),intent(inout) :: x  !! the optimization's point (halo's scaled variables); replaced by the solution

    type(eclipse_newton_type) :: nle
    real(wp),dimension(:),allocatable :: ks
    integer,dimension(:),allocatable :: irow, icol, irow_all, icol_all
    logical,dimension(:),allocatable :: is_active, below
    integer :: attempt, n, m, na, k, j, nnz, istat, mode, stage_iter
    real(wp) :: phi_min, damp, fnorm0, fnorm
    real(wp),dimension(:),allocatable :: fvec, x_save
    character(len=:),allocatable :: message
    integer(int64) :: count0, count1, rate

    associate (mission => d%mission)

    call system_clock(count0, rate)
    n = size(x)
    m = d%m_halo
    phi_min = mission%eclipse_opt_phi_min
    call mission%get_sparsity_pattern(irow, icol)
    mode = mission%solver_mode
    if (mission%eclipse_opt_newton_solver_mode > 0) mode = mission%eclipse_opt_newton_solver_mode
    if (all(mode /= [2,3,4,6,7])) mode = 7  ! (a sparse solver of nlesolver's own)
    allocate(ks(size(d%rows)), is_active(size(d%rows)), below(size(d%rows)))
    call mission%put_x_in_segments(x)
    call eval_rows(d, ks)
    is_active = ks < phi_min + mission%eclipse_opt_newton_active_tol
    block  ! TEMPORARY experiment: hold only the eclipse_opt_newton_max_active rows with the smallest KS
        integer :: kk, imin
        real(wp),dimension(:),allocatable :: kt
        if (mission%eclipse_opt_newton_max_active >= 0) then
            is_active = .false.
            kt = ks
            do kk = 1, min(mission%eclipse_opt_newton_max_active, size(ks))
                imin = minloc(kt, 1)
                is_active(imin) = .true.
                kt(imin) = huge(one)
            end do
        end if
    end block
    write(*,'(A,*(F10.6))') '   Newton finish: the rows'' KS at the start: ', ks

    do attempt = 1, 5

        nle%d => d
        nle%mask = is_active
        nle%active = pack([(k, k=1,size(d%rows))], is_active)
        na = size(nle%active)
        write(*,'(A,I0,A,I0,A)') '   Newton finish ', attempt, ': the continuity and ', na, &
                                 ' eclipse rows (KS = phi_min)'
        nnz = size(irow) + sum(d%rows(nle%active)%ncol)
        allocate(irow_all(nnz), icol_all(nnz))
        irow_all(1:size(irow)) = irow
        icol_all(1:size(irow)) = icol
        j = size(irow)
        do k = 1, na
            associate (row => d%rows(nle%active(k)))
                irow_all(j+1:j+row%ncol) = m + k
                icol_all(j+1:j+row%ncol) = row%col(1:row%ncol)
                j = j + row%ncol
            end associate
        end do

        ! the damping stages: from eclipse_opt_newton_damp down by 10 each stage to 0 (a damped
        ! stage has eclipse_opt_newton_damp_iter iterations; the undamped last one, up to
        ! eclipse_opt_newton_max_iter). A stage that ends worse than it started is undone.
        damp = mission%eclipse_opt_newton_damp
        allocate(fvec(m+na))
        do
            stage_iter = mission%eclipse_opt_newton_max_iter
            if (damp > zero) stage_iter = mission%eclipse_opt_newton_damp_iter
            nle%d => d
            nle%mask = is_active
            nle%active = pack([(k, k=1,size(d%rows))], is_active)
            call nle%initialize(n                = n,                      &
                                m                = m + na,                 &
                                max_iter         = stage_iter,             &
                                func             = newton_func,            &
                                grad_sparse      = newton_grad,            &
                                tol              = mission%nlesolver_tol,  &
                                tolx             = mission%nlesolver_tolx, &
                                step_mode        = 2,                      & ! backtracking (as halo's solve)
                                n_intervals      = 2,                      &
                                use_broyden      = .false.,                &
                                export_iteration = newton_export,          &
                                sparsity_mode    = mode,                   &
                                atol             = 1.0e-12_wp,             &
                                btol             = 1.0e-12_wp,             &
                                damp             = damp,                   &
                                itnlim           = 1000000,                &
                                irow             = irow_all,               &
                                icol             = icol_all,               &
                                lusol_method     = 0)
            call nle%status(istat=istat, message=message)
            if (istat /= 0) then
                write(*,'(A)') '   Newton finish: '//message
                error stop 'error initializing the Newton finish'
            end if
            call newton_func(nle, x, fvec)
            fnorm0 = norm2(fvec)
            x_save = x
            write(*,'(A,ES9.1,A,ES12.4)') '   Newton finish, damp ', damp, ': start norm ', fnorm0
            call nle%solve(x)
            call nle%status(istat=istat, message=message)
            call newton_func(nle, x, fvec)
            fnorm = norm2(fvec)
            write(*,'(A,ES12.4)') '   Newton finish: '//message//', norm ', fnorm
            if (.not. fnorm <= fnorm0) then
                write(*,'(A)') '   Newton finish: worse; the stage is undone'
                x = x_save
                fnorm = fnorm0
            end if
            if (damp <= zero) exit
            if (fnorm <= mission%nlesolver_tol) exit  ! (converged with damping)
            damp = damp / 10.0_wp
            if (damp < mission%eclipse_opt_newton_damp_min) damp = zero
        end do
        call nle%destroy()
        deallocate(irow_all, icol_all, fvec)

        ! the rows left free must still be at or above their bound:
        call mission%put_x_in_segments(x)
        call eval_rows(d, ks)
        below = (.not. is_active) .and. ks < phi_min - 1.0e-6_wp
        write(*,'(A,F10.6,A,F10.6,A,I0,A)') '   Newton finish: smallest KS ', minval(ks), ' (of the free rows: ', &
            minval(ks, mask=.not. is_active), '), ', count(below), ' free rows below phi_min'
        if (.not. any(below)) exit
        is_active = is_active .or. below

    end do
    call system_clock(count1)
    write(*,'(A,F10.2,A)') '   Newton finish wall time: ', real(count1-count0,wp)/real(rate,wp), ' s'

    end associate

    end subroutine newton_finish
!*****************************************************************************************

!*****************************************************************************************
!>
!  The Newton finish's functions: halo's continuity constraints, then `KS - phi_min` of
!  the active eclipse rows. A point that can't be evaluated (outside the splined
!  ephemeris, or an integrator failure) gets large values, so the line search backs off.

    subroutine newton_func(me, x, f)

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x  !! halo's scaled variables
    real(wp),dimension(:),intent(out)   :: f  !! the functions

    real(wp),dimension(:),allocatable :: ks
    integer :: m
    logical :: ok

    real(wp),parameter :: big = 1.0e10_wp

    select type (me)
    type is (eclipse_newton_type)
        m = me%d%m_halo
        if (.not. epochs_in_range(me%d, x)) then
            f = big
            return
        end if
        call me%d%mission%constraint_violations(x, f(1:m))  ! (also puts x in the segments)
        allocate(ks(size(me%d%rows)))
        call eval_rows(me%d, ks, ok=ok, only=me%mask)
        f(m+1:) = ks(me%active) - me%d%mission%eclipse_opt_phi_min
        if (.not. ok .or. any(me%d%mission%segs%failed)) f = big
    class default
        error stop 'newton_func: invalid class'
    end select

    end subroutine newton_func
!*****************************************************************************************

!*****************************************************************************************
!>
!  The Newton finish's Jacobian: halo's analytic one, then the active eclipse rows
!  (central differences, see [[eval_rows]]).

    subroutine newton_grad(me, x, g)

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x  !! halo's scaled variables
    real(wp),dimension(:),intent(out)   :: g  !! the Jacobian's nonzeros

    real(wp),dimension(:),allocatable :: jac, ks, jr
    integer,dimension(:),allocatable :: first_nz
    integer :: k, j, o, i
    logical :: ok

    select type (me)
    type is (eclipse_newton_type)
        call me%d%mission%compute_jacobian_stm(x, jac)  ! (also puts x in the segments)
        g(1:me%d%nnz_halo) = jac
        allocate(ks(size(me%d%rows)), jr(sum(me%d%rows%ncol)), first_nz(size(me%d%rows)))
        o = 0
        do k = 1, size(me%d%rows)
            first_nz(k) = o
            o = o + me%d%rows(k)%ncol
        end do
        call eval_rows(me%d, ks, jr, ok=ok, only=me%mask)
        j = me%d%nnz_halo
        do k = 1, size(me%active)
            i = me%active(k)
            g(j+1:j+me%d%rows(i)%ncol) = jr(first_nz(i)+1:first_nz(i)+me%d%rows(i)%ncol)
            j = j + me%d%rows(i)%ncol
        end do
    class default
        error stop 'newton_grad: invalid class'
    end select

    end subroutine newton_grad
!*****************************************************************************************

!*****************************************************************************************
!>
!  The Newton finish's iteration printout.

    subroutine newton_export(me, x, f, iter)

    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x     !! the variables
    real(wp),dimension(:),intent(in)    :: f     !! the functions
    integer,intent(in)                  :: iter  !! the iteration

    integer :: m

    select type (me)
    type is (eclipse_newton_type)
        m = me%d%m_halo
        if (size(f) > m) then
            write(*,'(A,I4,A,ES12.4,A,ES12.4)') '     iter ', iter, ': continuity norm ', norm2(f(1:m)), &
                                              ', eclipse rows max |KS - phi_min| ', maxval(abs(f(m+1:)))
        else
            write(*,'(A,I4,A,ES12.4)') '     iter ', iter, ': continuity norm ', norm2(f(1:m))
        end if
    class default
        error stop 'newton_export: invalid class'
    end select

    end subroutine newton_export
!*****************************************************************************************

!*****************************************************************************************
!>
!  Equilibrate the continuity Jacobian at `x` for the optimization (`eclipse_opt_equilibrate`):
!  alternately normalize its rows (weights `w`) and its columns (the variables' factors `r`,
!  the optimization's variables being `x/r`), a few times, so that the rows and the columns
!  all have 2-norms near 1. This lowers the condition number of the KKT systems of sqpopt's
!  QP subproblems (for 100 revolutions, J's condition number goes from 8e4 to 7e2). The
!  objective, `|x/r - x0/r|^2/2`, is then measured in these units.

    subroutine equilibrate(d, x)

    type(eclipse_opt_data),intent(inout) :: d
    real(wp),dimension(:),intent(in) :: x

    real(wp),dimension(:),allocatable :: jac, a, norms
    integer :: it, k

    call d%mission%compute_jacobian_stm(x, jac)
    allocate(a(size(jac)), norms(max(size(d%r), size(d%w))))
    d%r = one
    d%w = one
    do it = 1, 10
        a = jac * d%r(d%icol) * d%w(d%irow)
        norms = zero
        do k = 1, size(a)
            norms(d%irow(k)) = norms(d%irow(k)) + a(k)**2
        end do
        where (norms(1:size(d%w)) > zero) d%w = d%w / sqrt(norms(1:size(d%w)))
        a = jac * d%r(d%icol) * d%w(d%irow)
        norms = zero
        do k = 1, size(a)
            norms(d%icol(k)) = norms(d%icol(k)) + a(k)**2
        end do
        where (norms(1:size(d%r)) > zero) d%r = d%r / sqrt(norms(1:size(d%r)))
    end do
    write(*,'(A,2ES10.2,A,2ES10.2)') '   equilibrated: column factors ', minval(d%r), maxval(d%r), &
                                     ', row weights ', minval(d%w), maxval(d%w)

    end subroutine equilibrate
!*****************************************************************************************

!*****************************************************************************************
!>
!  The largest change of each of the optimization's variables in one sqpopt iteration
!  (`problem%set_max_step`), from the limits in physical units (`eclipse_opt_max_step_*`):
!  a variable's unscaled value is `y*r*xscale` (plus its offset), so its limit is
!  `limit/(r*xscale)`. A variable without a limit gets `sqpopt_infinity`.

    function step_limits(d) result(lim)

    type(eclipse_opt_data),intent(in) :: d
    real(wp),dimension(:),allocatable :: lim  !! (allocatable: sized by the problem, so not on the stack)

    integer :: i, k, n_full, cc
    integer,dimension(:),allocatable :: removed
    real(wp) :: phys

    call d%mission%define_problem_size(n=n_full, full_problem=.true.)
    removed = d%mission%removed_columns()
    allocate(lim(size(d%r)))
    lim = sqpopt_infinity
    k = 0
    do i = 1, n_full
        if (any(removed==i)) cycle
        k = k + 1
        cc = mod(i-1,7) + 1
        select case (cc)
        case (1);   phys = d%mission%eclipse_opt_max_step_epoch
        case (2:4); phys = d%mission%eclipse_opt_max_step_position
        case default; phys = d%mission%eclipse_opt_max_step_velocity
        end select
        if (phys > zero) lim(k) = phys / (d%r(k)*d%mission%xscale(k))
    end do

    end function step_limits
!*****************************************************************************************

!*****************************************************************************************
!>
!  If all the node epochs of `x` (halo's scaled variables) are within the range of the
!  splined ephemeris (2 periods before the reference epoch to 2 periods after the last rev,
!  see [[initialize_the_mission]]), less a margin of half a period. A trial point of the
!  optimization can be far off.

    function epochs_in_range(d, x) result(ok)

    type(eclipse_opt_data),intent(in) :: d
    real(wp),dimension(:),intent(in) :: x
    logical :: ok

    real(wp) :: lo, hi, t
    integer :: i

    lo = -1.5_wp * d%mission%period
    hi = (d%mission%n_revs + 1.5_wp) * d%mission%period
    ok = .true.
    do i = 1, size(x)
        if (.not. d%mission%x_is_epoch(i)) cycle
        t = x(i)*d%mission%xscale(i) + d%mission%xoffset(i)  ! [days from the reference epoch]
        if (t < lo .or. t > hi .or. t /= t) then
            ok = .false.
            return
        end if
    end do

    end function epochs_in_range
!*****************************************************************************************

!*****************************************************************************************
!>
!  The KS aggregate of values `ks`: \( k_{lo} - \frac{1}{\rho} \ln \sum_k e^{-\rho (k_k - k_{lo})} \).
!  Of the segments' KS rows, it is exactly the KS of all their samples.

    pure function ks_aggregate(ks, rho) result(agg)

    real(wp),dimension(:),intent(in) :: ks
    real(wp),intent(in) :: rho
    real(wp) :: agg

    real(wp) :: kmin

    kmin = minval(ks)
    agg = kmin - log(sum(exp(-rho*(ks - kmin)))) / rho

    end function ks_aggregate
!*****************************************************************************************

!*****************************************************************************************
!>
!  With `single`: the columns of the single eclipse row (the union of the segment rows'),
!  and where each segment row's column is in them.

    subroutine union_columns(d)

    type(eclipse_opt_data),intent(inout) :: d

    integer :: k, j, o, i, nu

    if (allocated(d%ucol)) deallocate(d%ucol)
    if (allocated(d%upos)) deallocate(d%upos)
    allocate(d%ucol(sum(d%rows%ncol)), d%upos(sum(d%rows%ncol)))
    nu = 0
    o = 0
    do k = 1, size(d%rows)
        do j = 1, d%rows(k)%ncol
            i = findloc(d%ucol(1:nu), d%rows(k)%col(j), 1)
            if (i == 0) then
                nu = nu + 1
                d%ucol(nu) = d%rows(k)%col(j)
                i = nu
            end if
            d%upos(o+j) = i
        end do
        o = o + d%rows(k)%ncol
    end do
    d%ucol = d%ucol(1:nu)

    end subroutine union_columns
!*****************************************************************************************

!*****************************************************************************************
!>
!  Multiply the eclipse rows' partials by the optimization's variable scaling `r`.

    subroutine scale_row_columns(d, jac)

    type(eclipse_opt_data),intent(in) :: d
    real(wp),dimension(:),intent(inout) :: jac

    integer :: k, j, o

    o = 0
    do k = 1, size(d%rows)
        do j = 1, d%rows(k)%ncol
            jac(o+j) = jac(o+j) * d%r(d%rows(k)%col(j))
        end do
        o = o + d%rows(k)%ncol
    end do

    end subroutine scale_row_columns
!*****************************************************************************************

!*****************************************************************************************
!>
!  The KS aggregate of the shadow function at a row's samples, for the given start epoch,
!  final epoch [days] and initial rotating-frame state (a copy of the segment is propagated,
!  with output every `eclipse_opt_sample_step` seconds).

    function row_ks(mission, seg, row, t0, tf, x0_rotating, ok) result(ks)

    class(mission_type),intent(in) :: mission
    type(segment),intent(in) :: seg
    type(eclipse_row),intent(in) :: row
    real(wp),intent(in) :: t0, tf
    real(wp),dimension(6),intent(in) :: x0_rotating
    logical,intent(out) :: ok  !! false if the segment has none of the row's samples (e.g., a trial point
                               !! that moved its epochs far): the row can't be evaluated there
    real(wp) :: ks

    type(segment) :: s
    real(wp),dimension(:),allocatable :: phi
    real(wp),dimension(6) :: rv
    real(wp) :: step, rho, pmin
    integer :: i, k, n, np

    step = mission%eclipse_opt_sample_step
    rho  = mission%eclipse_opt_ks_rho
    s = seg
    call s%traj_inertial%destroy()
    call s%traj_rotating%destroy()
    call s%traj_se_rotating%destroy()
    call s%set_input(t0, tf, x0_rotating)
    call s%propagate(mode=2, tstep=step)
    if (s%failed) then
        ok = .false.
        ks = zero
        return
    end if

    ! the samples: the points at k*step from the start (the last point, at tf, is only one
    ! of them if it falls on the grid)
    n = size(s%traj_inertial%et)
    allocate(phi(row%k_hi - row%k_lo + 1))
    np = 0
    do k = row%k_lo, row%k_hi
        i = k + 1
        if (i > n) exit
        if (abs(abs(s%traj_inertial%et(i) - s%traj_inertial%et(1)) - k*step) > 1.0e-3_wp) exit
        rv = [s%traj_inertial%x(i),  s%traj_inertial%y(i),  s%traj_inertial%z(i), &
              s%traj_inertial%vx(i), s%traj_inertial%vy(i), s%traj_inertial%vz(i)]
        np = np + 1
        phi(np) = get_sun_fraction(s%eph, s%traj_inertial%et(i), rv, mission%eclipse_opt_bubble)
    end do
    ok = np > 0
    if (.not. ok) then
        ks = zero
        return
    end if
    pmin = minval(phi(1:np))
    ks = pmin - log(sum(exp(-rho*(phi(1:np) - pmin)))) / rho

    end function row_ks
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the shadow function along the whole trajectory (every `eclipse_opt_detect_step`
!  seconds), and make the rows of the segments that come near the shadow (where the shadow
!  function, with the constraints' bubble, is below `eclipse_opt_detect_phi`).

    subroutine find_rows(mission, x, rows, phimin, shadow_min)

    class(mission_type),intent(inout) :: mission
    real(wp),dimension(:),intent(in) :: x
    type(eclipse_row),dimension(:),allocatable,intent(out) :: rows
    real(wp),intent(out) :: phimin      !! the smallest shadow function (no bubble) along the trajectory
    real(wp),intent(out) :: shadow_min  !! the time in the shadow [min] (approximate: the number of points
                                        !! in the shadow times the step)

    integer :: n_segs, iseg, i, n, n_full, nk, kmax, nr
    real(wp),dimension(:),allocatable :: seg_phimin, seg_lo, seg_hi
    integer,dimension(:),allocatable :: seg_nshadow, full2red, removed
    logical,dimension(:),allocatable :: seg_near
    real(wp) :: dt, phi0, phib, step
    real(wp),dimension(6) :: rv

    call mission%define_problem_size(n_segs=n_segs)
    call mission%define_problem_size(n=n_full, full_problem=.true.)
    call mission%put_x_in_segments(x)

    allocate(seg_phimin(n_segs), seg_lo(n_segs), seg_hi(n_segs), seg_nshadow(n_segs), seg_near(n_segs))
    seg_phimin = huge(one); seg_lo = huge(one); seg_hi = -huge(one); seg_nshadow = 0; seg_near = .false.

!$OMP PARALLEL DO PRIVATE(iseg,i,n,dt,phi0,phib,rv) SCHEDULE(dynamic)
    do iseg = 1, n_segs
        associate (seg => mission%segs(iseg))
            call seg%traj_inertial%destroy()
            call seg%traj_rotating%destroy()
            call seg%traj_se_rotating%destroy()
            call seg%propagate(mode=2, tstep=mission%eclipse_opt_detect_step)
            n = size(seg%traj_inertial%et)
            do i = 1, n
                rv = [seg%traj_inertial%x(i),  seg%traj_inertial%y(i),  seg%traj_inertial%z(i), &
                      seg%traj_inertial%vx(i), seg%traj_inertial%vy(i), seg%traj_inertial%vz(i)]
                phi0 = get_sun_fraction(seg%eph, seg%traj_inertial%et(i), rv, zero)
                if (mission%eclipse_opt_bubble > zero) then
                    phib = get_sun_fraction(seg%eph, seg%traj_inertial%et(i), rv, mission%eclipse_opt_bubble)
                else
                    phib = phi0
                end if
                seg_phimin(iseg) = min(seg_phimin(iseg), phi0)
                if (phi0 < zero) seg_nshadow(iseg) = seg_nshadow(iseg) + 1
                if (phib < mission%eclipse_opt_detect_phi) then
                    dt = abs(seg%traj_inertial%et(i) - seg%traj_inertial%et(1))  ! [sec from the start]
                    seg_near(iseg) = .true.
                    seg_lo(iseg) = min(seg_lo(iseg), dt)
                    seg_hi(iseg) = max(seg_hi(iseg), dt)
                end if
            end do
            call seg%traj_inertial%destroy()
            call seg%traj_rotating%destroy()
            call seg%traj_se_rotating%destroy()
        end associate
    end do
!$OMP END PARALLEL DO

    phimin = minval(seg_phimin)
    shadow_min = real(sum(seg_nshadow),wp) * mission%eclipse_opt_detect_step / 60.0_wp

    ! the map from the full problem's columns to x:
    removed = mission%removed_columns()
    allocate(full2red(n_full))
    full2red = 0
    nk = 0
    do i = 1, n_full
        if (any(removed==i)) cycle
        nk = nk + 1
        full2red(i) = nk
    end do

    ! the rows: the samples across the near-shadow part of each flagged segment, padded by one
    ! sample on each side (and not past the segment's last full step)
    step = mission%eclipse_opt_sample_step
    allocate(rows(count(seg_near)))
    nr = 0
    do iseg = 1, n_segs
        if (.not. seg_near(iseg)) cycle
        associate (seg => mission%segs(iseg))
            kmax = floor(abs(seg%data%tf - seg%data%t0)*day2sec/step - 1.0e-6_wp)
            nr = nr + 1
            rows(nr)%iseg = iseg
            rows(nr)%k_lo = max(0, floor(seg_lo(iseg)/step) - 1)
            rows(nr)%k_hi = min(kmax, ceiling(seg_hi(iseg)/step) + 1)
            call row_columns(iseg, full2red, rows(nr))
        end associate
    end do

    end subroutine find_rows
!*****************************************************************************************

!*****************************************************************************************
!>
!  The variables a row of segment `iseg` depends on: the epoch and state of its start node.

    subroutine row_columns(iseg, full2red, row)

    integer,intent(in) :: iseg
    integer,dimension(:),intent(in) :: full2red
    type(eclipse_row),intent(inout) :: row

    integer :: irev, l, js, ks, k, fc

    irev = (iseg-1)/8
    l    = iseg - 8*irev
    js   = 4*irev + start_node(l)
    ks   = 4*irev + tf_node(l)
    row%tf_with_t0 = ks == js
    row%ncol = 0
    do k = 1, 7
        fc = 7*(js-1) + k
        if (full2red(fc) > 0) then
            row%ncol = row%ncol + 1
            row%col(row%ncol)  = full2red(fc)
            row%kind(row%ncol) = k - 1
        end if
    end do

    end subroutine row_columns
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add the new rows: a segment already in the list gets the union of the sample ranges.

    subroutine merge_rows(rows, new)

    type(eclipse_row),dimension(:),allocatable,intent(inout) :: rows
    type(eclipse_row),dimension(:),intent(in) :: new

    integer :: k, j
    logical :: found

    do k = 1, size(new)
        found = .false.
        do j = 1, size(rows)
            if (rows(j)%iseg == new(k)%iseg) then
                rows(j)%k_lo = min(rows(j)%k_lo, new(k)%k_lo)
                rows(j)%k_hi = max(rows(j)%k_hi, new(k)%k_hi)
                found = .true.
                exit
            end if
        end do
        if (.not. found) rows = [rows, new(k)]
    end do

    end subroutine merge_rows
!*****************************************************************************************

!*****************************************************************************************
!>
!  Print how far the optimization moved the solution.

    subroutine report_change(mission, x0, x)

    class(mission_type),intent(inout) :: mission
    real(wp),dimension(:),intent(in) :: x0, x

    integer :: i, n_full, k, cc
    integer,dimension(:),allocatable :: removed
    real(wp) :: dt_max, dr_max, dv_max, d

    call mission%define_problem_size(n=n_full, full_problem=.true.)
    removed = mission%removed_columns()
    dt_max = zero; dr_max = zero; dv_max = zero
    k = 0
    do i = 1, n_full
        if (any(removed==i)) cycle
        k = k + 1
        d = abs(x(k) - x0(k)) * mission%xscale(k)  ! unscaled
        cc = mod(i-1,7) + 1
        select case (cc)
        case (1);   dt_max = max(dt_max, d)
        case (2:4); dr_max = max(dr_max, d)
        case (5:7); dv_max = max(dv_max, d)
        end select
    end do
    write(*,'(A,F10.3,A,F12.3,A,F10.3,A)') '   largest change of a node: ', dt_max*24.0_wp, ' h in epoch, ', &
                                           dr_max, ' km in position, ', dv_max*1000.0_wp, ' m/s in velocity'

    end subroutine report_change
!*****************************************************************************************

    end submodule halo_eclipse_opt_submodule
!*****************************************************************************************
