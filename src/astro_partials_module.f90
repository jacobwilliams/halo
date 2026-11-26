!*****************************************************************************************
!>
!  Analytic partial derivatives for trajectory sensitivities (state transition
!  matrices): gravity gradients, and the time derivatives of the two-body rotating
!  frame transformation.
!
!  These are general astrodynamics routines (they do not depend on the rest of HALO),
!  intended to be moved into the Fortran Astrodynamics Toolkit eventually:
!
!  * [[pines_gravity]]: spherical harmonic gravity acceleration and its gradient
!    (the 3x3 matrix of second partials of the potential), with Pines' nonsingular
!    formulation. It uses the same unnormalized coefficient arrays as the toolkit's
!    `geopotential_model_pines` (`cnm`, `snm`), so it can be used with that class.
!  * [[pointmass_gravity_gradient]], [[third_body_gravity_gradient]].
!  * [[two_body_rotating_frame_matrices]]: the 6x6 inertial-to-rotating state
!    transformation of the toolkit's `two_body_rotating_frame` (and its inverse),
!    and their exact time derivatives.
!
!### Pines' formulation
!
!  With \( r = |\mathbf{r}| \) and the direction cosines
!  \( \hat{\mathbf{s}} = \mathbf{r}/r = (s,t,u) \), the potential is
!
!  \( U = \sum_{n=0}^{N} \rho_n(r) \sum_{m=0}^{\min(n,M)} A_{n,m}(u)
!         \left[ C_{n,m} r_m(s,t) + S_{n,m} i_m(s,t) \right] \)
!
!  where \( \rho_n = (\mu/r)(R/r)^n \), \( A_{n,m} = d^m P_n / du^m \) are the derived
!  Legendre functions, and \( r_m + i\, i_m = (s + i t)^m \). Every term is a
!  polynomial in \( s,t,u \), so there are no singularities (e.g., at the poles).
!
!  Treating \( U \) as a function of \( r \) and of \( \hat{\mathbf{s}} \) (as
!  independent variables), with \( \mathbf{g} = \partial U/\partial \hat{\mathbf{s}} \),
!  \( G = \partial^2 U/\partial \hat{\mathbf{s}}^2 \) and the projection
!  \( P = I - \hat{\mathbf{s}}\hat{\mathbf{s}}^T \)
!  (since \( \partial \hat{\mathbf{s}}/\partial \mathbf{r} = P/r \)):
!
!  \( \mathbf{a} = U_r \hat{\mathbf{s}} + \frac{1}{r} P \mathbf{g} \)
!
!  \( \frac{\partial \mathbf{a}}{\partial \mathbf{r}} =
!     U_{rr} \hat{\mathbf{s}}\hat{\mathbf{s}}^T
!   + \frac{1}{r}\left( P \mathbf{g}_r \hat{\mathbf{s}}^T + \hat{\mathbf{s}} \mathbf{g}_r^T P \right)
!   + \frac{U_r}{r} P
!   + \frac{1}{r^2}\left[ P G P - \left( P \mathbf{g} \hat{\mathbf{s}}^T + \hat{\mathbf{s}} \mathbf{g}^T P \right)
!                         - (\hat{\mathbf{s}}\cdot\mathbf{g}) P \right] \)
!
!### References
!  * S. Pines, "Uniform Representation of the Gravitational Potential and its
!    Derivatives", AIAA Journal, Vol. 11, No. 11, 1973, pp. 1508-1511.
!  * J. L. Spencer, "Pines' Nonsingular Gravitational Potential Derivation,
!    Description and Implementation", NASA CR-147478, 1976.

    module astro_partials_module

    use halo_kinds_module, only: wp

    implicit none

    private

    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one  = 1.0_wp
    real(wp),parameter :: two  = 2.0_wp
    real(wp),parameter :: three = 3.0_wp

    public :: pines_gravity
    public :: pointmass_gravity_gradient
    public :: third_body_gravity_gradient
    public :: two_body_rotating_frame_matrices

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Spherical harmonic gravity acceleration, and optionally its gradient
!  \( \partial \mathbf{a} / \partial \mathbf{r} \), with Pines' formulation
!  (see the module documentation). Includes the central (point mass) term.
!
!  The coefficients are unnormalized, as stored by the toolkit's
!  `geopotential_model_pines` (`cnm(n,m)`, `snm(n,m)`, for `n>=1`, `m>=0`;
!  degree 1 is ignored, as in the toolkit's `gravpot`).
!
!  All vectors are in the body-fixed frame of the coefficients.

    pure subroutine pines_gravity(r,nmax,mmax,re,mu,cnm,snm,acc,grad)

    implicit none

    real(wp),dimension(3),intent(in)             :: r     !! position vector [km]
    integer,intent(in)                           :: nmax  !! degree to use (>=0)
    integer,intent(in)                           :: mmax  !! order to use (0 to `nmax`)
    real(wp),intent(in)                          :: re    !! body reference radius [km]
    real(wp),intent(in)                          :: mu    !! body gravitational parameter [km^3/s^2]
    real(wp),dimension(:,0:),intent(in)          :: cnm   !! unnormalized C coefficients `cnm(n,m)`
    real(wp),dimension(:,0:),intent(in)          :: snm   !! unnormalized S coefficients `snm(n,m)`
    real(wp),dimension(3),intent(out)            :: acc   !! acceleration [km/s^2]
    real(wp),dimension(3,3),intent(out),optional :: grad  !! gravity gradient \( \partial \mathbf{a}/\partial \mathbf{r} \) [1/s^2]

    real(wp),dimension(0:nmax,0:nmax+2) :: a  !! derived Legendre functions A(n,m) = d^m P_n / du^m
    real(wp),dimension(0:nmax) :: rm, im  !! real and imaginary parts of (s + i t)^m
    real(wp),dimension(3) :: shat  !! direction cosines (s,t,u)
    real(wp),dimension(3) :: g     !! dU/dshat
    real(wp),dimension(3) :: gr    !! d2U/dshat dr
    real(wp),dimension(3,3) :: gss !! d2U/dshat2
    real(wp),dimension(3) :: bn_g     !! gradient of B_n w.r.t. shat
    real(wp),dimension(3,3) :: bn_h   !! Hessian of B_n w.r.t. shat
    real(wp),dimension(3,3) :: p   !! projection I - shat*shat^T
    real(wp),dimension(3) :: pg, pgr
    real(wp) :: rmag, rinv, s, t, u
    real(wp) :: ur, urr  !! dU/dr, d2U/dr2
    real(wp) :: rho, rho1, rho2  !! rho_n and its first and second derivatives w.r.t. r
    real(wp) :: bn  !! B_n = sum_m A(n,m) D(n,m)
    real(wp) :: c, sn  !! C(n,m), S(n,m)
    real(wp) :: d, ds, dt, dss, dst  !! D(n,m) = C r_m + S i_m, and its derivatives w.r.t. s and t
    real(wp) :: sg  !! shat . g
    integer :: n, m, mtop, i, j
    logical :: need_grad

    need_grad = present(grad)

    rmag = norm2(r)
    rinv = one / rmag
    shat = r * rinv
    s = shat(1)
    t = shat(2)
    u = shat(3)

    ! derived Legendre functions A(n,m) (zero for m > n):
    a = zero
    a(0,0) = one
    if (nmax >= 1) then
        a(1,0) = u
        a(1,1) = one
    end if
    do n = 2, nmax
        a(n,0) = (real(2*n-1,wp)*u*a(n-1,0) - real(n-1,wp)*a(n-2,0)) / real(n,wp)
        do m = 1, n
            ! from P'_n - P'_{n-2} = (2n-1) P_{n-1}, differentiated m-1 times:
            a(n,m) = real(2*n-1,wp)*a(n-1,m-1) + a(n-2,m)
        end do
    end do

    ! (s + i t)^m:
    rm(0) = one
    im(0) = zero
    do m = 1, nmax
        rm(m) = s*rm(m-1) - t*im(m-1)
        im(m) = s*im(m-1) + t*rm(m-1)
    end do

    ! accumulate U_r, U_rr, g, g_r, G_ss:
    ur  = zero
    urr = zero
    g   = zero
    gr  = zero
    gss = zero
    rho = mu * rinv  ! rho_0
    do n = 0, nmax

        if (n == 1) then  ! degree 1 is zero (center of mass at the origin)
            rho = rho * re * rinv
            cycle
        end if

        rho1 = -real(n+1,wp) * rho * rinv
        rho2 =  real((n+1)*(n+2),wp) * rho * rinv * rinv

        bn   = zero
        bn_g = zero
        bn_h = zero
        if (n == 0) then
            bn = one  ! the central term (C00 = 1)
        else
            mtop = min(n, mmax)
            do m = 0, mtop
                c  = cnm(n,m)
                sn = snm(n,m)
                d = c*rm(m) + sn*im(m)
                if (m >= 1) then
                    ds = real(m,wp) * (c*rm(m-1) + sn*im(m-1))
                    dt = real(m,wp) * (sn*rm(m-1) - c*im(m-1))
                else
                    ds = zero
                    dt = zero
                end if
                bn = bn + a(n,m)*d
                bn_g(1) = bn_g(1) + a(n,m)*ds
                bn_g(2) = bn_g(2) + a(n,m)*dt
                bn_g(3) = bn_g(3) + a(n,m+1)*d
                if (need_grad) then
                    if (m >= 2) then
                        dss = real(m*(m-1),wp) * (c*rm(m-2) + sn*im(m-2))
                        dst = real(m*(m-1),wp) * (sn*rm(m-2) - c*im(m-2))
                    else
                        dss = zero
                        dst = zero
                    end if
                    bn_h(1,1) = bn_h(1,1) + a(n,m)*dss
                    bn_h(1,2) = bn_h(1,2) + a(n,m)*dst
                    bn_h(2,2) = bn_h(2,2) - a(n,m)*dss  ! (D_tt = -D_ss)
                    bn_h(1,3) = bn_h(1,3) + a(n,m+1)*ds
                    bn_h(2,3) = bn_h(2,3) + a(n,m+1)*dt
                    bn_h(3,3) = bn_h(3,3) + a(n,m+2)*d
                end if
            end do
        end if

        ur  = ur  + rho1*bn
        urr = urr + rho2*bn
        g   = g   + rho*bn_g
        gr  = gr  + rho1*bn_g
        if (need_grad) gss = gss + rho*bn_h

        rho = rho * re * rinv  ! rho_{n+1}

    end do

    ! acceleration:
    sg  = dot_product(shat, g)
    pg  = g - sg*shat
    acc = ur*shat + rinv*pg

    if (need_grad) then
        gss(2,1) = gss(1,2)
        gss(3,1) = gss(1,3)
        gss(3,2) = gss(2,3)
        p = -outer(shat,shat)
        do i = 1, 3
            p(i,i) = p(i,i) + one
        end do
        pgr = gr - dot_product(shat,gr)*shat
        grad = urr*outer(shat,shat) &
             + rinv*(outer(pgr,shat) + outer(shat,pgr)) &
             + (ur*rinv)*p &
             + (rinv*rinv)*( matmul(p,matmul(gss,p)) - (outer(pg,shat) + outer(shat,pg)) - sg*p )
        ! symmetrize (it is symmetric analytically):
        do j = 1, 3
            do i = j+1, 3
                grad(i,j) = 0.5_wp*(grad(i,j) + grad(j,i))
                grad(j,i) = grad(i,j)
            end do
        end do
    end if

    end subroutine pines_gravity
!*****************************************************************************************

!*****************************************************************************************
!>
!  Gravity gradient of a point mass at the origin:
!  \( \partial \mathbf{a}/\partial \mathbf{r} = \frac{\mu}{r^3}(3\hat{\mathbf{r}}\hat{\mathbf{r}}^T - I) \).

    pure function pointmass_gravity_gradient(r,mu) result(grad)

    implicit none

    real(wp),dimension(3),intent(in) :: r   !! position vector [km]
    real(wp),intent(in)              :: mu  !! gravitational parameter [km^3/s^2]
    real(wp),dimension(3,3)          :: grad

    grad = point_gradient(r,mu)

    end function pointmass_gravity_gradient
!*****************************************************************************************

!*****************************************************************************************
!>
!  Gradient of the third-body acceleration
!  \( \mathbf{a} = \mu \left( \frac{\mathbf{d}}{|\mathbf{d}|^3} - \frac{\mathbf{r}_b}{|\mathbf{r}_b|^3} \right) \),
!  \( \mathbf{d} = \mathbf{r}_b - \mathbf{r} \), with respect to the spacecraft position
!  \( \mathbf{r} \) (the indirect term does not depend on it):
!  \( \frac{\mu}{|\mathbf{d}|^3}(3\hat{\mathbf{d}}\hat{\mathbf{d}}^T - I) \).
!  This is also the gradient of Battin's equivalent form (the same function).

    pure function third_body_gravity_gradient(r,rb,mu) result(grad)

    implicit none

    real(wp),dimension(3),intent(in) :: r   !! spacecraft position vector [km]
    real(wp),dimension(3),intent(in) :: rb  !! third-body position vector [km]
    real(wp),intent(in)              :: mu  !! third-body gravitational parameter [km^3/s^2]
    real(wp),dimension(3,3)          :: grad

    grad = point_gradient(r-rb,mu)  ! (a function of r - rb only)

    end function third_body_gravity_gradient
!*****************************************************************************************

!*****************************************************************************************
!>
!  \( \frac{\mu}{|\mathbf{x}|^3}(3\hat{\mathbf{x}}\hat{\mathbf{x}}^T - I) \)

    pure function point_gradient(x,mu) result(grad)

    real(wp),dimension(3),intent(in) :: x
    real(wp),intent(in)              :: mu
    real(wp),dimension(3,3)          :: grad

    real(wp) :: xmag2, xmag, c
    real(wp),dimension(3) :: xhat
    integer :: i

    xmag2 = dot_product(x,x)
    xmag  = sqrt(xmag2)
    xhat  = x / xmag
    c     = mu / (xmag2*xmag)
    grad  = (three*c) * outer(xhat,xhat)
    do i = 1, 3
        grad(i,i) = grad(i,i) - c
    end do

    end function point_gradient
!*****************************************************************************************

!*****************************************************************************************
!>
!  The 6x6 state transformations of a two-body rotating frame (as defined by the
!  toolkit's `two_body_rotating_frame`), and their exact time derivatives.
!
!  The frame's axes are \( \hat{\mathbf{x}} = \mathbf{r}/r \),
!  \( \hat{\mathbf{z}} = \mathbf{h}/h \) (\( \mathbf{h} = \mathbf{r}\times\mathbf{v} \)),
!  \( \hat{\mathbf{y}} = \hat{\mathbf{z}}\times\hat{\mathbf{x}} \), where
!  \( \mathbf{r},\mathbf{v} \) is the secondary body's state w.r.t. the primary body.
!  With \( C \) the matrix whose rows are the axes, and (as in the toolkit)
!  \( \dot{C}_f = -C[\boldsymbol{\omega}\times] \), \( \boldsymbol{\omega} = \mathbf{h}/r^2 \):
!
!  * inertial to rotating: \( M = \begin{bmatrix} C & 0 \\ \dot{C}_f & C \end{bmatrix} \),
!  * rotating to inertial: \( N = M^{-1} = \begin{bmatrix} C^T & 0 \\ \dot{C}_f^T & C^T \end{bmatrix} \).
!
!  (Both frames are assumed to have the same origin.) The time derivatives of
!  \( M \) and \( N \) are computed exactly from the time derivatives of
!  \( \mathbf{r} \) and \( \mathbf{v} \) (`rdot`, `vdot`). Note that \( \dot{C}_f \) is
!  the toolkit's definition, which is not exactly \( dC/dt \) (it neglects the
!  rotation of \( \mathbf{h} \)); the derivatives here are of the transformation as
!  it is actually defined.

    pure subroutine two_body_rotating_frame_matrices(rv,rdot,vdot,m,mdot,n,ndot)

    implicit none

    real(wp),dimension(6),intent(in)  :: rv    !! state of the secondary w.r.t. the primary [km, km/s]
    real(wp),dimension(3),intent(in)  :: rdot  !! time derivative of `rv(1:3)` [km/s]
    real(wp),dimension(3),intent(in)  :: vdot  !! time derivative of `rv(4:6)` [km/s^2]
    real(wp),dimension(6,6),intent(out),optional :: m     !! inertial to rotating
    real(wp),dimension(6,6),intent(out),optional :: mdot  !! dM/dt [1/s]
    real(wp),dimension(6,6),intent(out),optional :: n     !! rotating to inertial
    real(wp),dimension(6,6),intent(out),optional :: ndot  !! dN/dt [1/s]

    real(wp),dimension(3) :: r, v, h, hdot, xhat, yhat, zhat, xhatd, yhatd, zhatd, w, wdot
    real(wp),dimension(3,3) :: c, cd, cf, cfd
    real(wp) :: rmag, hmag

    r = rv(1:3)
    v = rv(4:6)
    rmag = norm2(r)
    h    = cross(r,v)
    hmag = norm2(h)
    hdot = cross(rdot,v) + cross(r,vdot)

    xhat = r / rmag
    zhat = h / hmag
    yhat = cross(zhat,xhat)
    xhatd = (rdot - xhat*dot_product(xhat,rdot)) / rmag
    zhatd = (hdot - zhat*dot_product(zhat,hdot)) / hmag
    yhatd = cross(zhatd,xhat) + cross(zhat,xhatd)

    c(1,:)  = xhat ; c(2,:)  = yhat ; c(3,:)  = zhat
    cd(1,:) = xhatd; cd(2,:) = yhatd; cd(3,:) = zhatd   ! dC/dt

    w    = h / (rmag*rmag)
    wdot = hdot/(rmag*rmag) - two*h*dot_product(r,rdot)/rmag**4
    cf   = -matmul(c, cross_matrix(w))                               ! toolkit's "cdot"
    cfd  = -matmul(cd, cross_matrix(w)) - matmul(c, cross_matrix(wdot))  ! d(cf)/dt

    if (present(m)) then
        m = zero
        m(1:3,1:3) = c
        m(4:6,1:3) = cf
        m(4:6,4:6) = c
    end if
    if (present(mdot)) then
        mdot = zero
        mdot(1:3,1:3) = cd
        mdot(4:6,1:3) = cfd
        mdot(4:6,4:6) = cd
    end if
    if (present(n)) then
        n = zero
        n(1:3,1:3) = transpose(c)
        n(4:6,1:3) = transpose(cf)
        n(4:6,4:6) = transpose(c)
    end if
    if (present(ndot)) then
        ndot = zero
        ndot(1:3,1:3) = transpose(cd)
        ndot(4:6,1:3) = transpose(cfd)
        ndot(4:6,4:6) = transpose(cd)
    end if

    end subroutine two_body_rotating_frame_matrices
!*****************************************************************************************

!*****************************************************************************************
!  small utilities

    pure function outer(a,b) result(c)
    real(wp),dimension(3),intent(in) :: a, b
    real(wp),dimension(3,3) :: c
    integer :: j
    do j = 1, 3
        c(:,j) = a*b(j)
    end do
    end function outer

    pure function cross(a,b) result(c)
    real(wp),dimension(3),intent(in) :: a, b
    real(wp),dimension(3) :: c
    c = [a(2)*b(3) - a(3)*b(2), &
         a(3)*b(1) - a(1)*b(3), &
         a(1)*b(2) - a(2)*b(1)]
    end function cross

    pure function cross_matrix(w) result(x)  !! [w x]: x*v = w cross v
    real(wp),dimension(3),intent(in) :: w
    real(wp),dimension(3,3) :: x
    x(1,:) = [ zero, -w(3),  w(2)]
    x(2,:) = [ w(3),  zero, -w(1)]
    x(3,:) = [-w(2),  w(1),  zero]
    end function cross_matrix

!*****************************************************************************************
    end module astro_partials_module
!*****************************************************************************************
