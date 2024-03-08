!*****************************************************************************************
!>
!  Various astrodynamics utilities.
!
!@todo Some of this should eventually be moved into the fortran astrodynamics toolkit

    module halo_utilities_module

    use fortran_astrodynamics_toolkit
    use iso_fortran_env
    use halo_kinds_module
    use splined_ephemeris_module
    use parameters_module

    implicit none

    private

    public :: from_j2000moon_to_j2000ssb
    public :: apparent_position
    public :: get_sun_fraction

    contains
!********************************************************************************

!********************************************************************************
    subroutine from_j2000moon_to_j2000ssb(eph, et, rv_moon, rv_ssb)

        !! convert from a j2000-moon frame to a j2000-ssb frame.

        class(jpl_ephemeris),intent(inout) :: eph
        real(wp),intent(in) :: et !! ephemeris time (sec)
        real(wp),dimension(6),intent(in) :: rv_moon !! j2000-moon state (km, km/s)
        real(wp),dimension(6),intent(out) :: rv_ssb !! j2000-ssb state (km, km/s)

        type(icrf_frame) :: f1, f2
        logical :: status_ok

        f1 = icrf_frame(b=body_moon)
        f2 = icrf_frame(b=body_ssb)
        call f1%transform(rv_moon,f2,et,eph,rv_ssb,status_ok) ! from f1 to f2
        if (.not. status_ok) error stop 'transformation error in from_j2000moon_to_j2000ssb'

    end subroutine from_j2000moon_to_j2000ssb
!********************************************************************************

!********************************************************************************
    subroutine get_pos(eph,et,b_target,b_obs,r,status_ok)

        !! ephemeris wrapper to just return position vector
        !! see also: [[ballistic_derivs]]

        class(jpl_ephemeris),intent(inout) :: eph
        real(wp),intent(in) :: et !! ephemeris time (sec)
        type(celestial_body),intent(in) :: b_target !! target body
        type(celestial_body),intent(in) :: b_obs !! observer body
        real(wp),dimension(3),intent(out) :: r !! j2000 state (km, km/s)
        logical,intent(out) :: status_ok !! true if no problems

        real(wp),dimension(6) :: rv !! r,v vector

        select type (eph)
        type is (jpl_ephemeris)
            call eph%get_rv(et,b_target,b_obs,rv,status_ok)
            r = rv(1:3)
        type is (jpl_ephemeris_splined)
            ! for this, we can just get position vector only
            call eph%get_r(et,b_target,b_obs,r,status_ok)
        class default
            error stop 'error in get_pos: invalid eph class'
        end select

    end subroutine get_pos
!********************************************************************************

!********************************************************************************
    subroutine apparent_position(eph, b_target, et, rv_obs_ssb, r_target, status_ok)

        !! Return the position of a target body relative to an observer,
        !! corrected for light time and stellar aberration.
        !!
        !! see the SPICELIB routine `spkapo` (with 'lt+s')

        class(jpl_ephemeris),intent(inout) :: eph !! the ephemeris
        type(celestial_body),intent(in) :: b_target !! target body
        real(wp),dimension(6),intent(in) :: rv_obs_ssb !! state of the observer
                                                       !! (j2000 frame w.r.t. solar system barycenter)
        real(wp),intent(in) :: et !! observer ephemeris time (sec)
        real(wp),dimension(3),intent(out) :: r_target !! apparant state of the target (j2000 frame)
                                                      !! Corrected for one-way light time and stellar aberration
        logical,intent(out) :: status_ok !! true if no problems

        real(wp),parameter :: c = 299792.458_wp !! speed of light in km/s

        real(wp),dimension(3) :: r_targ_ssb !! target body r wrt. ssb
        real(wp),dimension(6) :: rv_targ_ssb !! target body rv wrt. ssb
        real(wp) :: lt !! one-way light time [sec]

        ! Find the geometric position of the target body with respect to the
        ! solar system barycenter. Subtract the position of the observer
        ! to get the relative position. Use this to compute the one-way
        ! light time.
        call get_pos(eph,et,b_target,body_ssb,r_targ_ssb,status_ok)
        if (.not. status_ok) return
        r_targ_ssb = r_targ_ssb - rv_obs_ssb(1:3) ! relative pos of target
        lt = norm2(r_targ_ssb) / c ! light time

        ! To correct for light time, find the position of the target body
        ! at the current epoch minus the one-way light time. Note that
        ! the observer remains where he is.
        call get_pos(eph,et-lt,b_target,body_ssb,r_targ_ssb,status_ok)
        if (.not. status_ok) return
        r_targ_ssb = r_targ_ssb - rv_obs_ssb(1:3)

        ! At this point, r_targ_ssb contains the geometric or light-time
        ! corrected position of the target relative to the observer

        ! stellar aberration correction
        r_target = stellar_aberration(r_targ_ssb,rv_obs_ssb(4:6))

        contains

        function stellar_aberration ( pobj, vobs ) result(appobj)
            !!  Correct the apparent position of an object for stellar aberration.
            !!  see SPICELIB routine `STELAB`

            real(wp),dimension(3),intent(in) :: pobj
            real(wp),dimension(3),intent(in) :: vobs
            real(wp),dimension(3) :: appobj

            real(wp),dimension(3) :: u, vbyc,h
            real(wp) :: lensqr, sinphi, phi
            real(wp),parameter :: zero_tol = epsilon(1.0_wp) !! tolerance for zero

            u = unit(pobj)
            vbyc = vobs / c
            lensqr = dot_product ( vbyc, vbyc )
            if ( lensqr >= 1.0_wp) error stop 'velocity > speed of light'
            h = cross(u, vbyc)
            sinphi  = norm2 ( h )
            if ( abs(sinphi) > zero_tol ) then  ! if (sinphi /= 0)
                ! rotate the position of the object by phi
                ! radians about h to obtain the apparent position.
                phi = asin ( sinphi )
                call axis_angle_rotation ( pobj, h, phi, appobj )
            else
                ! observer is moving along the line of sight to the object,
                ! and no correction is required
                appobj = pobj
            end if

        end function stellar_aberration

    end subroutine apparent_position
!********************************************************************************

!********************************************************************************
!>
!  Compute the "sun fraction" using the circular cubic shadow model.

    function get_sun_fraction(eph, et, rv, rbubble) result (phi)

        class(jpl_ephemeris),intent(inout) :: eph !! the ephemeris
        real(wp),dimension(6),intent(in) :: rv !! state of the spacecraft !! (j2000-moon frame)
        real(wp),intent(in) :: et !! observer ephemeris time (sec)
        real(wp),intent(in) :: rbubble !! eclipse bubble [km]. see the reference.
        real(wp) :: phi !! circular cubic sun frac value.
                        !!
                        !!  * >0 no eclipse,
                        !!  * <0 eclipse,
                        !!  * =0 on the eclipse line

        logical :: status_ok !! true if no problems
        real(wp),dimension(3) :: r_sun !! apparent state of the sun (j2000-ssb frame)
        real(wp),dimension(3) :: r_earth !! apparent state of the earth (j2000-ssb frame)
        real(wp),dimension(6) :: rv_ssb !! state of the spacecraft !! (j2000-ssb frame)

        call from_j2000moon_to_j2000ssb(eph, et, rv, rv_ssb) ! state of spacecraft in j2000-ssb
        call apparent_position(eph, body_sun,   et, rv_ssb, r_sun,   status_ok) ! apparent position of sun in j2000
        call apparent_position(eph, body_earth, et, rv_ssb, r_earth, status_ok) ! apparent position of earth in j2000
        call cubic_shadow_model(r_sun, rad_sun, r_earth, rad_earth, phi, rbubble) ! compute sun fraction value

    end function get_sun_fraction
!********************************************************************************

!********************************************************************************
!>
!  The "circular cubic" shadow model.
!
!### Reference
!  * J. Williams, et. al, "A new eclipse algorithm for use in
!    spacecraft trajectory optimization", 2023, AAS 23-243

    subroutine cubic_shadow_model(rsun, radsun, rplanet, radplanet, sunfrac, rbubble)

    real(wp),dimension(3),intent(in)   :: rsun      !! apparent position vector of sun wrt spacecraft [km]
    real(wp), intent(in)               :: radsun    !! radius of sun [km]
    real(wp),dimension(3), intent(in)  :: rplanet   !! apparent position vector of eclipsing body wrt spacecraft [km]
    real(wp), intent(in)               :: radplanet !! radius of the eclipsing body [km]
    real(wp), intent(out)              :: sunfrac   !! value of the function (>0 no eclipse,
                                                    !! <0 eclipse, =0 on the shadow line)
    real(wp),intent(in),optional       :: rbubble   !! eclipse bubble radius. if present, then `sunfrac` is
                                                    !! the value along an arc length of `rbubble`
                                                    !! in the direction of the max eclipse line.

    real(wp),dimension(3) :: r   !! radius vector from eclipsing body to spacecraft
    real(wp),dimension(3) :: rsb !! radius vector from the sun to the eclipsing body
    real(wp) :: tmp              !! temp value
    real(wp) :: alpha            !! [deg]
    real(wp) :: alpha0           !! [deg]
    real(wp) :: sin_alpha0       !! `sin(alpha0)`
    real(wp) :: rsbmag           !! magnitude of radius vector from the sun to the eclipsing body
    real(wp) :: rmag             !! magnitude of `r`
    logical :: compute_bubble    !! use the `rbubble` inputs to adjust `alpha`

    compute_bubble = present(rbubble)
    if (compute_bubble) compute_bubble = rbubble > zero

    r      = -rplanet
    rmag   = norm2(r)
    if (rmag<radplanet) then
        ! if inside the body, just return value from the surface
        r    = radplanet * unit(r)
        rmag = radplanet
    end if
    rsb        = rplanet - rsun
    alpha      = safe_acosd(dot_product(unit(r),unit(rsb)))
    if (compute_bubble) alpha = rad2deg*abs(wrap_angle(alpha*deg2rad-abs(rbubble)/rmag))
    rsbmag     = norm2(rsb)
    tmp        = (radsun + radplanet) / rsbmag
    sin_alpha0 = (one/rmag)*(radplanet*sqrt((one/tmp)**2-one)+sqrt(rmag**2-radplanet**2))*tmp
    alpha0     = safe_asind(sin_alpha0)
    sunfrac    = (alpha**2/(alpha0**2-alpha0**3/270.0_wp))*(one-alpha/270.0_wp)-one

    contains

        pure real(wp) function safe_asind(x)
            !! `asind` with range checking
            real(wp),intent(in) :: x
            safe_asind = asind(min(one,max(-one,x)))
        end function safe_asind

        pure real(wp) function safe_acosd(x)
            !! `acosd` with range checking
            real(wp),intent(in) :: x
            safe_acosd = acosd(min(one,max(-one,x)))
        end function safe_acosd

    end subroutine cubic_shadow_model

end module halo_utilities_module
