!*****************************************************************************************
!>
!  Various astrodynamics utilities.
!
!  The lighting routines it used to have (`apparent_position`, `cubic_shadow_model`,
!  and the frame conversion to the solar system barycenter) are now in the Fortran
!  Astrodynamics Toolkit's `lighting_module`.

    module halo_utilities_module

    ! (the toolkit's get_sun_fraction is more general [any eclipsing body, with the state
    !  relative to it]; the one here is for halo: Moon-centered state, eclipses by the Earth)
    use fortran_astrodynamics_toolkit, fat_get_sun_fraction => get_sun_fraction
    use iso_fortran_env
    use halo_kinds_module
    use splined_ephemeris_module
    use parameters_module

    implicit none

    private

    public :: get_sun_fraction

    contains
!********************************************************************************

!********************************************************************************
!>
!  Compute the "sun fraction" for eclipses by the Earth, using the circular
!  cubic shadow model (the apparent positions of the Sun and the Earth).

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

        call from_j2000body_to_j2000ssb(body_moon, eph, et, rv, rv_ssb) ! state of spacecraft in j2000-ssb
        call apparent_position(eph, body_sun,   et, rv_ssb, r_sun,   status_ok) ! apparent position of sun in j2000
        call apparent_position(eph, body_earth, et, rv_ssb, r_earth, status_ok) ! apparent position of earth in j2000
        call cubic_shadow_model(r_sun, rad_sun, r_earth, rad_earth, phi, rbubble) ! compute sun fraction value

    end function get_sun_fraction
!********************************************************************************

end module halo_utilities_module
