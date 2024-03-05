    module parameters_module

    use halo_kinds_module,             only: wp
    use fortran_astrodynamics_toolkit, only: zero

    implicit none

    public

    public :: wp ! real kind

    ! constants:
    real(wp),parameter :: mu_earth = 3.9860043543609593E+05_wp !! \( \mu_{earth} \)
    real(wp),parameter :: mu_moon  = 4.9028000661637961E+03_wp !! \( \mu_{moon} \)
    real(wp),parameter :: mu_sun   = 1.3271244004193938E+11_wp !! \( \mu_{sun} \)
    real(wp),parameter :: r_moon   = 1737.4_wp      !! radius of hte Moon
    integer,parameter  :: n_eoms   = 6              !! size of EOM derivative vector [x,y,z,vx,vy,vz]
    real(wp),parameter :: rtol     = 1.0e-12_wp     !! integrator tols
    real(wp),parameter :: atol     = 1.0e-12_wp     !! integrator tols
    integer,parameter  :: maxnum   = 10000          !! integrator max steps
    integer,parameter  :: grav_n   = 8              !! max grav degree
    integer,parameter  :: grav_m   = 8              !! max grav order

    ! scale values (good values for the 4500 Rp case)
    real(wp),dimension(6),parameter,public :: xscale_x0 = &
        [1.0e+05_wp,1.0e+05_wp,1.0e+05_wp,2.0e+00_wp,2.0e+00_wp,2.0e+00_wp]
    real(wp),dimension(6),parameter,public :: fscale_xf = &
        [1.0e+04_wp,1.0e+04_wp,1.0e+04_wp,1.0e+02_wp,1.0e+02_wp,1.0e+02_wp]

    type,public :: patch_point
        !! a CR3BP Halo state, to be used as
        !! an initial guess for the full force model.
        real(wp) :: t = zero                !! Time from Periapse (days from periapsis)
        real(wp),dimension(6) :: rv = zero  !! Rx (km)
                                            !! Ry (km)
                                            !! Rz (km)
                                            !! Vx (km/s)
                                            !! Vy (km/s)
                                            !! Vz (km/s) [moon-centered earth-moon rotating frame]
    end type patch_point

    end module parameters_module
