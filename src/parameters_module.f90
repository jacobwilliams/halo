    module parameters_module

    use halo_kinds_module,             only: wp
    use fortran_astrodynamics_toolkit, only: zero

    implicit none

    private

    public :: wp ! real kind

    ! constants:
    integer,parameter,public  :: n_eoms = 6  !! size of EOM derivative vector [x,y,z,vx,vy,vz]

    real(wp),public :: mu_earth   = 3.9860043543609593E+05_wp !! \( \mu_{earth} \)
    real(wp),public :: mu_moon    = 4.9028000661637961E+03_wp !! \( \mu_{moon} \)
    real(wp),public :: mu_sun     = 1.3271244004193938E+11_wp !! \( \mu_{sun} \)
    real(wp),public :: mu_jupiter = 1.266865349218008E+08_wp !! \( \mu_{jupiter} \) -- from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
    real(wp),public :: rad_moon   = 1737.4_wp      !! radius of the Moon
    real(wp),public :: rad_sun    = 695700.0_wp    !! sun radius (km)
    real(wp),public :: rad_earth  = 6378.0_wp      !! earth radius (km)
    integer ,public :: maxnum     = 10000          !! integrator max steps
    integer ,public :: grav_n     = 8              !! max grav degree
    integer ,public :: grav_m     = 8              !! max grav order
    real(wp),dimension(:),allocatable,public :: xscale_x0 ! scale values for opt vars [must be size 6]
    real(wp),dimension(:),allocatable,public :: fscale_xf ! scale values for constraints [must be size 6]

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
