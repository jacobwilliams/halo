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
    integer, public :: grav_frame = 1 !! 1: iau_moon, 2: moon_pa (splined)
    real(wp),dimension(:),allocatable,public :: xscale_x0 ! scale values for opt vars [must be size 6]
    real(wp),dimension(:),allocatable,public :: fscale_xf ! scale values for constraints [must be size 6]
    real(wp),public :: fscale_rdot = 1.0_wp  !! scale value for the rdot=0 constraint
    logical,public :: use_battin_gravity = .false. !! use Battin gravity formulation
    character(len=:),allocatable,public :: mkspk_path  !! for mkspk: path to mkspk executable. e.g. `kernel/mkspk`
    integer,public :: object_id = -50000                            !! for mkspk: OBJECT_ID
    character(len=:),allocatable,public :: object_name       !! for mkspk: OBJECT_NAME
    character(len=:),allocatable,public :: leapseconds_file  !! for mkspk: LEAPSECONDS_FILE
    character(len=:),allocatable,public :: segment_id        !! for mkspk: SEGMENT_ID
    integer,public :: polynom_degree = 9                     !! for mkspk: POLYNOM_DEGREE
    integer,public :: output_spk_type = 9                    !! for mkspk: OUTPUT_SPK_TYPE
    real(wp),public :: min_export_time_step = 0.0_wp         !! minimum time (sec) step for exporting to mkspk file
                                                             !! [to try and avoid spice interpolation issues for very small time steps]

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
