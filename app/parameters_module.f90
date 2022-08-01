    module parameters_module

    use iso_fortran_env,   only: wp => real64
    use numbers_module,    only: zero

    implicit none

    public

    public :: wp

    ! constants:
    real(wp),parameter :: mu_earth = 3.9860043543609593E+05_wp !! \( \mu_{earth} \)
    real(wp),parameter :: mu_moon  = 4.9028000661637961E+03_wp !! \( \mu_{moon} \)
    real(wp),parameter :: mu_sun   = 1.3271244004193938E+11_wp !! \( \mu_{sun} \)
    real(wp),parameter :: r_moon   = 1737.4_wp  !! radius of hte Moon

    ! variables loaded from the config file:
    real(wp) :: rp = 0.0_wp                     !! peripais radius
    character(len=:),allocatable :: N_or_S      !! 'N' or 'S'
    character(len=:),allocatable :: L1_or_L2    !! 'L1' or 'L2'

    logical :: generate_plots = .true.  !! to generate the python plots
    logical :: generate_trajectory_files = .true.  !! to export the txt trajectory files.

    integer :: year   = 0 !! epoch of first point (first periapsis crossing)
    integer :: month  = 0 !! epoch of first point (first periapsis crossing)
    integer :: day    = 0 !! epoch of first point (first periapsis crossing)
    integer :: hour   = 0 !! epoch of first point (first periapsis crossing)
    integer :: minute = 0 !! epoch of first point (first periapsis crossing)
    real(wp) :: sec   = 0.0_wp !! epoch of first point (first periapsis crossing)
    real(wp) :: et_ref = zero  !! ephemeris time reference epoch
                               !! (all times relative to this)
                               !! this is computed from the year,month,day,hour,minute,sec values
                               !! [sec]

    integer :: n_revs = 10  !! Number of revs in the Halo.

    character(len=:),allocatable :: ephemeris_file !! the JPL ephemeris file to load
    ! [note: these are build using the get_third_party script in FAT.
    !  the files are platform (and maybe compiler specific)
    !  'data/eph/JPLEPH.421'               !! JPL DE421 ephemeris file [mac,gfortran]
    !  'data/eph/JPLEPH_intel_mac.421'     !! [mac,ifort]
    !  'data/eph/JPLEPH_windows_ifort.421' !! [windows,ifort]

    ! environment and integrator parameters:

    character(len=:),allocatable :: gravfile  !! spherical harmonic gravity coeff file (Moon)
    !! example: 'data/grav/gggrx_0020pm_sha.tab'

    character(len=:),allocatable :: patch_point_file   !! Halo CR3BP patch point solution file
    !! example: 'data/L2_halos.json'

    integer,parameter  :: n_eoms = 6            !! size of EOM derivative vector [x,y,z,vx,vy,vz]
    real(wp),parameter :: rtol = 1.0e-10_wp     !! integrator tols
    real(wp),parameter :: atol = 1.0e-10_wp     !! integrator tols
    integer,parameter  :: maxnum = 10000        !! integrator max steps
    integer,parameter  :: grav_n = 8            !! max grav degree
    integer,parameter  :: grav_m = 8            !! max grav order

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

    ! the initial guess for the Halo:
    ! (4500 L2 S case)
    !
    ! Note: this should be read from the patch point file... for now, just hard code it
    ! This is for the "three_patch_point_solution_forward" version of the transcription.

    type(patch_point),public :: periapsis    !! Patchpoint Periapse State
    type(patch_point),public :: quarter      !! Patchpoint 1/4 Rev State
    type(patch_point),public :: apoapsis     !! Patchpoint 1/2 Rev State

    real(wp)        :: period  = 0.0_wp !! Halo period [days]
    real(wp),public :: period8 = 0.0_wp !! 1/8 of Halo period [days]

    ! fix these for now (good values for the 4500 Rp case)
    ! TODO ... set the xscale based on the actual values ...
    real(wp),dimension(6),parameter,public :: xscale_x0 = &
        [1.0e+05_wp,1.0e+05_wp,1.0e+05_wp,2.0e+00_wp,2.0e+00_wp,2.0e+00_wp]
    real(wp),dimension(6),parameter,public :: fscale_xf = &
        [1.0e+04_wp,1.0e+04_wp,1.0e+04_wp,1.0e+02_wp,1.0e+02_wp,1.0e+02_wp]

    logical,public  :: fix_initial_time = .false. !! to fix the initial epoch in the mission

    end module parameters_module
