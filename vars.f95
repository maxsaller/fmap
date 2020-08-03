module variables

    implicit none

    ! CONSTANTS
    double complex, parameter :: eye = dcmplx(0.d0, 1.d0)   ! Complex i
    double precision, parameter :: pi = 3.14159265359d0     ! Pi is delicious
    double precision, parameter :: sol = 137.035999084d0    ! Speed of light
    double precision, parameter :: eps0 = 0.079577471546d0  ! Vacuum permittivity

    ! SYSTEM PROPERTIES
    integer, parameter :: S = 2                             ! Electronic states
    double precision, parameter :: mu = 1.034d0             ! Dipole moment                       
    double precision, parameter :: t0 = 0.197d0             ! KE matrix element
    double precision, parameter :: L = 236215.765578d0      ! Cavity length

    ! BATH VARIABLES
    double precision, allocatable :: omega(:)       ! Bath frequencies
    double precision, allocatable :: lambda(:)      ! Electron-photon couplings
    double precision, allocatable :: lam_mat(:,:)   ! Lambda matrix

    ! INPUT PARAMETERS
    integer :: F                                ! Nuclear DoFs
    integer :: ntraj                            ! Number of trajectories
    integer :: tsteps                           ! Number of time-steps
    integer :: intensity_points                 ! Number of intensity snapshots
    integer :: cavitysteps                      ! Number of cavity spatial steps
    double precision :: dt                      ! Time-step duration
    character(len=6) :: Aop                     ! Electronic sampling type
    character(len=6) :: Bop                     ! Electronic sampling type
    character(len=4) :: electronic              ! Electronic sampling type
    character(len=14):: intgt                   ! Type of integrator to use

    ! TRAJECTORY VARIABLES
    double precision, allocatable :: xn(:)      ! Nuclear position
    double precision, allocatable :: pn(:)      ! Nuclear momentum
    double precision, allocatable :: XE(:)      ! Mapping variable position
    double precision, allocatable :: PE(:)      ! Mapping variable momentum

    ! POTENTIAL AND FORCE MATRICES
    double precision :: V0                      ! State-independent potential
    double precision, allocatable :: V(:,:)     ! Potential energy matrix
    double precision, allocatable :: G0(:)      ! State-independent force

    ! LAPACK PARAMETERS
    integer :: info, lenwork                    ! Integer parameters for LAPACK
    double precision, allocatable :: work(:)    ! Work array for LAPACK

    ! OBSERVABLES AND OPERATORS
    double precision :: pop_norm                ! Normalisation for mapp. pop.
    double precision :: pop_0(S)                ! Time-zero mapping population
    double precision :: pop_t(S)                ! Time-t mapping population
    double precision, allocatable :: Cpop(:,:,:)! Mapping population c. f.
    double precision, allocatable :: Cpopm(:,:,:)! Mapp. pop. c. f. mixed states
    double precision, allocatable :: EFI(:,:,:) ! Electric field intensity

end module variables
