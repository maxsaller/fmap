module variables

    implicit none

    ! CONSTANTS
    double complex, parameter :: eye = dcmplx(0.d0, 1.d0)
    double precision, parameter :: pi = 3.14159265359d0
    double precision, parameter :: sol = 137.03599908381668d0

    ! PARALLEL PARAMETERS
    integer :: threads                          ! Number of max OMP threads

    ! INPUT PARAMETERS
    integer :: F                                ! Nuclear DoFs
    integer :: S                                ! Electronic states
    integer :: ntraj                            ! Number of trajectories
    integer :: tsteps                           ! Number of time-steps
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
    double precision :: V0                       ! State-independent potential
    double precision, allocatable :: V(:,:)      ! Potential energy matrix
    double precision, allocatable :: G0(:)       ! State-independent force

    ! SYSTEM PROPERTIES
    double precision :: L                        ! Length of cavity
    double precision :: mu                       ! Dipole operator
    double precision, allocatable :: eps(:)      ! State energies

    ! BATH VARIABLES
    double precision, allocatable :: omega(:)   ! Bath frequencies
    double precision, allocatable :: c(:)       ! Electron-phonon couplings

    ! POPULATION OBSERVABLES
    double precision, allocatable :: pop_0(:)   ! Time-zero populations
    double precision, allocatable :: pop_t(:)   ! Time-t populations
    double precision, allocatable :: Qop_0(:)   ! Time-zero improved operator
    double precision, allocatable :: Qop_t(:)   ! Time-t improved operator
    double precision, allocatable :: Cpop(:,:,:)! Population correlation fn.
    double precision, allocatable :: Cimp(:,:,:)! Improved operator corr. fn.

    ! BATH OBSERVABLES
    double precision, allocatable :: zeta(:,:)  ! Utility array for intensity
    double precision, allocatable :: I_pop(:,:) ! Cavity intensity traditional
    double precision, allocatable :: I_imp(:,:) ! Cavity intensity improved
    double precision, allocatable :: NP_pop(:,:)! Number of photons per mode trad.
    double precision, allocatable :: NP_imp(:,:)! Number of photons per mode impr.
    double precision, allocatable :: SNP_pop(:) ! Total number of photons trad.
    double precision, allocatable :: SNP_imp(:) ! Total number of photons impr.

    ! LAPACK PARAMETERS
    integer :: info, lenwork                    ! Integer parameters for LAPACK
    double precision, allocatable :: work(:)    ! Work array for LAPACK

end module variables
