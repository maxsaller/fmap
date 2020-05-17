module variables

    implicit none

    ! CONSTANTS
    double precision, parameter :: pi = 3.14159265359d0
    double complex, parameter :: eye = dcmplx(0.d0, 1.d0)

    ! PARALLEL PARAMETERS
    integer :: threads                          ! Number of max OMP threads

    ! INPUT PARAMETERS
    integer :: F                                ! Nuclear DoFs
    integer :: S                                ! Electronic states
    integer :: ntraj                            ! Number of trajectories
    integer :: tsteps                           ! Number of time-steps
    double precision :: dt                      ! Time-step duration
    double precision :: beta                    ! Inverse temperature
    double precision :: kondo                   ! Bath kondo parameter
    double precision :: delta                   ! Constant electronic coupling
    double precision :: omegac                  ! Bath cutoff frequency
    double precision :: epsilon                 ! Electronic energy bias
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
    double precision, allocatable :: G(:,:,:)    ! Force tensor

    ! BATH VARIABLES
    double precision, allocatable :: omega(:)   ! Bath frequencies
    double precision, allocatable :: c(:)       ! Bath coupling coefficients

    ! OBSERVABLES
    double precision, allocatable :: pop_0(:)   ! Time-zero populations
    double precision, allocatable :: pop_t(:)   ! Time-t populations
    double precision, allocatable :: Qop_0(:)   ! Time-zero improved operator
    double precision, allocatable :: Qop_t(:)   ! Time-t improved operator
    double precision, allocatable :: Cpop(:,:,:)! Population correlation fn.
    double precision, allocatable :: Cimp(:,:,:)! Improved operator corr. fn.

    ! LAPACK PARAMETERS
    integer :: info, lenwork                    ! Integer parameters for LAPACK
    double precision, allocatable :: work(:)    ! Work array for LAPACK

end module variables
