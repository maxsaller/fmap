module variables

    implicit none

    ! CONSTANTS
    integer, parameter :: S = 2                              ! # of el.st.
    double complex, parameter :: eye = dcmplx(0.d0, 1.d0)    ! Complex i
    double precision, parameter :: pi = 3.14159265359d0      ! Pie is delicious
    double precision, parameter :: sol = 137.03599908381668d0! Speed of light
    double precision, parameter :: eps0 = 0.079577471546d0   ! Vacuum permittivity

    ! PARALLEL PARAMETERS
    integer :: threads                          ! Number of max OMP threads

    ! INPUT PARAMETERS
    integer :: F                                ! Nuclear DoFs
    integer :: ntraj                            ! Number of trajectories
    integer :: tsteps                           ! Number of time-steps
    integer :: cavitysteps                      ! Number of cavity spatial steps
    double precision :: epsilon                 ! Energy bias
    double precision :: coupling                ! Coupling matrix element
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
    double precision, allocatable :: Cpop(:,:,:)! Population CF
    !double precision, allocatable :: Cmix(:,:,:)! Mixed population CF

    ! BATH OBSERVABLES

    ! LAPACK PARAMETERS
    integer :: info, lenwork                    ! Integer parameters for LAPACK
    double precision, allocatable :: work(:)    ! Work array for LAPACK

end module variables
