! Simple correlation functions via mapping approach methods !
! including LSC-IVR, PBME and a number of mLSC approaches.  !
! Integration via diagonalisation or velocity verlet.       !
! Copyright: Max Saller 2020                                !
program fmap

    use variables
    implicit none
    integer :: t,ts,pctg
    character(len=1) :: cr = char(13)

    ! HEADER
    write(6,"(a50,//,23x,a4,23x,//,2x,a46,2x,//,a50,/)") repeat("#",50),&
    "FMAP","Correlation functions via the mapping approach",repeat("#",50)

    ! READ INPUT FILE
    call read_input()

    ! ALLOCATE ARRAYS
    call allocate_arrays()

    ! SET UP BATH
    call spectral_density()

    ! MONTE CARLO LOOP
    do t = 1,ntraj

        ! REPORT TRAJECTORY PROGRESS
        if ( t > 1 .and. mod(t,(ntraj/100)) == 0 ) then
            pctg = floor(1.d2*t/dble(ntraj))
            write(6, "(a22,i3,a1,1x,52a)", advance="no") &
            "RUNNING TRAJECTORIES: ", pctg, "%", &
            "[", repeat("#", pctg/2), repeat(" ", 50-pctg/2), "]", cr
            flush(6)
        end if

        ! SAMPLE INITIAL CONDITIONS
        call sample_nuclear()
        call sample_electronic()

        ! CALCULATE TIME ZERO OPERATORS
        call time_zero_ops()

        ! TRAJECTORY LOOP
        do ts = 1,tsteps

            ! MAKE DYNAMICS STEP
            if ( intgt == "diagonalise" .or. intgt == "diag" ) then
                call step_diag()
            else if ( intgt == "velocityverlet" .or. intgt == "vv") then
                call step_vverlet()
            else
                write(6,*) "ERROR: integrator type must be either ",&
                           "'diagonalise'/'diag' or 'velocityverlet'/'vv'!"
                stop
            end if

            ! CALCULATE TIME t OPERATORS AND ACCUMULATE OBSERVABLES
            call time_t_ops()
            call accumulate_obs(ts)

        end do

    end do

    ! AVERAGE AND OUTPUT OBSERVABLES
    call average_obs() 

    ! DEALLOCATE ARRAYS
    call deallocate_arrays()

    ! FINISH UP
    write(6,"(/'SIMULATION COMPLETE!')")

end program fmap


! Reads the "input" file and prints input arguments to log
subroutine read_input()

    use variables
    implicit none
    integer :: iost=1
    character(len=80) :: dum

    open(11, file="input", action="read")

    read(11, '(A)') dum
    read(11, *) dum, S
    read(11, *) dum, F
    read(11, *) dum, epsilon
    read(11, *) dum, delta
    read(11, *) dum, kondo
    read(11, *) dum, omegac
    read(11, *) dum, beta
    read(11, '(A)') dum
    read(11, *)

    read(11, '(A)') dum
    read(11, *) dum, ntraj
    read(11, *) dum, tsteps
    read(11, *) dum, dt
    read(11, *) dum, intgt
    read(11, '(A)') dum
    read(11, *)

    read(11, '(A)') dum
    read(11, *) dum, Aop
    read(11, *) dum, Bop
    read(11, *) dum, electronic
    read(11, '(A)') dum

    rewind(11)
    write(6,"('INPUT ARGUMENTS PARSED:'/)")
    do while ( .true. )
        read(11,'(A)',iostat=iost) dum
        if ( iost < 0 ) exit
        write(6,*) "   ", dum
    end do
    write(6,*)

    close(11)

end subroutine read_input


! Allocates all relevant arrays (and zeroes where required)
subroutine allocate_arrays()

    use variables
    implicit none

    ! LAPACK WORK ARRAY
    allocate( work(1) )
    call dsyev("V", "U", S, 0.d0, S, 0.d0, work, -1, info)
    lenwork = int(work(1))
    deallocate(work)
    allocate( work(lenwork) )

    ! TRAJECTORY ARRAYS
    allocate( xn(F) )
    allocate( pn(F) )
    allocate( XE(S) )
    allocate( PE(S) )

    ! BATH ARRAYS
    allocate( c(F) )
    allocate( omega(F) )

    ! POTENTIAL AND FORCE MATRICES
    allocate( G0(F) )
    allocate( V(S,S) )
    allocate( G(F,S,S) )

    ! OBSERVABLE ARRAYS
    allocate( pop_0(S) )
    allocate( pop_t(S) )
    allocate( Qop_0(S) )
    allocate( Qop_t(S) )
    allocate( Cpop(tsteps,S,S) )
    allocate( Cimp(tsteps,S,S) )
    Cpop(:,:,:) = 0.d0
    Cimp(:,:,:) = 0.d0

end subroutine allocate_arrays


! Deallocates all relevant arrays
subroutine deallocate_arrays()

    use variables
    implicit none

    ! LAPACK WORK ARRAY
    deallocate(work)

    ! TRAJECTORY ARRAYS
    deallocate( xn )
    deallocate( pn )
    deallocate( XE )
    deallocate( PE )

    ! BATH ARRAYS
    deallocate( c )
    deallocate( omega )

    ! POTENTIAL AND FORCE MATRICES
    deallocate( V )
    deallocate( G )
    deallocate( G0 )

    ! OBSERVABLE ARRAYS
    deallocate( Cpop )
    deallocate( Cimp )
    deallocate( pop_0 )
    deallocate( pop_t )
    deallocate( Qop_0 )
    deallocate( Qop_t )

end subroutine deallocate_arrays


! Sets bath frequencies and coupling coefficients
subroutine spectral_density()

    use variables
    implicit none
    integer :: i
    double precision :: eta

    eta = 0.5d0 * pi * kondo

    ! Ohmic bath, discretization a la Craig&Mano
    do i = 1,F
        omega(i) = -omegac * log( (i-0.5d0) / dble(F) )
        c(i) = omega(i) * dsqrt( (2 * eta * omegac) / (pi * F) )
    end do

end subroutine spectral_density

! Samples nuclear positions and momenta
subroutine sample_nuclear()

    use variables
    implicit none
    integer :: i
    double precision :: tv, r1, r2, xn_stdev, pn_stdev

    do i = 1,F
        tv = tanh(0.5d0 * beta * omega(i))
        xn_stdev = dsqrt( 1.d0 / (2.d0 * omega(i) * tv) ) 
        pn_stdev = dsqrt( omega(i) / (2.d0 * tv) )
        call gauss(r1, r2)
        xn(i) = xn_stdev * r1
        pn(i) = pn_stdev * r2
    end do

end subroutine sample_nuclear


! Samples mapping variable positions and momenta
subroutine sample_electronic()

    use variables
    implicit none
    integer :: i
    double precision :: r1, r2, XE_stdev, PE_stdev

    if ( electronic == "phi" ) then
        XE_stdev = 1.d0 / dsqrt(2.d0) 
        PE_stdev = 1.d0 / dsqrt(2.d0) 
        do i = 1,S
            call gauss(r1, r2)
            XE(i) = XE_stdev * r1
            PE(i) = PE_stdev * r2
        end do
    else if ( electronic == "phi2" ) then
        XE_stdev = 1.d0 / 2.d0
        PE_stdev = 1.d0 / 2.d0 
        do i = 1,S
            call gauss(r1, r2)
            XE(i) = XE_stdev * r1
            PE(i) = PE_stdev * r2
        end do
    else
        write(6,*) "ERROR: electronic_sampling must be either 'phi' or 'phi2'!"
        stop
    end if

end subroutine sample_electronic


! Calculates time-zero operators
subroutine time_zero_ops()

    use variables
    implicit none
    integer :: i,j
    double precision :: zpe

    ! TRADITIONAL POPULATION OPERATORS
    if ( Aop == "seo" ) then
        zpe = 0.5d0
    else if ( Aop == "wigner" ) then
        zpe = 1.d0
    else
        write(6,*) "ERROR: A-operator type must be 'seo' or 'wigner'!"
        stop
    end if

    do i = 1,S
        pop_0(i) = 0.5d0 * ( XE(i)**2 + PE(i)**2 - zpe )
    end do

    ! IMPROVED POPULATION OPERATORS
    do i = 1,S
        Qop_0(i) = 0.5d0 * dble(S-1) * ( XE(i)**2 + PE(i)**2 )
        do j = 1,S
            if ( j /= i ) then
                Qop_0(i) = Qop_0(i) - 0.5d0 * ( XE(j)**2 + PE(j)**2 )
            end if
        end do
    end do

end subroutine time_zero_ops


! Calculates time-t opeartors
subroutine time_t_ops()

    use variables
    implicit none
    integer :: i,j
    double precision :: zpe

    ! TRADITIONAL POPULATION OPERATORS
    if ( Bop == "seo" ) then
        zpe = 0.5d0
    else if ( Bop == "wigner" ) then
        zpe = 1.d0
    else
        write(6,*) "ERROR: B-operator type must be 'seo' or 'wigner'!"
        stop
    end if

    do i = 1,S
        pop_t(i) = 0.5d0 * ( XE(i)**2 + PE(i)**2 - zpe )
    end do

    ! IMPROVED POPULATION OPERATORS
    do i = 1,S
        Qop_t(i) = 0.5d0 * dble(S-1) * ( XE(i)**2 + PE(i)**2 )
        do j = 1,S
            if ( j /= i ) then
                Qop_t(i) = Qop_t(i) - 0.5d0 * ( XE(j)**2 + PE(j)**2 )
            end if
        end do
    end do

end subroutine time_t_ops


! Makes a single trajectory step using velocity verlet
subroutine step_vverlet()

    use variables
    implicit none
    integer :: i,j
    double precision :: hdt, qdt

    hdt = 0.5d0*dt
    qdt = 0.5d0*hdt

    call potential_force()

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

    ! HALF STEP IN MAPPING MOMENTA
    PE = PE - hdt * matmul(V, XE)

    ! FULL STEP IN NUCLEAR POSITIONS
    do i = 1,F
        xn(i) = xn(i) + dt * pn(i)
    end do

    ! FULL STEP IN MAPPING POSITIONS
    call potential_force()
    XE = XE + dt * matmul(V, PE)

    ! HALF STEP IN MAPPING MOMENTA
    call potential_force()
    PE = PE - hdt * matmul(V, XE)

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

end subroutine step_vverlet


! Makes a single trajectory step using diagonalisation of the potential matrix
subroutine step_diag

    use variables
    implicit none
    integer :: i,j
    double complex :: propagator(S) 
    double precision :: hdt, qdt, eval(S), evec(S,S), XEU(S), PEU(S)

    hdt = 0.5d0*dt
    qdt = 0.5d0*hdt

    call potential_force()

    ! HALF STEP IN MAPPING VARIABLES
    evec(:,:) = V(:,:)
    call dsyev("V", "U", S, evec, S, eval, work, lenwork, info)
    XEU = matmul(XE, evec)
    PEU = matmul(PE, evec)
    propagator = exp(-eye * hdt * eval) * dcmplx(XEU, PEU)
    XE = matmul(evec, real(propagator))
    PE = matmul(evec, aimag(propagator))

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

    ! FULL STEP IN NUCLEAR POSITIONS
    do i = 1,F
        xn(i) = xn(i) + dt * pn(i)
    end do

    ! RECALCULATE POTENTIAL AND FORCE WITH NEW NUCLEAR POSITIONS
    call potential_force()

    ! HALF STEP IN NUCLEAR MOMENTA
    do i = 1,F
        pn(i) = pn(i) - hdt * G0(i)
        do j = 1,S
            pn(i) = pn(i) - qdt * (G(i,j,j)*(XE(j)**2 + PE(j)**2))
        end do
    end do

    ! HALF STEP IN MAPPING VARIABLES
    evec(:,:) = V(:,:)
    call dsyev("V", "U", S, evec, S, eval, work, lenwork, info)
    XEU = matmul(XE, evec)
    PEU = matmul(PE, evec)
    propagator = exp(-eye * hdt * eval) * dcmplx(XEU, PEU)
    XE = matmul(evec, real(propagator))
    PE = matmul(evec, aimag(propagator))

end subroutine step_diag


! Calculates potential energy matrices
subroutine potential_force()

    use variables
    implicit none
    integer :: i,j
    double precision :: tr, trace

    ! STATE-INDEPENDENT POTENTIAL AND FORCE
    V0 = 0.d0
    do i = 1,F
        V0 = V0 + 0.5d0 * omega(i)**2 * xn(i)**2
        G0(i) = omega(i)**2 * xn(i)
    end do

    ! POTENTIAL ENERGY MATRIX AND FORCE TENSOR
    V(1,1) = epsilon
    V(1,2) = delta
    V(2,1) = delta
    V(2,2) = -epsilon
    G(:,:,:) = 0.d0
    do i = 1,F
        V(1,1) = V(1,1) + c(i) * xn(i)
        V(2,2) = V(2,2) - c(i) * xn(i)
        G(i,1,1) = c(i)
        G(i,2,2) = -c(i)
    end do

    ! SHIFT TRACE OF V and G to V0 and G0
    tr = trace(V,S)/dble(S)
    V0 = V0 + tr
    do i = 1,S
        V(i,i) = V(i,i) - tr
    end do
    
    do i = 1,F
        tr = trace(G(i,:,:),S)/dble(S)
        G0(i) = G0(i) + tr
        do j = 1,S
            G(i,j,j) = G(i,j,j) - tr
        end do
    end do

end subroutine potential_force


! Accumulates observables
subroutine accumulate_obs(ts)

    use variables
    implicit none
    integer :: i,j
    integer, intent(in) :: ts
    double precision :: norm

    ! TRADITIONAL POPULATION OPERATORS
    if ( Aop == "seo" .and. Bop == "seo" ) then
        norm = 16.d0
    else if ( Aop == "wigner" .and. Bop == "wigner" ) then
        write(6,*) "ERROR: Having both the A- and B-operator be of type ",&
                   "'wigner' does not make sense! At least one operator ",&
                   "must be projected onto onto the SEO subspace!"
        stop
    else
        norm = 4.d0
    endif

    do i = 1,S
        do j = 1,S
            Cpop(ts,i,j) = Cpop(ts,i,j) + norm * pop_0(i) * pop_t(j)
        end do
    end do

    ! IMPROVED POPULATION OPERATORS
    if ( electronic == "phi" ) then
        norm = 4.d0
    else if ( electronic == "phi2" ) then
        norm = 16.d0
    end if
    do i = 1,S
        do j = 1,S
            Cimp(ts,i,j) = Cimp(ts,i,j) + &
            ( S + norm * Qop_t(j) + norm * Qop_0(i)*Qop_t(j) ) / dble(S**2)
        end do
    end do    

end subroutine accumulate_obs


! Averages and outputs observable arrays
subroutine average_obs()

    use variables
    implicit none
    integer :: i,j,k

    write(6,"(//'AVERAGING OBSERVABLES:')")

    ! TRADITIONAL POPULATION OPERATORS
    open(11, file="Cpop.out", action="write", status="unknown")
    do i = 1,tsteps
        Cpop(i,:,:) = Cpop(i,:,:) / dble(ntraj)
        write(11,'(F10.4,2x,4(ES13.6,2x))') &
        dble(i)*dt, Cpop(i,1,1), Cpop(i,1,2), Cpop(i,2,1), Cpop(i,2,2)
    end do
    close(11)
    write(6,"('- Saved population autocorrelation functions to Cpop.out')")

    ! IMPROVED POPULATION OPERATORS
    open(11, file="Cimp.out", action="write", status="unknown")
    do i = 1,tsteps
        Cimp(i,:,:) = Cimp(i,:,:) / dble(ntraj)
        write(11,'(F10.4,2x,4(ES13.6,2x))') &
        dble(i)*dt, Cimp(i,1,1), Cimp(i,1,2), Cimp(i,2,1), Cimp(i,2,2)
    end do
    close(11)
    write(6,"('- Saved improved population operator corr. fn. to Cimp.out')")

end subroutine average_obs


! Returns two numbers sampled from the standard normal distribution
! which are obtained via the Box-Muller transform of a RANLUX uniform call
subroutine gauss(r1, r2)

    use variables
    implicit none
    real(4) :: yield(2)
    double precision, intent(out) :: r1, r2

    call RANLUX(yield, 2)

    r1 = dsqrt(-2.d0 * log(yield(1))) * cos(2.d0 * pi * yield(2))
    r2 = dsqrt(-2.d0 * log(yield(1))) * sin(2.d0 * pi * yield(2))

end subroutine gauss


! Calculates the trace of a square matrix
double precision function trace(A,D)

    implicit none
    integer :: i
    integer, intent(in) :: D
    double precision, intent(in) :: A(D,D)
    
    trace = 0.d0

    do i = 1, D
        trace = trace + A(i,i)
    end do

end function trace
