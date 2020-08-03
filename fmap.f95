! Simple correlation functions via mapping approach methods !
! including LSC-IVR, PBME and a number of mLSC approaches.  !
! Integration via diagonalisation or velocity verlet.       !
! Copyright: Max Saller 2020                                !
program fmap

    use variables
    implicit none
    integer :: t,ts,pctg,id,i,count
    character(len=1) :: cr = char(13)

    ! HEADER
    write(6,"(a50,//,23x,a4,23x,//,2x,a46,2x,//,a50,/)") repeat("#",50),&
    "FMAP","Correlation functions via the mapping approach",repeat("#",50)

    ! INITIALIZE RANDOM NUMBER GENERATOR FROM SYSTEM CLOCK
    call system_clock(count)
    call RLUXGO(4,count,0,0)

    ! READ INPUT FILE
    call read_input()

    ! ALLOCATE ARRAYS
    call allocate_arrays()

    ! SET UP BATH
    call system_bath_properties()

    ! MONTE CARLO LOOP
    do t = 1,ntraj

        ! REPORT TRAJECTORY PROGRESS
        if ( ntraj > 100 .and. t > 1 .and. mod(t,(ntraj/100)) == 0 ) then
            pctg = floor(1.d2*t/dble(ntraj))
            write(6, "(a22,i3,a1,1x,52a)", advance="no") &
            "RUNNING TRAJECTORIES: ", pctg, "%", &
            "[", repeat("#", pctg/2), repeat(" ", 50-pctg/2), "]", cr
            flush(6)
        end if

        ! SAMPLE INITIAL CONDITIONS
        call sample_nuclear()
        call sample_electronic()

        ! CALCULATE TIME ZERO OPERATORS (AND TIME-T OPERATORS AT T=0)
        call time_zero_ops()
        call time_t_ops()
        call accumulate_obs(1)
        call calc_field_intensity(1)

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
            call accumulate_obs(ts+1)
            if ( ts > 1 .and. mod(ts, tsteps/intensity_points) == 0) then
                call calc_field_intensity(ts+1)
            end if

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
    character(len=50) :: dum

    open(11, file="input", action="read")

    read(11, '(A)') dum
    read(11, *) dum, F
    read(11, '(A)') dum
    read(11, *)

    read(11, '(A)') dum
    read(11, *) dum, ntraj
    read(11, *) dum, tsteps
    read(11, *) dum, dt
    read(11, *) dum, intensity_points
    read(11, *) dum, cavitysteps
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

    if ( mod(tsteps, intensity_points) /= 0 ) then
        write(6,*) "ERROR: tsteps/intensity_points should be an integer!"
        stop
    end if

end subroutine read_input


! Allocates all relevant arrays (and zeroes where required)
subroutine allocate_arrays()

    use variables
    implicit none
    double precision :: evec(S,S), eval(S)

    ! LAPACK WORK ARRAY
    allocate( work(1) )
    call dsyev("V", "U", S, evec, S, eval, work, -1, info)
    lenwork = int(work(1))
    deallocate(work)
    allocate( work(lenwork) )

    ! TRAJECTORY ARRAYS
    allocate( xn(F) )
    allocate( pn(F) )
    allocate( XE(S) )
    allocate( PE(S) )

    ! BATH ARRAYS
    allocate( omega(F) )
    allocate( lambda(F) )
    allocate( lam_mat(F, cavitysteps) )

    ! POTENTIAL AND FORCE MATRICES
    allocate( G0(F) )
    allocate( V(S,S) )

    ! OBSERVABLES ARRAYS
    allocate( Cpop(tsteps+1,S,S) )
    allocate( Cpopm(tsteps+1,S,2*S) )
    allocate( EFI(S,intensity_points+1,cavitysteps) )

    ! ZERO ARRAYS
    Cpop(:,:,:) = 0.d0
    Cpopm(:,:,:) = 0.d0
    EFI(:,:,:) = 0.d0

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
    deallocate( omega )
    deallocate( lambda )
    deallocate( lam_mat )

    ! POTENTIAL AND FORCE MATRICES
    deallocate( V )
    deallocate( G0 )

    ! OBSERVABLES ARRAYS
    deallocate( EFI )
    deallocate( Cpop )
    deallocate( Cpopm )

end subroutine deallocate_arrays


! Sets bath frequencies and coupling coefficients
subroutine system_bath_properties()

    use variables
    implicit none
    integer :: i,j
    double precision :: r, pf
    
    open(11, file="freq.out", status="unknown", action="write")

    pf = sqrt(2.d0 / L / eps0)

    ! FIELD FREQUENCIES AND ELECTRON-PHONON COUPLINGS (AT L/2)
    do i = 1,F
        omega(i) = pi * sol * dble(2 * i - 1) / L
        lambda(i) = pf * (-1)**(i+1)
        write(11, *) omega(i), lambda(i), lambda(i)*omega(i)*mu
    end do

    close(11)

    ! DISCRETIZED LAMBDA MATRIX
    do i = 1, cavitysteps
        r = (i-1) * L / dble(cavitysteps-1)
        do j = 1,F
            lam_mat(j,i) = sin( omega(j) / sol * r )
        end do
    end do
    lam_mat(:,:) = ( lam_mat(:,:) * pf )**2

end subroutine system_bath_properties


! Samples nuclear positions and momenta
subroutine sample_nuclear()

    use variables
    implicit none
    integer :: i
    double precision :: r1, r2, xn_stdev, pn_stdev

    do i = 1,F
        xn_stdev = dsqrt( 1.d0 / (2.d0 * omega(i)) ) 
        pn_stdev = dsqrt( omega(i) / 2.d0 )
        call gauss(r1, r2)
        xn(i) = xn_stdev * r1
        pn(i) = pn_stdev * r2
    end do

end subroutine sample_nuclear


! Samples mapping variable positions and momenta as well as normalisation const.
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
        pop_norm = 4.d0
    else if ( electronic == "phi2" ) then
        XE_stdev = 1.d0 / 2.d0
        PE_stdev = 1.d0 / 2.d0 
        do i = 1,S
            call gauss(r1, r2)
            XE(i) = XE_stdev * r1
            PE(i) = PE_stdev * r2
        end do
        pop_norm = 16.d0
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

end subroutine time_zero_ops


! Calculates time-t operators
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

end subroutine time_t_ops


! Calculate eletric field intensity
subroutine calc_field_intensity(ts)

    use variables
    implicit none
    integer :: i,j,index
    integer, intent(in) :: ts
    double precision :: zpe, Esq

    index = 1 + (ts-1) / (tsteps/intensity_points)

    do i = 1,cavitysteps
        do j = 1, F
            Esq = omega(j)**2*lam_mat(j, i)*xn(j)**2 - 0.5d0*omega(j)*lam_mat(j, i)
            
            EFI(1,index,i) = EFI(1,index,i) + pop_norm * ( pop_t(1) + pop_t(2) ) * &
                             ( pop_0(1) + pop_0(2) )/sqrt(2.d0) * Esq

            EFI(2,index,i) = EFI(2,index,i) + pop_norm * ( pop_t(1) + pop_t(2) ) * &
                             ( pop_0(1) - pop_0(2) )/sqrt(2.d0) * Esq
                             
        end do
    end do

end subroutine calc_field_intensity

! Accumulates observables
subroutine accumulate_obs(ts)

    use variables
    implicit none
    integer :: i,j
    integer, intent(in) :: ts
    double precision :: oost

    oost = 1.d0 / sqrt(2.d0)

    ! TRADITIONAL POPULATION CORRELATION FUNCTION
    do i = 1,S
        do j = 1,S
            Cpop(ts,i,j) = Cpop(ts,i,j) + pop_norm * pop_0(i) * pop_t(j)
        end do
    end do

    ! MIXED STATE POPULATION CORRELATION FUNCTIONS
    Cpopm(ts,1,1) = Cpopm(ts,1,1) + pop_norm * oost * (pop_0(1) + pop_0(2)) * &
                                                       pop_t(1)
    Cpopm(ts,1,2) = Cpopm(ts,1,2) + pop_norm * oost * (pop_0(1) + pop_0(2)) * &
                                                       pop_t(2)
    Cpopm(ts,2,1) = Cpopm(ts,2,1) + pop_norm * oost * (pop_0(1) - pop_0(2)) * &
                                                       pop_t(1)
    Cpopm(ts,2,2) = Cpopm(ts,2,2) + pop_norm * oost * (pop_0(1) - pop_0(2)) * &
                                                       pop_t(2)
    
    Cpopm(ts,1,3) = Cpopm(ts,1,3) + pop_norm * oost * (pop_0(1) + pop_0(2)) * &
                                               oost * (pop_t(1) + pop_t(2))
    Cpopm(ts,1,4) = Cpopm(ts,1,4) + pop_norm * oost * (pop_0(1) + pop_0(2)) * &
                                               oost * (pop_t(1) - pop_t(2))
    Cpopm(ts,2,3) = Cpopm(ts,2,3) + pop_norm * oost * (pop_0(1) - pop_0(2)) * &
                                               oost * (pop_t(1) + pop_t(2))
    Cpopm(ts,2,4) = Cpopm(ts,2,4) + pop_norm * oost * (pop_0(1) - pop_0(2)) * &
                                               oost * (pop_t(1) - pop_t(2))

end subroutine accumulate_obs


! Averages and outputs observable arrays
subroutine average_obs()

    use variables
    implicit none
    integer :: i,j,k
    character(len=120) :: fmt

    write(6,"(//'AVERAGING OBSERVABLES:')")
    
    ! MAPPING POPULATION CORRELATION FUNCTION
    open(11, file="Cpop.out", action="write", status="unknown")
    fmt = '(F10.4'//repeat(',2x,ES13.5',S**2)//')'
    write(11,'("# time C_11(t) C_12(t) C_21(t) C_22(t)")')
    do i = 1, tsteps+1
        write(11,fmt) (i-1)*dt, &
                     Cpop(i,1,1)/dble(ntraj), &
                     Cpop(i,1,2)/dble(ntraj), &
                     Cpop(i,2,1)/dble(ntraj), &
                     Cpop(i,2,2)/dble(ntraj)
    end do
    write(6,"('- Wrote mapping populations to Cpop.out')")
    close(11)

    ! MIXED STATE  POPULATION CORRELATION FUNCTION
    open(11, file="Cpopm.out", action="write", status="unknown")
    fmt = '(F10.4'//repeat(',2x,ES13.5',2*S**2)//')'
    write(11,'("# time C_11(t) C_12(t) C_21(t) C_22(t)")')
    do i = 1, tsteps+1
        write(11,fmt) (i-1)*dt, &
                      Cpopm(i,1,1)/dble(ntraj), &
                      Cpopm(i,1,2)/dble(ntraj), &
                      Cpopm(i,1,3)/dble(ntraj), &
                      Cpopm(i,1,4)/dble(ntraj), &
                      Cpopm(i,2,1)/dble(ntraj), &
                      Cpopm(i,2,2)/dble(ntraj), &
                      Cpopm(i,2,3)/dble(ntraj), &
                      Cpopm(i,2,4)/dble(ntraj)
    end do
    write(6,"('- Wrote mixed state mapping populations to Cpopm.out')")
    close(11)

    ! Electric field intensity
    do k = 1, S
        write(fmt,'(I1)') k
        open(11, file="EFI_init"//trim(fmt)//".out", action="write", status="unknown")
        do i = 0, intensity_points
            write(11,'(F10.4)', advance="no") i*tsteps/intensity_points*dt
            do j = 1,cavitysteps
                write(11, '(2x,ES13.5)', advance='no') EFI(k,i+1,j)/dble(ntraj)
            end do
            write(11,'(1x)')
        end do
        close(11)
    end do
    write(6,"('- Wrote electric field intensities to EFI_initX.out')")


end subroutine average_obs


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
        pn(i) = pn(i) - qdt * lambda(i) * omega(i) * mu * &
                (XE(1)**2 - XE(2)**2 + PE(1)**2 - PE(2)**2)
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
        pn(i) = pn(i) - qdt * lambda(i) * omega(i) * mu * &
                (XE(1)**2 - XE(2)**2 + PE(1)**2 - PE(2)**2)
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
        pn(i) = pn(i) - qdt * lambda(i) * omega(i) * mu * &
                (XE(1)**2 - XE(2)**2 + PE(1)**2 - PE(2)**2)
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
        pn(i) = pn(i) - qdt * lambda(i) * omega(i) * mu * &
                (XE(1)**2 - XE(2)**2 + PE(1)**2 - PE(2)**2)
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
    V(:,:) = 0.d0
    V(1,2) = -t0
    V(2,1) = -t0
    do i = 1,F
        V(1,1) = V(1,1) + omega(i) * lambda(i) * mu * xn(i)
        V(2,2) = V(2,2) - omega(i) * lambda(i) * mu * xn(i)
    end do

end subroutine potential_force


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
