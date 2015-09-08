! Copyright (c) 2013, Sandia Corporation.  Under the terms of Contract
! DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
! certain rights in this software. 

! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:

!  o  Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.

!  o  Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.

!  o  Neither the name of Sandia Corporation nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.


! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

PROGRAM CACTUS

    ! CACTUS (Code for Axial and Cross-flow TUrbine Simulation) v 1.0
    ! Lead Author: J. Murray
    ! Contributor: M. Barone
    !
    ! A lifting-line, vortex-lattice wind or water power rotor simulation code
    ! originally based on the legacy SNL VDART3 code (Strickland
    !  et al.)
    !
    ! -----------
    ! Notes: 
    !   Command line inputs:
    !    FILENAME: Full name of the input namelist file, if
    !              located in the calling directory. Full path otherwise.
    !  
    ! Note that all internal variables are normalized. Velocity
    !  variables are normalized by freestream velocity, distance
    !  variables are normalized by reference radius, velocity
    !  derivatives by (freestream velocity)/(reference radius),
    !  circulation variables by (freestream velocity)*(reference radius)

    use parameters
    use util
    use dystl
    use blade
    use wake
    use element
    use varscale
    use shear
    use iecgust
    use airfoil
    use configr
    use pidef
    use vortex
    use wakedata
    use time
    use wallsoln
    use regtest
    use output
    use tower
    use fnames
    use probesystem
    use compiler
!$  use omp_lib

    !IMPLICIT NONE !JCM: eventually...      

    integer :: ErrFlag
    logical :: FinalConv
    logical :: ConvFlag
    logical :: ContinueRevs
    logical :: ContinueNL
    integer :: iConv
    real :: NLTol
    real :: CPAve_last

!$  integer :: nthreads, tid

    write(*,'(A)') 'Starting CACTUS Execution.'
    write(*,'(A)') '--------------------------'

    ! print compiler info to stdout
    call print_compiler_info()
    write(*,*) ''


    ! Alert the user if OpenMP is enabled.
!$  write(*,*) 'OpenMP is Enabled.'

!$omp parallel private(nthreads,tid)
!$  tid = omp_get_thread_num()
!$  nthreads = omp_get_num_threads()
!$  if (tid == 0) then
!$      write(*,*) 'Executing with ', nthreads, ' threads.'
!$  end if
!$omp end parallel

    write(*,*) ''
    write(*,*) ''

    ! Check for correct number of command line arguments
    nargin=command_argument_count()
    if (nargin < 1) then
        write(6,'(A)') 'Please call the program with the name of the input file on the command line. Ex. CACTUS INPUTFILE.in' 
        stop
    end if

    ! Parse command line to get the base of the input filename
    Call get_FNBase

    ! Write the input files to stdout (for bookkeping)
    write(*,*) 'Input file'
    write(*,*) '--------------------------'
    Call file_to_stdout(InputFN)
    write(*,*) ''

    ! Set the filenames for the output files
    SFOutputFN=trim(FNBase)//'_Param.csv'
    RevOutputFN=trim(FNBase)//'_RevData.csv'
    TSOutputFN=trim(FNBase)//'_TimeData.csv'
    ELOutputFN=trim(FNBase)//'_ElementData.csv'  

    ! Pi definition
    pi = 4.0*atan(1.0)
    conrad = pi/180.0                                                  
    condeg = 180.0/pi 

    ! Namelist input file                                                       
    OPEN(4, FILE= InputFN)                                     

    ! Output files                                                      
    OPEN(8, FILE=SFOutputFN) 
    OPEN(9, FILE=RevOutputFN)
    OPEN(10, FILE=TSOutputFN)

    ! Initialize iteration parameters                                                       
    irev=0
    nt = 0
    Theta=0.0
    TimeN=0.0
    ntTerm=1 
    NLTol=1.0e-04                                                      
    CPAve_last=0.0
    CPSum=0.0
    CTRSum=0.0
    CFxSum=0.0
    CFySum=0.0
    CFzSum=0.0

    ! Error flags                                                                                                                                                                  
    ilxtp=0
    iuxtp=0

    ! Read inputs    
    ErrFlag = 0
    CALL input(ErrFlag)
    if (ErrFlag == 1) then
        write(6,'(A)') 'Input data fail. Exiting...'
        stop
    end if


    ! Setup output

    ! Write headers on standard output files
    Call csvwrite(8,Output_SFHead,Output_SFData,1,0)
    Call csvwrite(9,Output_RevHead,Output_RevData,1,0) 
    Call csvwrite(10,Output_TSHead,Output_TSData,1,0)

    ! Simple output for regression testing        
    if (RegTFlag == 1) then      
        RegOutputFN=trim(FNBase)//'_RegData.out'
        OPEN(7, FILE= RegOutputFN,  FORM= 'FORMATTED' )  
        Call WriteRegTOutput(0)
    end if

    ! Optional element load output
    if (Output_ELFlag == 1) then
        OPEN(11, FILE=ELOutputFN)
        Call csvwrite(11,Output_ELHead,Output_ELData,1,0)
    end if

    ! Optional wall model output
    if (WallOutFlag > 0) then
        if (FSFlag == 1) then
            FSOutputFN=trim(FNBase)//'_FSData.csv'
            OPEN(15, FILE=FSOutputFN)
            write(15,'(A)') trim(FSOutHead)
        end if
    end if

    ! Optional dynamic stall diagnostic output
    if (Output_DSFlag == 1) then
        DSOutputFN=trim(FNBase)//'_DSData.csv'
        OPEN(16, FILE=DSOutputFN)
        if (Output_DSType == 1) then
            Call csvwrite(16,Output_BVHead,Output_BVData,1,0)
        else if (Output_DSType == 2) then
            Call csvwrite(16,Output_LBHead,Output_LBData,1,0)
        end if
    end if

    ! If wake update interval set to a negative number, set next wake update iteration to -1 (no wake velocity updates will be performed)
    ! Otherwise, make the first update on the second iteration (when the wake first appears)
    if (iut < 0) then
        nsw=-1
    else
        nsw=2
    end if

    ! Set first wall update timestep
    if (GPFlag == 1 .or. WPFlag == 1 .or. FSFlag == 1) then
        nsWall=1
    end if

    ! Blade and strut geometry setup. 
    ! Global axes: x is oriented with the nominal freestream flow direction, y is in the vertically upward direction (opposite gravity),
    ! z direction from RHR (to the right when looking in the streamwise direction).
    CALL BGeomSetup()
    CALL SGeomSetup()
    if (Itower.EQ.1) Then
       ! Non-dimensionalize tower inputs
       tower_D = tower_D / Rmax
       tower_x = tower_x / Rmax
       tower_ybot = tower_ybot / Rmax
       tower_ytop = tower_ytop / Rmax
       CALL setup_tower()
    End if

    ! Set normalized turbine rotation rate
    wRotX=ut*RotX
    wRotY=ut*RotY  
    wRotZ=ut*RotZ
    delt=2.0*pi/nti ! change in theta per timestep

    ! Setup wall geometry and solution if necessary
    if (GPFlag == 1 .OR. WPFlag == 1 .OR. FSFlag == 1) then
        ! Wall Geometry setup
        Call WGeomSetup()

        ! Setup wall solution
        Call WSolnSetup()
    end if

    ! Normalization parameters for geometry and performance outputs
    romega=2.0*pi*Rmax*rpm/60.0                                       
    uinf=romega/ut 
    DT=DelT/ut                                      ! normalized simulation timestep (dt*Uinf/Rmax)      

    gustT = gusttime * uinf / Rmax
    gustA = gustamp / uinf

    rem=rho*uinf*Rmax/vis                                          

    if (Incompr .eq. 1) then            ! incompressible/compressible flow switch (used by dynamic stall model)
        Minf = 0.0
    else
        Minf=uinf/sqrt(1.4*1716.0*(tempr+459.6))
    end if

    ! Setup vortex core radius for bound vorticity based on max chord, and for trailing and spanwise wake based on 
    ! temporal and spatial discretization levels, respectively.
    if (ivtxcor > 0) then
        vRad_B = CrRef*VCRFB
        vRad_T = dSGeom*VCRFT
        dSWake = delt*(1.0+1.0/max(ut,1.0))   ! representative wake discretization size (wake line at Rmax in Uinf freestream)
        vRad_S = dSWake*VCRFS       
    else
        vRad_B = 0.0         
        vRad_T = 0.0  
        vRad_S = 0.0
    end if
    vRad2_B = vRad_B*vRad_B                             
    vRad2_T = vRad_T*vRad_T
    vRad2_S = vRad_S*vRad_S

    areat=at*Rmax**2                                ! frontal area (at is (frontal area) / Rmax^2 ) 
    dynpress=rho/2.0*uinf**2                        ! dynamic pressure        
    forcec=dynpress*areat                           ! force coeff normalization
    torquec=dynpress*areat*Rmax                     ! torque coeff normalization                         
    powerc=rho/2.0*areat*romega**3*0.7457/550.      ! normalization for power coeff using tip speed (kp), with conversion from lb-ft/s to kW. (Used to write output)

    ! Write flow properties output
    Output_SFData(1,1)=Rmax         ! length scale
    Output_SFData(1,2)=areat        ! frontal area
    Output_SFData(1,3)=rpm          ! turbine RPM
    Output_SFData(1,4)=Uinf         ! freestream U
    Output_SFData(1,5)=rho          ! density
    Output_SFData(1,6)=tempr        ! temperature
    Output_SFData(1,7)=vis          ! viscosity
    Output_SFData(1,8)=dynpress     ! dynamic pressure
    Output_SFData(1,9)=ut           ! tip speed ratio
    Output_SFData(1,10)=rem         ! machine Reynolds number based on U and Rmax
    Output_SFData(1,11)=FnR         ! Froude number based on radius 
    Call csvwrite(8,Output_SFHead,Output_SFData,0,1)

    ! close parameter file for writing
    CLOSE(8)
    
    ! Initialize needed arrays                             
    do i=1,ne                                                      
        gs(1,i)=0.0   
        gb(i)=0.0                                                    
        ogb(i)=0.0 
        AOA(i)=0.0                                                                                                                                                             
    end do
    Call UpdateAOALast(ne)

    ! Initialize dynamic stall model
    if (DSFlag==1) then
        Call dystl_init_BV()
    else if (DSFlag==2) then
        Call dystl_init_LB()
    end if

    ! CPU time markers
    Call cpu_time(t0)                                                                                                           
!$  t0 = omp_get_wtime()
    Time1=t0                                                                             

    write(*,*) 'Simulation Status'
    write(*,*) '--------------------------'

    ! Do revolutions until convergence or MaxRevs    
    ContinueRevs = .TRUE.   
    FinalConv = .FALSE.             
    do while (ContinueRevs)  

        ! Increment revolution counter
        irev=irev+1                                                                                              

        ! Do timesteps
        do i=1,nti    

            ! Increment nt (total time step counter)
            nt=nt+1  

            ! Write diagnostic info to stdout if requested
            if (DiagOutFlag == 1) then
                write(6,'(A,I0)') 'Timestep: ', nt
            end if

            ! Set new wake element locations
            ! JCM: can move this function into expanded blade module when created...                                                                 
            CALL SetBoundWake()        

            ! Fixed-point iteration to converge non-linear system consisting of 
            ! blade element bound vorticity (potentially non linear with local AOA),
            ! and its own effect on local AOA. Wake remains constant during this iteration...
            iflg=0
            iter=0 
            ContinueNL=.TRUE.
            do while (ContinueNL) 

                ! Increment iterations                                                
                iter=iter+1 

                ! Initialize bound vorticity iteration:                                                      
                ! On first iteration, calc influence of all elements on blade AOA. On subsequent iterations of the bound vorticity, 
                ! only recalculate the bound vorticity component...
                if (iflg == 0) then 
                    CALL UpdateBladeVel(iflg) 
                    iflg=1
                else
                    CALL UpdateBladeVel(iflg) 
                end if

                If (Itower.EQ.1) Call UpdateTowerVelocity() !JCM: Matt to add this function...

                ! Regression test
                if (RegTFlag == 1) then
                    Reg_TS=nt
                    Reg_NLIter=iter
                end if

                ! Calculate blade loads and bound vorticity
                ! JCM: can move this function into expanded blade module when created...                                                           
                CALL BladeLoads(NLTol,iConv)                                                  

                if ((iConv == 0) .OR. (iter == MaxNLIters)) then   
                    ContinueNL=.FALSE.      
                end if

            end do

            if (iConv == 1) then                                      
                ! Non-linear iteration didn't converge. Write final output and exit.
                WRITE(6,610)                                                  
                CALL WriteFinalOutput()   
                stop                                                       
            end if

            ! Update strut loads
            CALL UpdateStrutLoads()

            ! Update influence on wall RHS and calc new wall panel strengths
            if (GPFlag == 1 .OR. WPFlag == 1 .or. FSFlag == 1) then
                Call UpdateWall() 
            end if

            ! Write output ----

            ! Collect timestep results and output
            Call EndTS()

            ! Write diagnostic info to stdout if requested
            if (DiagOutFlag == 1) then
                ! Use machine level time step output (norm. time, revolution, torque coeff., power coeff.)
                write(6,'(A)') 'Norm. Time, Theta (rad), Revolution, Torque Coeff., Power Coeff.'
                write(6,'(2E13.5,F8.0,2E13.5)') Output_TSData(1,1),Output_TSData(1,2),Output_TSData(1,3),Output_TSData(1,4),Output_TSData(1,5)
            end if

            ! Write current wake data
            if (WakeElementOutFlag > 0) then
                ! Write wake element data
                if ((NT >= WakeElementOutStartTimestep) .AND. (NT < WakeElementOutEndTimestep .OR. WakeElementOutEndTimestep == -1)) then
                    if (MOD(NT-1, WakeElementOutIntervalTimesteps) == 0) then
                        Call WriteWakeElementData()
                    end if
                end if
            end if

            if (WakeGridOutFlag > 0) then
                ! Write wake grid data
                if ((NT >= WakeGridOutStartTimestep) .AND. (NT < WakeGridOutEndTimestep .OR. WakeGridOutEndTimestep == -1)) then
                    if (MOD(NT-1, WakeGridOutIntervalTimesteps) == 0) then
                        Call WriteWakeGridData()
                    end if
                end if
            end if

            ! Write wall system/ground plane data
            if (WallOutFlag > 0) then
                ! Write wake grid data
                if ((NT >= WallOutStartTimestep) .AND. (NT < WallOutEndTimestep .OR. WallOutEndTimestep == -1)) then
                    if (MOD(NT-1, WallOutIntervalTimesteps) == 0) then
                        Call WriteWallData()
                    end if
                end if
            end if

            ! Write probe data
            if (ProbeFlag > 0) then
                ! Write wake grid data
                if ((NT >= ProbeOutStartTimestep) .AND. (NT < ProbeOutEndTimestep .OR. ProbeOutEndTimestep == -1)) then
                    if (MOD(NT-1, ProbeOutIntervalTimesteps) == 0) then
                        call write_probes() ! in probesystem module
                    end if
                end if
            end if

            ! State Updates ----  

            ! Update current wake convection velocities (excluding wake to be shed from the current blade)                                                                    
            CALL UpdateWakeVel()      ! Use all wake points to update the wake node velocities

            ! Convect the wake (including wake to be shed from the current blade location)
            CALL conlp()                                                                      

            ! Update bound vorticity timestep filter (used to reject grid scale temporal modes)
            if (TSFilFlag == 1) then
                Call UpdateTSFilter(ne)
            end if

            ! Shed new wake                         
            CALL shedvor()  

            ! Update states for the LB dynamic stall model, if used   
            if (DSFlag == 2) then                                 
                Call LB_UpdateStates(nb,nbe)   
            end if

            ! Update last AOA values
            Call UpdateAOALast(ne)                                   

            ! Rotate turbine geometry
            Call RotateTurbine

            ! Update time and phase
            TimeN=TimeN+dt
            Theta=Theta+delt

            ! Regression test
            if (RegTFlag == 1) then
                CALL WriteRegTOutput(2)
                if (nt == 2) then
                    stop
                end if
            end if

        end do ! Time steps 

        ! Calculate revolution average performance
        CPAve_last=CPAve                                                                                                                                
        CALL EndRev()                                          

        ! Write diagnostic info to stdout if requested
        if (DiagOutFlag == 1) then
            ! Write rev average power
            write(6,'(A,F13.2)') 'Revolution Average Power Coeff.: ', CPAve
            write(6,'(A,F13.2)') ' Rev Wall Time (sec): ', dtime
            write(6,'(A,F13.2)') ' Total Elapsed Wall Time (sec): ', etime
            write(6,'(A)') ' '
        end if

        ! If nr revs have been performed, then done. Otherwise, if initial convergence is hit, set final convergence params (if desired) and continue
        if (irev == nr) then
            ContinueRevs = .FALSE.
        else if (irev > 1) then

            ! Define convergence as convergence of revolution average power coeff. Note, negative values of convrg bypass the convergence check.
            ! Additionally, when using final convergence, the user can specify an intermediate revolution number at which to switch to
            ! final convergence (if initial convergence level hasn't already been hit, or is bypassed).
            if (abs((CPAve-CPAve_last)/CPAve) < convrg) then
                ConvFlag=.TRUE.
            else if (ifc == 1 .AND. .NOT. FinalConv .AND. irev == nric) then
                ConvFlag=.TRUE.
            else
                ConvFlag=.FALSE.
            end if

            if (ConvFlag) then
                ! If final convergence has already been performed (or not desired), then done. Otherwise, set parameters for final convergence and continue
                if (FinalConv .OR. (ifc == 0)) then
                    ContinueRevs = .FALSE.
                else
                    FinalConv = .TRUE.
                    ! Reset params for final convergence: Set higher (final) number of time steps, and associated wake update iteration period and 
                    ! more stringent convergence level (*f values) for final convergence                                                      
                    nti=ntif                                                                                  
                    iut=iutf 
                    ! If wake update interval set to a negative number, set next wake update iteration to -1 (no wake velocity updates will be performed).
                    ! Otherwise, restart wake updates on next iteration.
                    if (iut < 0) then
                        nsw=-1
                    else
                        nsw=nt+1        
                    end if
                    convrg=convrgf 
                    delt=2.0*pi/nti
                    DT=DelT/ut
                end if
            end if
        end if

    end do ! Revs                                                         

    ! Write output
    CALL WriteFinalOutput()    

610 FORMAT('0','***** NON-LINEAR ITERATION LOOP DID NOT CONVERGE IN 10 ITERATIONS. PROGRAM TERMINATED. *****')                                                    
End PROGRAM CACTUS
