PROGRAM CACTUS

        ! Development code for wind/water turbine performance calculation based on VDART3.
        !       J. Murray, 3/2010
        !       J. Murray, 3/2010: Name changed from PoTHyDyn (Potential flow Turbine HydroDynamics) to CACTUS (Code for Axial and Crossflow Turbine Simulation) 
        !               as the acronym was deemed less comical...
        ! -----------
        ! Notes: 
        !       J. Murray, 3/10/10: Uses f2kcli to add F2003-like command line interface.
        !               F2KCLI: Fortran 200x Command Line Interface
        !                       copyright Interactive Software Services Ltd. 2001
        !                       http://www.winteracter.com/f2kcli
        !               Command line inputs:
        !                       FILENAME: Full name of the input namelist file, if located in the calling directory. Full path otherwise.                       
        
        ! Note that all internal variables are normalized. Velocity variables are normalized by freestream velocity, distance variables
        ! are normalized by reference radius, velocity derivatives by (freestream velocity)/(reference radius), circulation variables 
        ! by (freestream velocity)*(reference radius)
        
        use parameters
        use dystl
        use blade
        use element
        use varscale
        use cltab
        use shear
        use airfoil
        use configr
        use xwake
        use pidef
        use test
        use ioption
        use vortex
        use uwake
        use wakedata
        use time
        use freestream
        use wallsoln
        use f2kcli
        use regtest
        use output       
        
        !IMPLICIT NONE !JCM: eventually...      
    
        include 'csvwrite.inc'
        
        integer :: ErrFlag, nargin, FNLength, status, DecInd
        logical :: FinalConv
        logical :: ConvFlag
        logical :: ContinueRevs
        logical :: ContinueNL
        logical :: back
        integer :: iConv
        integer :: WakeLineInd(4)
        real :: NLTol
        real :: cpave, cpave_last
        real :: delt, delty, deltb, deltr
        
        real :: hmachl(2)
        real :: hmachm(2)
                                                 
        character(80) :: InputFN, SFOutputFN, RevOutputFN, TSOutputFN, ELOutputFN, RegOutputFN, FNBase                                   
                                                 
        ! Pi definition
        pi = 4.0*atan(1.0)
        conrad = pi/180.0                                                  
        condeg = 180.0/pi 
        
        ! Parse command line for the name of the input file
        nargin=command_argument_count()
        if (nargin < 1) then
                write(6,*) 'Please call the program with the name of the input file on the command line. Ex. CACTUS INPUTFILE.in' 
                stop            
        end if
        Call get_command_argument(1,InputFN,FNLength,status)
        back=.true.
        FNBase=InputFN((index(InputFN,'/',back)+1):len(InputFN))
        DecInd=index(FNBase,'.',back)     
        if (DecInd > 1) then
                FNBase=FNBase(1:(DecInd-1))       
        end if       
        SFOutputFN=trim(FNBase)//'_Param.csv'
        RevOutputFN=trim(FNBase)//'_RevData.csv'
        TSOutputFN=trim(FNBase)//'_TimeData.csv'
        ELOutputFN=trim(FNBase)//'_ElementData.csv'        
        
        ! Namelist input file                                                       
        OPEN(4, FILE= InputFN)                                     
        
        ! Output files                                                      
        OPEN(8, FILE=SFOutputFN) 
        OPEN(9, FILE=RevOutputFN)
        OPEN(10, FILE=TSOutputFN)

        ! Initialize iteration parameters                                                       
        irev=0
        nt = 0
        ntTerm=1                                                                                                                                                       
        FitStartRev=1  
        NLTol=1.0e-04                                                      
        cpave=0.0
        cpave_last=0.0
        
        ! Error flags                                                                                                                                                                  
        ilxtp=0
        iuxtp=0
        ierr0=0                                                           
        ierr1=0
        
        ! Read inputs    
        ErrFlag = 0
        CALL input(ErrFlag)
        if (ErrFlag == 1) then
                write(6,*) 'Input data fail. Exiting...'
                stop
        end if
        
        ! Simple output for regression testing        
        if (RegTFlag == 1) then      
                 RegOutputFN=trim(FNBase)//'_RegData.out'
                 OPEN(7, FILE= RegOutputFN,  FORM= 'FORMATTED' )  
                 Call WriteRegTOutput(0)
        end if
        
        ! Optional element load output
        if (Output_ELFlag == 1) then
                OPEN(11, FILE=ELOutputFN)
        end if
                                                                            
        ! If wake update interval set to a negative number, set next wake update iteration to -1 (no wake velocity updates will be performed)
        ! Otherwise, make the first update on the second iteration (when the wake first appears)
        if (iut < 0) then
                nsw=-1
        else
                nsw=2
        end if 
        
        ! Initialize fixed wake grid data with arbitrary x station intervals, y and z calculated with constant intervals...
        npw=0 ! Number of points in the fixed wake is used to determine whether to use the fixed wake calculation, initialized to zero and incremented in swivel as necessary
        ifw=1 ! Initial value of ifw set to 1. Will be modified in swivel when free wake x exceeds xfw(ifw)
        xfw=[1.01,1.5,2.0,3.0,5.0,7.0,10.0,15.0,20.0,25.0,35.0,50.0]                        
        jfw=5 ! Desired number of fixed wake points in the y direction (constant)
        kfw=5 ! Desired number of fixed wake points in the z direction (constant)
        
        if (GeomFlag == 1) then
                yfw(1)=-0.5 
                yfw(jfw)=hr+0.5
        else
                yfw(1)=-1.5 
                yfw(jfw)=1.5
        end if
        jfwm1=(jfw-1) 
        dyfw=(yfw(jfw)-yfw(1))/jfwm1
        do i = 2, jfwm1    
                yfw(i) = yfw(i-1)+dyfw    
        end do    
        
        zfw(1)=-1.5 
        zfw(kfw)=1.5
        kfwm1=(kfw-1) 
        dzfw=(zfw(kfw)-zfw(1))/kfwm1
        do i = 2, kfwm1    
                zfw(i) = zfw(i-1)+dzfw    
        end do         
        
        ! Blade Geometry setup  
        delt=2.0*pi/nti 
        deltb=2.0*pi/nb
        if (GeomFlag == 1) then 
                ! Create VAWT geometry                                                                                                                              
                delty=hr/nbe                                                 
                CALL BGeomSetup_v(delty,delt,deltb,at)
                
                ! VAWT axis and normalized rotation rate
                RotX=0.0
                RotY=1.0
                RotZ=0.0
                wRotX=0.0
                wRotY=ut  
                wRotZ=0.0
                
                ! Wakelines to output for VAWT (in WriteWakeData)
                WakeLineInd=[2,4,8,10]  
        else
                ! Create HAWT geometry                                                                                                                                   
                deltr=(1.0-hubrr)/nbe                                                       
                CALL BGeomSetup_h(deltr,hubrr,delt,deltb,at) 
                
                ! HAWT axis and normalized rotation rate
                RotX=1.0
                RotY=0.0
                RotZ=0.0
                wRotX=ut
                wRotY=0.0  
                wRotZ=0.0
                
                ! Wakelines to output for HAWT (in WriteWakeData)
                WakeLineInd=[3,6,9,12]
        end if
        
!---------- VAWT/HAWT specific geometry creation code above here -----------
        
        ! Setup wall geometry and solution if necessary
        if (GPFlag == 1) then
                ! Wall Geometry setup
                Call WGeomSetup()
                
                ! Setup wall solution
                Call WSolnSetup()
        end if

        ! Normalization parameters for geometry and performance outputs
        romega=2.0*pi*Rmax*rpm/60.0                                       
        uinf=romega/ut 
        DT=DelT/ut                                      ! normalized simulation timestep (dt*Uinf/Rmax)                                               
        rem=rho*uinf*Rmax/vis                                          
	
	if (Incompr .eq. 1) then			! incompressible/compressible flow switch (used by dynamic stall model)
		Minf = 0.0
	else
	        Minf=uinf/sqrt(1.4*1716.0*(tempr+459.6))
	end if

        areat=at*Rmax**2                                ! frontal area (at is (frontal area) / Rmax^2 ) 
        dynpress=rho/2.0*uinf**2                        ! dynamic pressure                  
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
        Output_SFData(1,10)=rem          ! machine Reynolds number based on U and Rmax
        Call csvwrite(8,Output_SFHead,Output_SFData,0)
        
        ! Dynamic stall setup                                                
        k1pos = 1.0   ! effect magnitude for CL increasing                      
        k1neg = 0.5   ! effect magnitude for CL decreasing                                      
        do i=1,nsect                                                   
                diff=0.06-tc(i)                                                   
                smachl(i)=0.4+5.0*diff                                            
                hmachl(i)=0.9+2.5*diff                                            
                gammaxl(i)=1.4-6.0*diff                                           
                dgammal(i)=gammaxl(i)/(hmachl(i)-smachl(i))                       
                smachm(i)=0.2                                                     
                hmachm(i)=0.7+2.5*diff                                            
                gammaxm(i)=1.0-2.5*diff                                           
                dgammam(i)=gammaxm(i)/(hmachm(i)-smachm(i))                       
        end do 
                                                                                                                   
        ! Initialize needed arrays                             
        do i=1,ne                                                      
                gs(1,i)=0.0                                                       
                ogb(i)=0.0 
                AOA(i)=0.0                                                                                                                                                             
        end do 

        ! CPU time markers
        t0 = secnds(0.0)                                                                                                           
        Time1 = secnds(t0)                                                                             
           
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
                        
                        ! Get current geometry                                                                 
                        CALL bgeom(i)        
                     
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
                                ! Set dynamic stall alpha old to current blade AOA (old AOA left fixed during bound vorticity iteration)
                                if (iflg == 0) then 
                                        do k=1,ne                                                                                                             
                                                alfold(k)=AOA(k)                                                         
                                        end do  
                                 
                                        CALL bivel(iflg) 
                                        iflg=1
                                else
                                        CALL bivel(iflg) 
                                end if
                                
                                ! Regression test
                                if (RegTFlag == 1) then
                                        Reg_TS=nt
                                        Reg_NLIter=iter
                                end if
                                
                                ! Calculate blade loads and bound vorticity                                                           
                                CALL bvort(i,NLTol,iConv)                                                  
        
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
                            
                        ! Update freestream, bound and wake vorticity influence on wall RHS and calc new wall panel strengths
                        if (GPFlag == 1) then
                                Call UpdateWall() 
                        end if 
  
                        ! JCM: Write current wake data for viewing in Matlab
                        ! JCM: Should change this to write all wake lines in a csv. 
                        ! For each timestep: first col is blade, second is elem, then location, velocity. Each row is a wake line element...
                        ! Similar function for wall elements...
                        ! Make this output an option in the input file...
                        WakeOut=0
                        if (WakeOut == 1) then
                                Call WriteWakeData(i,DelT,WakeLineInd)
                        end if  
  
                        ! Update current wake convection velocities (excluding wake to be shed from the current blade)                                                                    
                        nfpw=ifw*jfw*kfw  ! Number of fixed wake points                                                
                        if ((npw <= nfpw) .OR. (ifwg == 0)) then                     
                                CALL wivel()      ! Use all wake points to update the wake node velocities
                        else                        
                                CALL swivel()     ! Use fixed wake grid to update the wake node velocities by interpolation                      
                        end if  
                                                                               
                        ! If new wall and wake velocities were calculated on this timestep, set the next update timestep
                        if (nt .eq. nsw) then                                         
                                nsw=nt+iut 
                        end if                                                       
                                                                               
                        ! Convect the wake (including wake to be shed from the current blade location)
                        CALL conlp()                                                                      
                        
                        ! Shed new wake                         
                        CALL shedvor()  
                                                                 
                        ! Regression test
                        if (RegTFlag == 1) then
                                CALL WriteRegTOutput(2)
                                if (nt == 2) then
                                        stop
                                end if                        
                        end if
                                                
                end do ! Time steps 
                
                ! Calculate revolution average performance
                cpave_last=cpave                                                                                                                                
                CALL endrev(cpave)                                          
                                                           
                ! If nr revs have been performed, then done. Otherwise, if initial convergence is hit, set final convergence params (if desired) and continue
                if (irev == nr) then
                        ContinueRevs = .FALSE.
                else if (irev > 1) then
                        
                        ! Define convergence as convergence of revolution average power coeff. Additionally, when using final convergence, the user can
                        ! specify an intermediate revolution number at which to switch to final convergence (if initial convergence level hasn't already been hit).
                        if (abs((cpave-cpave_last)/cpave) < convrg) then
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
                                        ! Recreate the geometry arrays with the new theta resolution (Note: VAWT/HAWT specific since we are recreating geometry)
                                        if (GeomFlag == 1) then                                                       
                                                CALL BGeomSetup_v(delty,delt,deltb,at)
                                        else
                                                CALL BGeomSetup_h(deltr,hubrr,delt,deltb,at) 
                                        end if
                                        ! Reset CPF and KPF fit in endrev (output)
                                        FitStartRev=irev
                                end if
                        end if
                end if
                
        end do ! Revs                                                         
                                                                                                                                                     
        ! Write output
        CALL WriteFinalOutput()    
                                                        
610  FORMAT('0','***** BIVEL/BVORT LOOP DID NOT CONVERGE IN 10 ITERATIONS. PROGRAM TERMINATED. *****')                                                       
End                                                              
