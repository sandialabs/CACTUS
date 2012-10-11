PROGRAM CACTUS

        ! Development code for wind/water turbine performance calculation based on VDART3.
        !       J. Murray, 3/2010
        ! -----------
        ! Notes: 
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
        use shear
        use airfoil
        use configr
        use pidef
        use ioption
        use vortex
        use wakedata
        use time
        use freestream
        use wallsoln
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
        real :: NLTol
        real :: cpave, cpave_last
        real :: delt, delty, deltb, deltr
                                                      
        character(80) :: InputFN, SFOutputFN, RevOutputFN, TSOutputFN, ELOutputFN, RegOutputFN, WakeOutputFN, WakeDefOutputFN, GPOutputFN, FSOutputFN, FNBase                                   
                                                 
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
                                                                            
        ! Optional wake line data output
        if (WakeOutFlag > 0) then
                WakeOutputFN=trim(FNBase)//'_WakeData.csv'
                OPEN(12, FILE=WakeOutputFN)
                
                if (WakeOutFlag > 1) then
                        ! wake deficit surface output
                        WakeDefOutputFN=trim(FNBase)//'_WakeDefData.csv'
                        OPEN(13, FILE=WakeDefOutputFN)
                end if
        end if                                                                    
                                                                            
        ! Optional wall model output
        if (WallOutFlag > 0) then
                if (GPFlag == 1) then
                    GPOutputFN=trim(FNBase)//'_GPData.csv'
                    OPEN(14, FILE=GPOutputFN)
                end if
                
                if (FSFlag == 1) then
                    FSOutputFN=trim(FNBase)//'_FSData.csv'
                    OPEN(15, FILE=FSOutputFN)
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
        if (GPFlag == 1 .OR. FSFlag == 1) then
                nsWall=1
        end if
        
        ! Blade Geometry setup. 
        ! Global axes: x is oriented with the nominal freestream flow direction, y is in the vertically upward direction (opposite gravity),
        ! z direction from RHR (to the right when looking in the streamwise direction).
        CALL BGeomSetup()
       
        ! Set normalized turbine rotation rate
        wRotX=ut*RotX
        wRotY=ut*RotY  
        wRotZ=ut*RotZ
        delt=2.0*pi/nti 
        
        ! Setup wall geometry and solution if necessary
        if (GPFlag == 1 .OR. FSFlag == 1) then
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

        ! Setup vortex core radius for bound vorticity based on max chord, and for trailing and spanwise wake based on 
        ! temporal and spatial discretization levels, respectively.
        if (ivtxcor == 0) then
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
        Call csvwrite(8,Output_SFHead,Output_SFData,0)
                                                                                                                 
        ! Initialize needed arrays                             
        do i=1,ne                                                      
                gs(1,i)=0.0                                                       
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
        Time1=t0                                                                             
           
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
                                write(6,*) 'Timestep: ', nt
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
                            
                        ! Write diagnostic info to stdout if requested
                        if (DiagOutFlag == 1) then
                                ! Use machine level time step output (norm. time, revolution, torque coeff., power coeff.)
                                write(6,*) 'Norm. Time, Revolution, Torque Coeff., Power Coeff.'
                                write(6,'(E13.5,F8.0,2E13.5)') Output_TSData(nt,1),Output_TSData(nt,2),Output_TSData(nt,3),Output_TSData(nt,4)
                        end if    
                            
                        ! Update influence on wall RHS and calc new wall panel strengths
                        if (GPFlag == 1 .OR. FSFlag == 1) then
                                Call UpdateWall() 
                        end if 
  
                        ! Write current wake data for viewing in Matlab
                        if (WakeOutFlag > 0) then
                                Call WriteWakeData()
                        end if  
  
                        ! Write current wall data for viewing in Matlab
                        if (WallOutFlag > 0) then
                                Call WriteWallData()
                        end if    
  
                        ! Update current wake convection velocities (excluding wake to be shed from the current blade)                                                                    
                        CALL UpdateWakeVel()      ! Use all wake points to update the wake node velocities

                        ! State Updates ----

                        ! Convect the wake (including wake to be shed from the current blade location)
                        CALL conlp()                                                                      
                        
                        ! Shed new wake                         
                        CALL shedvor()  

                        ! Update states for the LB dynamic stall model, if used   
                        if (DSFlag == 2) then                                 
                                Call LB_UpdateStates(nb,nbe)   
                        end if    

                        ! Update last AOA values
                        Call UpdateAOALast(ne)                                   

                        ! Rotate turbine geometry
                        Call RotateTurbine(delt)

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

                ! Write diagnostic info to stdout if requested
                if (DiagOutFlag == 1) then
                        ! Write rev average power
                        write(6,*) 'Revolution Average Power Coeff.: ', cpave
                        write(6,*) ' '
                end if                                           
                                                           
                ! If nr revs have been performed, then done. Otherwise, if initial convergence is hit, set final convergence params (if desired) and continue
                if (irev == nr) then
                        ContinueRevs = .FALSE.
                else if (irev > 1 .AND. convrg > 0.0) then
                        
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
                                        ! Reset CPF and KPF fit in endrev (output)
                                        FitStartRev=irev
                                end if
                        end if
                end if
                
        end do ! Revs                                                         
                                                                                                                                                     
        ! Write output
        CALL WriteFinalOutput()    
                                                        
610  FORMAT('0','***** NON-LINEAR ITERATION LOOP DID NOT CONVERGE IN 10 ITERATIONS. PROGRAM TERMINATED. *****')                                                    
End                                                              
