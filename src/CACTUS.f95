PROGRAM CACTUS

	! Development code for wind/water turbine performance calculation based on VDART3.
	! 	J. Murray, 3/2010
	!	J. Murray, 3/2010: Name changed from PoTHyDyn (Potential flow Turbine HydroDynamics) to CACTUS (Code for Axial and Crossflow Turbine Simulation) 
	!		as the acronym was deemed less comical...
	! -----------
	! Notes: 
	!	J. Murray, 3/10/10: Uses f2kcli to add F2003-like command line interface.
	! 		F2KCLI: Fortran 200x Command Line Interface
	! 			copyright Interactive Software Services Ltd. 2001
	!			http://www.winteracter.com/f2kcli
	!		Command line inputs:
	!			FILENAME: Full name of the input namelist file, if located in the calling directory. Full path otherwise.			
	
	! Note that all internal variables are normalized. Velocity variables are normalized by freestream velocity, distance variables
	! are normalized by reference radius, velocity derivatives by (freestream velocity)/(reference radius), circulation variables 
	! by (freestream velocity)*(reference radius)
	
	use parameters
	use dystl
	use wakeloc
	use element
	use vel
	use veo
	use gam
	use rvrr
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
	
	!IMPLICIT NONE !JCM: eventually...	
	
	integer :: ErrFlag, nargin, FNLength, status, DecInd
	logical :: FinalConv
	logical :: ConvFlag
	logical :: ContinueRevs
	logical :: ContinueNL
	logical :: back
	integer :: iConv
	integer :: WakeLineInd(4)
	real :: cpsum, NLTol
	real :: delt, delty, deltb, deltr
	
	real :: hmachl(2)
	real :: hmachm(2)
                                                 
        character(80) :: InputFN, OutputFN, DOutputFN, SOutputFN                                    
                                                 
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
	DecInd=index(InputFN,'.',back)
	if (DecInd > 1) then
		OutputFN=InputFN(1:(DecInd-1))//'.out'
		DOutputFN=InputFN(1:(DecInd-1))//'_Data.out'
		SOutputFN=InputFN(1:(DecInd-1))//'_Short.out'
	else
		OutputFN=trim(InputFN)//'.out'
		DOutputFN=trim(InputFN)//'_Data.out'
		SOutputFN=trim(InputFN)//'_Short.out'
	end if
	
	! Namelist input file                                                       
	OPEN(4, FILE= InputFN)                                     
	
	! Output files                                                      
	OPEN(6, FILE= OutputFN,  FORM= 'FORMATTED' )                                    
	OPEN(9, FILE= DOutputFN,   FORM= 'FORMATTED' )                                         
	OPEN(12, FILE= SOutputFN,  FORM= 'FORMATTED' )   
	
	! Read the Date and Time.
	CALL dattim(DMY,HMS)
	
	! Initialize iteration parameters							
	irev=0
	nt = 0
	ntTerm=1                                                            	                                                                                                         
	cpsum = 0.0   
	power = 0.0                                                    
	FitStartRev=1  
	NLTol=1.0e-04                                                      
	
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
		CALL BGeomSetup_v(nti,nbe,nb,delty,delt,deltb,at)
		
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
		CALL BGeomSetup_h(nti,nbe,nb,deltr,hubrr,delt,deltb,at) 
		
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
	uMPH=uinf*3600./5280.                                              
	rem=rho*uinf*CrRef*Rmax/vis                                          
	Minf=uinf/sqrt(1.4*1716.0*(tempr+459.6))                                                                         							          
	areat=at*Rmax**2     				! frontal area (at is (frontal area) / Rmax^2 )                                            
	trqcon=rho/2.0*areat*Rmax*uinf**2       	! torque coeff normalization                          
	powerc=rho/2.0*areat*romega**3*0.7457/550.  	! normalization for power coeff using tip speed (kp), with conversion from lb-ft/s to kW. (Used to write output)
	
	! Dynamic stall setup                                                
	k1pos = 1.0                                                                                                                                                         
	k1neg = 0.5                                          
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

	! Setup temp arrays				
	do i=1,ne                                                      
		gs(1,i)=0.0                                                       
		ogb(i)=0.0                                                        
		alfold(i)=0.0                                                                                                         
	end do 

	! List input data
	CALL listit(0,6)
	! Write header on data output file
	write (9,901) jbtitle      
	write (9,902) trqcon
	
	! CPU time markers
	t0 = secnds(0.0)                                                                                                           
	Time1 = secnds(t0)                                                                             
	   
	! Do revolutions until convergence or MaxRevs    
	ContinueRevs = .TRUE.	
	FinalConv = .FALSE.		
	do while (ContinueRevs)  
		
		! Increment revolution counter
		irev=irev+1
		
		! Reset power coeff. sum
		cpsum=0.0	                                                                                                    
	
		! Do timesteps
		do i=1,nti    
		         
			! Increment nt (total time step counter)
			nt=nt+1  
			
			! Get current geometry                                                                 
			CALL bgeom(i,nt,nbe,nb)        
			
			iflg=0                                         
			iter=0 
			! Calculate wake and wall induced velocity at the current blade locations                                                    
			CALL bivel(nt,ntTerm,nbe,nb,ne,iflg)    
                                             
			! Fixed-point iteration to converge non-linear system consisting of 
			! blade element bound vorticity (potentially non linear with local AOA),
			! and its own effect on local AOA. Wake remains constant during this iteration...
			ContinueNL=.TRUE.
			do while (ContinueNL) 
			        
				! Increment iterations                                                
				iter=iter+1                                                       
				
				! Calculate blade loads and bound vorticity                                                           
				CALL bvort(i,NLTol,iConv)                                                  
	
				iflg=1 
				! Calculate wake induced velocity at the current blade locations   
				CALL bivel(nt,ntTerm,nbe,nb,ne,iflg)                      
			
				if ((iConv == 0) .OR. (iter == MaxNLIters)) then   
					ContinueNL=.FALSE.
				end if
			
			end do

			if (iConv == 1) then                                      
				! Non-linear iteration didn't converge. Write final output and exit.
				WRITE(6,610)  
				                                                    
				CALL perf(i,cpl) !JCM Test                                                    
				                                                    
				CALL WriteFinalOutput()   
				stop                                                       
                        end if                                 
			
			! Calculate current performance parameters
			! JCM: note that in addition to calculating performance, perf actually updates the bound vorticity strength once more, following the 
			! last update of the blade velocities in bivel from the loop above... This should probably be reorganized somehow to avoid confusion... 
			! Possibly, could include in bvort with an output flag to be set on the last iteration...
			CALL perf(i,cpl)                                          
			! Update power coeff. sum
			cpsum=cpsum+cpl
			     
			! Update freestream, bound and wake vorticity influence on wall RHS and calc new wall panel strengths
			if (GPFlag == 1) then
				Call UpdateWall(nt,ntTerm,nbe,nb,ne,iut,nsw) 
			end if 
  
  			! JCM: Write current wake data for viewing in Matlab
			WakeOut=0
			if (WakeOut == 1) then
				Call WriteWakeData(i,DelT,WakeLineInd)
			end if	
  
			! Update current wake convection velocities (excluding wake to be shed from the current blade)                                                                    
			nfpw=ifw*jfw*kfw  ! Number of fixed wake points                                                
			if ((npw <= nfpw) .OR. (ifwg == 0)) then                     
                        	CALL wivel(nt,ntTerm,nbe,nb,ne,iut,nsw,npw)      ! Use all wake points to update the wake node velocities
			else                        
				CALL swivel(nt,ntTerm,nbe,nb,ne,iut,nsw,npw)     ! Use fixed wake grid to update the wake node velocities by interpolation                      
			end if  
			                                                       
			! If new wall and wake velocities were calculated on this timestep, set the next update timestep
			if (nt .eq. nsw) then                                         
				nsw=nt+iut 
      			end if                                                       
			                                                       
			! Convect the wake (including wake to be shed from the current blade location)
			CALL conlp(nt,ntTerm,ne,delt,ut)                                                                      
			
			! Shed new wake		    		
			CALL shedvor(nt,nbe,nb)                                           
						
		end do ! Time steps 
		
		! Calculate revolution average performance					                                                                                          
		CALL endrev(cpsum)                                          
                                                           
		! If nr revs have been performed, then done. Otherwise, if initial convergence is hit, set final convergence params (if desired) and continue
		if (irev == nr) then
			ContinueRevs = .FALSE.
		else if (irev > 1) then
			
			! Define convergence as convergence of revolution average power coeff. Additionally, when using final convergence, the user can
			! specify an intermediate revolution number at which to switch to final convergence (if initial convergence level hasn't already been hit).
			if (abs((cp(irev)-cp(irev-1))/cp(irev)) < convrg) then
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
					! Recreate the geometry arrays with the new theta resolution (Note: VAWT/HAWT specific since we are recreating geometry)
					if (GeomFlag == 1) then	                                                      
						CALL BGeomSetup_v(nti,nbe,nb,delty,delt,deltb,at)
					else
						CALL BGeomSetup_h(nti,nbe,nb,deltr,hubrr,delt,deltb,at) 
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
901  FORMAT('RESULTS FOR JOB: ',A80)               
902  FORMAT('TORQUE CONVERSION FACTOR=',F15.0)                       
End                                                              
