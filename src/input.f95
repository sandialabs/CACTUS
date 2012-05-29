SUBROUTINE input(ErrFlag) 

	use parameters
	
	use dystl
	use element
	use blade
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
	
	integer, parameter :: InBufferNumSectionTables = 100
	integer, parameter :: InBufferNumSegPerBlade = 100
	integer, parameter :: InBufferNumSeg = 100
        integer, parameter :: MaxReadLine = 1000    
        integer, parameter :: MaxTempAOA = 1000   
        
	integer i, ii, jj, kk
	integer ErrFlag
	logical NotDone
        character(MaxReadLine) :: ReadLine
        integer :: CI, EOF
        real :: temp, temp1(MaxTempAOA,4)
               
	       
	! Temp buffers for airfoil section data inputs. If more sections or segments become necessary, its probably time to 
	! create a more generic geometry file with all the relevant geometry for one blade defined, and read this file in (using a loop).
	character(MaxReadLine) :: AFDPath(InBufferNumSectionTables)	! Airfoil section data path
	integer :: iSection(InBufferNumSegPerBlade)	! section index for each element
	real :: ChR(InBufferNumSegPerBlade)	! chord ratio for each element
	real :: bTwist(InBufferNumSegPerBlade)	! twist for each element (currently only used for HAWT geometry)
	integer :: WLI(InBufferNumSeg)     ! wake line index buffer
                
	! Namelist input file declaration
	NAMELIST/ConfigInputs/RegTFlag,DiagOutFlag,GeomFlag,GPFlag,rho,vis,tempr,hFSRef,slex,nr,convrg,nti,iut,ivtxcor,VCRFB,VCRFT,VCRFS,ifc,convrgf,nric,ntif,iutf,ixterm,xstop,Output_ELFlag,Incompr,DSFlag,PRFlag,k1pos,k1neg
	NAMELIST/XFlowInputs/jbtitle,Rmax,RPM,Ut,ChR,hr,eta,nb,nbe,nSect,AFDPath,iSection,hAG,Istraight,Istrut,sThick,Cdpar,CTExcrM,WakeOutFlag,WLI,BladeFileFlag
	NAMELIST/AxFlowInputs/jbtitle,R,HubR,RPM,Ut,Tilt,ChR,bCone,bi,bTwist,eta,nb,nbe,nSect,AFDPath,iSection,hAG,CTExcrM,WakeOutFlag,WLI
	
	! Input Defaults
        RegTFlag = 0 
        DiagOutFlag = 0
        Output_ELFlag = 0 
        WakeOutFlag = 0     
	nb = 2
	nbe = 5 
	nSect = 1  
	ifc = 0   
	nr = 10
	convrg = -1
	convrgf = .0001
	nti = 16
	ntif = 16
	ivtxcor = 0
	ixterm = 0
	xstop = 5.0
	iut = 0
	iutf = 0 
	nric=-1
        VCRFB=1.0       
	VCRFT=1.0
        VCRFS=1.0       
	iSection(:)=1 ! all set to 1 
	bTwist(:)=0.0 ! all set to 0
	ChR(:)=0.0 ! all set to 0
        WLI(:)=0 ! all set to 0      
	hAG=0.0
	Tilt=0.0
	bCone=0.0 
	HubR=0.0                                                    
	Incompr=0
	Istrut=0
	sThick=0.0                         
        Cdpar=0.0
        CTExcrM=0.0
        DSFlag=1
        PRFlag=1
        k1pos = 1.0                      
        k1neg = 0.5
        BladeFileFlag = 0
	                                                                              
	! Config Namelist input
	read(4, nml=ConfigInputs)                                                                              
	                                                                              
	! If GeomFlag is 1, import VAWT data, otherwise import HAWT data...                                           
	if (GeomFlag == 1) then
		! VAWT Namelist input
		read(4, nml=XFlowInputs)
	else
		! HAWT Namelist input
		read(4, nml=AxFlowInputs)
		Rmax=R
	end if
	
	! Set array bounds based on inputs
	! Geometry
	MaxBlades = nb
	MaxSegPerBlade = nbe
	MaxSegEndPerBlade = MaxSegPerBlade+1
	MaxSegEnds = MaxSegEndPerBlade*MaxBlades
        MaxSeg = MaxSegPerBlade*MaxBlades       
	! Airfoil Data
	MaxAirfoilSect = nSect
	MaxReVals = 20
	MaxAOAVals = 1000
	! Wake advancement
	MaxRevs = nr
	MaxTimeStepPerRev = nti
	MaxWakeNodes = MaxRevs * MaxTimeStepPerRev   
	! Non-linear convergence iteration
	MaxNLIters = 10
        ! Outputs
        MaxTimeSteps = MaxRevs * MaxTimeStepPerRev       
	! Wake outputs
        NWakeInd=0       
        NotDone=.TRUE.
        i=1
        do while (NotDone .AND. i<=MaxSeg)   
                if (WLI(i) > 0) then
                       NWakeInd=NWakeInd+1
                else
                       NotDone=.FALSE. 
                end if
                i=i+1
        end do
        
	! Array construction
        CALL blade_cns(MaxWakeNodes,MaxSegEnds)
	CALL element_cns(MaxTimeStepPerRev,MaxSegEnds,MaxSegEndPerBlade)     
	CALL airfoil_cns(MaxAOAVals,MaxReVals,MaxAirfoilSect)
	CALL wakedata_cns()
	CALL freestream_cns(MaxWakeNodes,MaxSegEnds)
	CALL dystl_cns(MaxAirfoilSect,MaxReVals,MaxSegEnds)
        CALL output_cns(MaxRevs, MaxTimeSteps, MaxSeg, MaxBlades)       
        
	! Write from buffer...
	iSect(1:(nbe+1))=iSection(1:(nbe+1))
        If (BladeFileFlag) Then
           Call ReadBladeFile(MaxSegEndPerBlade,yB,rr,cr)
        Else
           cr(1:(nbe+1))=ChR(1:(nbe+1))
        End If
	btw(1:(nbe+1))=bTwist(1:(nbe+1))     
        WakeLineInd(1:NWakeInd)=WLI(1:NWakeInd)       
	
	! Set hub radius ratio
	hubrr=HubR/Rmax
	
	! Set ground plane location for wall solution
	GPy=-hAG/Rmax
	
	! Normalize ground shear inputs
	yref = hFSRef/Rmax  ! location of freestream (99% maybe) normalized to radius... 
	ygc  = hAG/Rmax   ! Ground clearance normalized to radius  
	
	! JCM: Doing this real to int business floors the real tip speed ratio (ut) to the next lowest int
	! and uses this as the default wake velocity update interval... Probably OK...
	if (iutf == 0) iutf = ut  !? integer = real ?
	if (iut == 0) iut = ut  !? integer = real ?
	
	! Check chord to radius ratio. If a scalar has been input for cr, replicate for the entire blade.
	if (cr(2) == 0) then
		do i = 1,(nbe+1)
			cr(i)=cr(1)
		end do 
	end if
	                  				
	ne = (nbe+1)*nb ! Total number of blade segment ends (over all blades)
	
	! iSect setup
	if (nSect > 1) then                                       

		! Replicate input first blade iSect setup for all subsequent blades		
		nbse = nbe+1         
		do i = 2, nb            
			nei = 1+(i-1)*nbse
			iSect(nei) = iSect(1)           
			do j = 1, nbe   
				nej = nei+j        
				jp1 = j+1      
				iSect(nej) = iSect(jp1)      
			end do
		end do
		
		! Catch invalid iSect inputs and default them to use section 1, write warning to std out...
		do i = 1, ne
			if (iSect(i) < 1 .OR. iSect(i) > nsect) then 
				iSect(i) = 1     
				write (6,601) i
			end if 
		end do  
		
	else
		! All blades use first and only section
		iSect(:)=1
	end if

	
        ! Airfoil Data Tables: Read CL, CD, CM vs AOA from data files
        ! Format Example:
        ! Title: AFTitle
        ! Thickness to chord ratio: 0.2
        ! Zero Lift AOA (deg): 0.0
        ! Reverse camber direction: 0
        !
        ! Reynolds Number: 1e6
        ! BV Dyn. Stall Model - Positive stall AOA (deg): 10
        ! BV Dyn. Stall Model - Negative stall AOA (deg): -10
        ! LB Dyn. Stall Model - Lift Coeff. Slope at Zero Lift AOA (per radian): 6.28
        ! LB Dyn. Stall Model - Positive Critical Lift Coeff.: 1.3
        ! LB Dyn. Stall Model - Negative Critical Lift Coeff.: -1.3
        ! AOA (deg) CL CD Cm
        ! ... ... ... ...
        ! 
        ! Reynolds Number: 5e6
        ! ...                 
        do kk = 1, nsect
        
                ! Open input file for this section
                open(15, file=AFDPath(kk))
                EOF=0
                
                ! Find title block 
                NotDone=.TRUE.
                do while (NotDone)
                    read(15,'(A)') ReadLine
                    CI=index(ReadLine,':')
                    if (CI>0) then
                        NotDone=.FALSE.
                    end if
                end do
                
                ! Read title and airfoil thickness
                if (len_trim(ReadLine)>CI) then
                    aftitle(kk) = ReadLine(CI+1:len_trim(ReadLine))
                else
                    aftitle(kk) = 'No Title'
                end if
                read(15,'(A)') ReadLine
                read(ReadLine(index(ReadLine,':')+1:),*) tc(kk) 
                read(15,'(A)') ReadLine
                read(ReadLine(index(ReadLine,':')+1:),*) alzer(kk)
                alzer(kk)=alzer(kk)*conrad
                read(15,'(A)') ReadLine
                read(ReadLine(index(ReadLine,':')+1:),*) camb(kk)                            
                
                ! Reverse camber direction if desired                                  
                if (camb(kk) == 1) then
                    alzer(kk) = -alzer(kk)
                end if                                  
                                                  
                ! Find first Re block
                NotDone=.TRUE.
                do while (NotDone)
                    read(15,'(A)',IOSTAT=EOF) ReadLine
                    CI=index(ReadLine,':')
                    if (CI>0 .OR. EOF<0) then
                        NotDone=.FALSE.
                    end if
                end do                                                  
                                                  
                ! Read data for each Re value 
                i=0                                                                
                do while (EOF >= 0  .AND. (i < MaxReVals)) 
                    
                    i=i+1
                    ! Read Re and dyn. stall data                                
                    read(ReadLine(index(ReadLine,':')+1:),*) TRE(i,kk)
                    read(15,'(A)') ReadLine                          
                    read(ReadLine(index(ReadLine,':')+1:),*) alstlp(i,kk)
                    alstlp(i,kk)=alstlp(i,kk)*conrad
                    read(15,'(A)') ReadLine                          
                    read(ReadLine(index(ReadLine,':')+1:),*) alstln(i,kk)
                    alstln(i,kk)=alstln(i,kk)*conrad
                    read(15,'(A)') ReadLine                          
                    read(ReadLine(index(ReadLine,':')+1:),*) CLaData(i,kk)
                    read(15,'(A)') ReadLine                          
                    read(ReadLine(index(ReadLine,':')+1:),*) CLCritPData(i,kk)
                    read(15,'(A)') ReadLine                          
                    read(ReadLine(index(ReadLine,':')+1:),*) CLCritNData(i,kk)
                    
                    ! Reverse camber direction if desired
                    if (camb(kk) == 1) then
                        temp = alstlp(i,kk)
                        alstlp(i,kk) = -alstln(i,kk)
                        alstln(i,kk) = -temp   
                        temp = CLCritPData(i,kk)
                        CLCritPData(i,kk) = -CLCritNData(i,kk)
                        CLCritNData(i,kk) = -temp 
                    end if
                    
                    ! Read AOA data
                    read(15,'(A)') ReadLine
                    NotDone=.TRUE.
                    ii=0
                    do while (NotDone)
                        read(15,'(A)',IOSTAT=EOF) ReadLine
                        if (EOF>=0 .AND. len_trim(ReadLine)>0) then
                            if (ii == MaxAOAVals) then
                                write(6,*) 'Max. allowed AOA values exceeded in airfoil data file: ', aftitle(kk)
                                ErrFlag=1
                                NotDone=.FALSE.
                            else
                                ii=ii+1                        
                                read(ReadLine,*) ta(ii,i,kk),tcl(ii,i,kk),tcd(ii,i,kk),tcm(ii,i,kk) 
                            end if
                        else
                            NotDone=.FALSE.
                        end if
                    end do
                    ntbl(i,kk)=ii
                    
                    ! Check AOA limits
                    if (ta(1,i,kk) > -180.0 .OR. ta(ntbl(i,kk),i,kk) < 180.0) then
                        write(6,*) 'AOA data needs to be +/-180 deg in airfoil data file: ', aftitle(kk)
                        ErrFlag=1
                    end if
                    
                    ! Reverse camber direction if desired
                    if(camb(kk) == 1) then       
                        do ii = 1, ntbl(i,kk)
                                temp1(ii,1) = ta(ii,i,kk) 
                                temp1(ii,2) = tcl(ii,i,kk)
                                temp1(ii,3) = tcd(ii,i,kk)
                                temp1(ii,4) = tcm(ii,i,kk)
                        end do  
                        
                        do ii = 1, ntbl(i,kk)
                                jj = ntbl(i,kk)-(ii-1)
                                ta(ii,i,kk) = -temp1(jj,1)
                                tcl(ii,i,kk) = -temp1(jj,2)
                                tcd(ii,i,kk) = temp1(jj,3)
                                tcm(ii,i,kk) = -temp1(jj,4)
                        end do  
                    end if
                    
                    ! Find next Re block
                    NotDone=.TRUE.
                    if (EOF<0) then
                        NotDone=.FALSE.
                    end if
                    do while (NotDone)
                        read(15,'(A)',IOSTAT=EOF) ReadLine
                        CI=index(ReadLine,':')
                        if (CI>0 .OR. EOF<0) then
                            NotDone=.FALSE.
                        end if
                    end do 
                                        
                end do
                ! Set number of Re vals for this section
                nRET(kk)=i

                ! Close input file for this section
                close(15)
                
                ! Check data
                if (i == 0) then
                    write(6,*) 'Error reading airfoil data file: ', aftitle(kk)
                    ErrFlag=1
                end if
                if (EOF > 0) then
                    write(6,*) 'Warning: Max. allowed Re values exceeded in airfoil data file: ', aftitle(kk)
                end if
                      
        end do
									
Return  									                                        
601  format(' ','***airfoil section specified for blade segment ',i2,' is illegal. set to airfoil section 1***')                        
End 
	

! Subroutine: ReadBladeFile
! Reads in a blade geometry file containing radius and chord as a function of vertical position						
Subroutine ReadBladeFile(Nsegends,y,r,chord)

  Use ioption

  Implicit None

  Integer :: Nsegends, N, I
  Real, Dimension(Nsegends) :: y, r, chord
  Integer, Parameter :: inunit=10

  Open(unit=inunit,file=BladeFileName,form='formatted',status='unknown')
  Read(inunit,*) N
  If (N .NE. Nsegends) Then
     Write(6,'(''Number of blade coordinates in blade file does not match the input file.'')')
     STOP
  End If
  Read(inunit,*) (y(I), r(I), chord(I), I=1,N)
  Close(inunit)

  Write(*,*) (y(I), r(I), chord(I), I=1,N)
  Return

End Subroutine ReadBladeFile


