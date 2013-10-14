SUBROUTINE input(ErrFlag) 

	use parameters

	use dystl
	use element
	use blade
        use strut       
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
    

	integer, parameter :: InBufferNumSectionTables = 100
	integer, parameter :: InBufferNumWL = 100
    integer, parameter :: MaxReadLine = 1000    
    integer, parameter :: MaxTempAOA = 1000   

	integer i, ii, jj, kk
	integer ErrFlag
	logical NotDone, NotBlank
    character(MaxReadLine) :: ReadLine
    character(MaxReadLine) :: GeomFilePath    ! path to geometry input file 
    integer :: CI, EOF
    real :: temp, temp1(MaxTempAOA,4)
        

	! Temp buffers
	character(MaxReadLine) :: AFDPath(InBufferNumSectionTables)	! Airfoil section data path      
	integer :: WLI(InBufferNumWL)     ! wake line index buffer

    ! Namelist input file declaration
	NAMELIST/ConfigInputs/RegTFlag,DiagOutFlag,GPFlag,FSFlag,nr,convrg,nti,iut,iWall,ivtxcor,VCRFB,VCRFT,VCRFS,ifc,convrgf,nric,ntif,iutf,ixterm,xstop, &
        Output_ELFlag,Output_DSFlag,WallOutFlag,Incompr,DSFlag,PRFlag, &
        k1pos,k1neg,GPGridSF,FSGridSF,TSFilFlag,ntsf
	NAMELIST/CaseInputs/jbtitle,GeomFilePath,RPM,Ut,nSect,AFDPath, &
        hAG,dFS,rho,vis,tempr,hBLRef,slex,Cdpar,CTExcrM, &
        WakeOutFlag,WLI,Igust,gustamp,gusttime,gustX0, &
        Itower,tower_Npts,tower_x,tower_ybot,tower_ytop,tower_D,tower_CD

    ! Input Defaults
    RegTFlag = 0 
    DiagOutFlag = 0
    Output_ELFlag = 0 
    Output_DSFlag = 0
    WakeOutFlag = 0 
    WallOutFlag = 0
    GPFlag=0
    FSFlag=0    
    TSFilFlag=0
    ntsf=3
	nSect = 1  
	ifc = 0   
	nr = 10
	convrg = -1
	convrgf = -1
	nti = 20
	ntif = -1
	ivtxcor = 1
	ixterm = 0
	xstop = 5.0
	iut = 0
    iWall = 0       
	iutf = 0 
	nric=-1
    VCRFB=1.0       
	VCRFT=1.0
    VCRFS=1.0       
    WLI(:)=0 ! all set to 0      
	hAG=0.0
    dFS=0.0                                                          
	Incompr=0                        
    Cdpar=0.0
    CTExcrM=0.0
    DSFlag=1
    PRFlag=1
    k1pos = 1.0                      
    k1neg = 0.5
    GPGridSF = 1.0
    FSGridSF = 1.0
    Igust = 0
    gustamp = 0.0
    gusttime = 0.0
    gustX0 = 0.0
    Itower = 0
    tower_Npts = 10
    tower_x = 0.0
    tower_ybot = 0.0
    tower_ytop = 0.0
    tower_D = 0.05
    tower_CD = 1.0

	! Namelist input
	read(4, nml=ConfigInputs) 
    read(4, nml=CaseInputs)                                                                                    

    ! Read geometry file
    Call InputGeom(GeomFilePath)

	! Set array bounds based on inputs
	! Geometry
	MaxBlades = nb
	MaxSegPerBlade = nbe
	MaxSegEndPerBlade = MaxSegPerBlade+1
	MaxSegEnds = MaxSegEndPerBlade*MaxBlades
    MaxSeg = MaxSegPerBlade*MaxBlades  
    MaxStruts = NStrut     
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
	CALL element_cns(MaxSegEnds,MaxSegEndPerBlade)     
	CALL airfoil_cns(MaxAOAVals,MaxReVals,MaxAirfoilSect)
	CALL wakedata_cns()
	CALL dystl_cns(MaxAirfoilSect,MaxReVals,MaxSegEnds)
    CALL output_cns(MaxSeg,MaxBlades,MaxStruts,DSFlag)       

	! Write from buffer...    
    WakeLineInd(1:NWakeInd)=WLI(1:NWakeInd)       

	! Set ground plane location for wall solution
	GPy=-hAG/Rmax

    ! Set depth and Froude number for free surface solution
    FSy=dFS/Rmax
    g=32.174  ! gravity, ft/s^2
    A=Rmax*(2.0*pi*rpm/60.0)**2/g   ! Accel ratio: w^2*R/g
    FnR=sqrt(A/ut**2)   ! Froude number based on turbine radius  FnR=Uinf/sqrt(g*R)	

	! Normalize ground shear inputs
	yref = hBLRef/Rmax  ! location of boundary layer edge (U/Uinf = 99% maybe) normalized to radius... 
	ygc  = hAG/Rmax   ! Ground clearance normalized to radius  

    ! Floor the tip speed ratio to the next lowest int
    ! and use this as the default update interval...
	if (iutf == 0) iutf = floor(ut)
	if (iut == 0) iut = floor(ut)
    if (iWall == 0) iWall = floor(ut)

    ! Default ntif to nti if nothing was input
    if (ntif .eq. -1) then
    	ntif=nti
	end if

    ! Set number of RHS evaluations to average for the free surface calculation (should cover approx 1 revolution)
    NFSRHSAve=nti/iWall

	ne = (nbe+1)*nb ! Total number of blade segment ends (over all blades)

    ! Timestep filter setup
    KTF=1.0/real(ntsf)

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
                ! Check for carriage return (len_trim doesn't consider this a blank)
                NotBlank=.TRUE.
                if (len_trim(ReadLine)==0) then
                    NotBlank=.FALSE.
                else if (len_trim(ReadLine)==1) then
                    if (ichar(ReadLine(len_trim(ReadLine):len_trim(ReadLine))) == 13) then
                        NotBlank=.FALSE.
                    end if
                end if
                if (EOF>=0 .AND. NotBlank) then
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
601 format(' ','***airfoil section specified for blade segment ',i2,' is illegal. set to airfoil section 1***')                        
End SUBROUTINE input




