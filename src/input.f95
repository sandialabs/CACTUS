SUBROUTINE input(ErrFlag) 

	use parameters
	
	use dystl
	use element
	use blade
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
        use regtest
        use output            
	
	integer, parameter :: InBufferNumSectionTables = 100
	integer, parameter :: InBufferNumSegPerBlade = 100
	
	integer i
	integer iend
	integer ErrFlag
	
	! Temp buffers for airfoil section data inputs. If more sections or segments become necessary, its probably time to 
	! create a more generic geometry file with all the relevant geometry for one blade defined, and read this file in (using a loop).
	character*80 :: AFDPath(InBufferNumSectionTables)	! Airfoil section data path
	integer :: iSection(InBufferNumSegPerBlade)	! section index for each element
	real :: ChR(InBufferNumSegPerBlade)	! chord ratio for each element
	real :: bTwist(InBufferNumSegPerBlade)	! twist for each element (currently only used for HAWT geometry)
	
                
	! Namelist input file declaration
	NAMELIST/ConfigInputs/RegTFlag,GeomFlag,GPFlag,rho,vis,tempr,hFSRef,slex,nr,convrg,nti,iut,ivtxcor,ifwg,ifc,convrgf,nric,ntif,iutf,ixterm,xstop,Output_ELFlag,Incompr
	NAMELIST/XFlowInputs/jbtitle,Rmax,RPM,Ut,CrRef,ChR,hr,eta,nb,nbe,nSect,AFDPath,iSection,hAG,Istraight,Istrut,sThick
	NAMELIST/AxFlowInputs/jbtitle,R,HubR,RPM,Ut,Tilt,CrRef,ChR,bCone,bi,bTwist,eta,nb,nbe,nSect,AFDPath,iSection,hAG
	
	! Input Defaults
        RegTFlag = 0 
        Output_ELFlag=0      
	nb = 2
	nbe = 5 
	nSect = 1 
	ifwg = 0 
	ifc = 0   
	nr = 10
	convrg = .0001
	convrgf = .0001
	nti = 16
	ntif = 16
	ivtxcor = 0
	ixterm = 0
	xstop = 5.0
	iut = 0
	iutf = 0 
	nric=-1
	CrRef=0
	iSection(:)=1 ! all set to 1 
	bTwist(:)=0.0 ! all set to 0
	ChR(:)=0.0 ! all set to 0
	hAG=0.0
	Tilt=0.0
	bCone=0.0 
	HubR=0.0                                                    
	Incompr=0
	Istrut=0
	sThick=0.0                         
	                                                                              
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
	MaxReVals = 15
	MaxAOAVals = 400
	! Wake advancement
	MaxRevs = nr
	MaxTimeStepPerRev = nti
	MaxWakeNodes = MaxRevs * MaxTimeStepPerRev   
	! Fixed wake grid 
	MaxFixWakeX = 12
	MaxFixWakeY = 20
	MaxFixWakeZ = 20
	! Non-linear convergence iteration
	MaxNLIters = 10
        ! Outputs
        MaxTimeSteps = MaxRevs * MaxTimeStepPerRev       
	
	! Array construction
        CALL blade_cns(MaxWakeNodes,MaxSegEnds)
	CALL element_cns(MaxTimeStepPerRev,MaxSegEnds,MaxSegEndPerBlade)     
	CALL cltab_cns(MaxAOAVals,MaxReVals,MaxAirfoilSect,MaxSegEnds)
	CALL airfoil_cns(MaxAirfoilSect)
	CALL xwake_cns(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ)
	CALL uwake_cns(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ)
	CALL wakedata_cns(MaxFixWakeY,MaxFixWakeZ)
	CALL freestream_cns(MaxWakeNodes,MaxSegEnds)
	CALL dystl_cns(MaxAirfoilSect,MaxReVals,MaxSegEnds)
        CALL output_cns(MaxRevs, MaxTimeSteps, MaxSeg, MaxBlades)       
        
	! Write from buffer...
	iSect(1:(nbe+1))=iSection(1:(nbe+1))
	cr(1:(nbe+1))=ChR(1:(nbe+1))
	btw(1:(nbe+1))=bTwist(1:(nbe+1))
	
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
	if (CrRef == 0) CrRef=cr(1)
	if (cr(2) == 0) then
		do i = 1,(nbe+1)
			cr(i)=cr(1)
		end do 
	end if
	
	if (ivtxcor == 0) then
		vrad = cr(1)       
	else
		vrad = 0.0  
	end if                                    
	vrad2 = vrad*vrad
	                  				
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

	
        ! Airfoil Stall Data Tables                         

	! Read stall angle data, and section CL and CD data from section data file 
	do kk = 1, nsect 
	
		! Open input file for this section
		open(15, file=AFDPath(kk))
	
		! Read in title on first line of stall data block and (char,real,real) data on second line
		read(15,420) aftitle(kk)                                  
		read(15,*) camber(kk),tc(kk),alzer(kk)
		
		! Create numerical flag for camber in (1), out (-1), or none (0) from text flag
		if ((camber(kk)(1:2) == 'in') .OR. (camber(kk)(1:2) == 'IN')) then
			camb(kk) = 1 
		else if ((camber(kk)(1:3) == 'out') .OR. (camber(kk)(1:3) == 'OUT')) then
			camb(kk) = -1         
		else
			camb(kk) = 0
		end if	
		
		! Read in stall angles (positive and negative) at up to MaxReVals Re values, starting on line 3 of stall data block
		nstl(kk)=0  ! number of stall angle data points 
		iend=0
		i=1
		do while ((iend /= 1) .AND. (i <= MaxReVals)) 
			read(15,*) iend,restl(i,kk),alstlp(i,kk),alstln(i,kk)
			nstl(kk) = i  
			i=i+1   
		end do                                                  

	        ! If airfoil is cambered out, transfer stall information 
	        ! to proper sign. Calculate slopes for interpolation.     
		alzer(kk) = alzer(kk)*conrad
		if (camb(kk) == -1) alzer(kk) = -alzer(kk) 
		
		istp = nstl(kk) 
		do i = 1, istp  
			if(camb(kk) == -1) then
				temp = alstlp(i,kk)
				alstlp(i,kk) = -alstln(i,kk)
				alstln(i,kk) = -temp   
			end if                                             
			alstlp(i,kk) = alstlp(i,kk)*conrad 
			alstln(i,kk) = alstln(i,kk)*conrad 
			if (i > 1) then
				dapdre(i-1,kk) = (alstlp(i,kk)-alstlp(i-1,kk))/(restl(i,kk)-restl(i-1,kk))         
				dandre(i-1,kk) = (alstln(i,kk)-alstln(i-1,kk))/(restl(i,kk)-restl(i-1,kk))  
			end if   
		end do 

		
		! Read CL and CD data for this airfoil as a func. of Re
		
		! Read section title and first Re number for this airfoil section. 
		read(15,1001) tre(1,kk), dftitle(kk)
	
		! Read in section lift and drag data for up to 15 Re numbers
		iend=0
		i=1
		do while ((iend < 2) .AND. (i <= MaxReVals))
			! Read Re value if it hasn't already been read (if not the first Re block of airfoil)
			if (i > 1) then                                           
				read(15,*) tre(i,kk)
			end if
			
			! Read up to MaxAOAVals lift and drag values as a function of AOA for this Re number
			ii=1
			ntb=0  
			iend=0 
			do while ((iend == 0) .AND. (ii <= MaxAOAVals))
				ntb = ntb+1  
				! if iend flag is 1, then done with this Re number
				! if iend flag is 2, done with entire section table 
				read(15,*) iend, ta(ii,i,kk), tcl(ii,i,kk),tcd(ii,i,kk) 
				ii=ii+1
			end do                
	                                                                       
			! Adjust the airfoil data so alpha goes from -180 to 180             
			! and the cambered data is of the proper sign for this rotor blade   					
			if(camb(kk) == -1) then
				
				npt = ntb/2 
				do ii = 1, npt  
					jj    = ntb-(ii-1)
					temp1 = ta(ii,i,kk) 
					temp2 = tcl(ii,i,kk)
					temp3 = tcd(ii,i,kk)
					ta(ii,i,kk)  = -ta(jj,i,kk)
					tcl(ii,i,kk) = -tcl(jj,i,kk)
					tcd(ii,i,kk) =  tcd(jj,i,kk) 
					ta(jj,i,kk)  = -temp1 
					tcl(jj,i,kk) = -temp2
					tcd(jj,i,kk) =  temp3
				end do  
				
				if ((npt*2) /= ntb) then
					ii            = npt+1    
					ta(ii,i,kk)  = -ta(ii,i,kk)
					tcl(ii,i,kk) = -tcl(ii,i,kk)
				end if
				
			else if(camb(kk) == 0) then
				
				npt = ntb*2-1 
				do ii = 1, ntb  
					jj = npt-(ii-1)
					k1 = ntb-(ii-1) 
					ta(jj,i,kk)  = ta(k1,i,kk)
					tcl(jj,i,kk) = tcl(k1,i,kk) 
					tcd(jj,i,kk) = tcd(k1,i,kk)
				end do  
				     
				npts = ntb-1  
				do ii = 1, npts 
					jj = ntb+ii   
					k1 = npts-(ii-1)
					ta(k1,i,kk)  = -ta(jj,i,kk)  
					tcl(k1,i,kk) = -tcl(jj,i,kk)  
					tcd(k1,i,kk) =  tcd(jj,i,kk)
				end do   
				     
				ntb = npt 
				  
			end if  
			 
			ntbl(i,kk) = ntb  
		
			i=i+1

		end do

		! Close input file for this section
		close(15)

	end do
									
Return  	
420  format(a80)  								                                        
601  format(' ','***airfoil section specified for blade segment ',i2,' is illegal. set to airfoil section 1***') 
1001 format(e10.1,a64)                       
End 
							
