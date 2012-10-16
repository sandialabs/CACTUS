MODULE element

	! Blade element geometry data
        
        use util
        
        ! Blade geometry input structure
        ! JCM: Currently used only to store input geometry file data and blade loads outputs. Should eventually 
        ! replace the arrays below that are used for the internal calculation, similar to the strut module. This will 
        ! require a change to the blade/element/wake iterators throughout the entire code...
        type BladeType
            integer :: NElem
            real, allocatable :: QCx(:)
            real, allocatable :: QCy(:)
            real, allocatable :: QCz(:)
            real, allocatable :: nx(:)
            real, allocatable :: ny(:)
            real, allocatable :: nz(:)
            real, allocatable :: tx(:)
            real, allocatable :: ty(:)
            real, allocatable :: tz(:)
            real, allocatable :: CtoR(:)
            real, allocatable :: AreaR(:)
            integer, allocatable :: iSect(:)
            
            ! Current total blade output (nomalized by machine scale parameters)
            real :: CP  ! Power coefficient due to this blade
            real :: CTR ! Torque coefficient due to this blade
            real :: CFx ! Fx coefficient due to this blade
            real :: CFy ! Fy coefficient due to this blade
            real :: CFz ! Fz coefficient due to this blade
        end type
            
        type(BladeType), allocatable :: Blades(:)   ! Input blade geometry
	
	real, allocatable :: xBE(:)		! X location for each blade segment end (quarter chord) 	
	real, allocatable :: yBE(:)		! Y location for each blade segment end (quarter chord)		
	real, allocatable :: zBE(:)		! Z location for each blade segment end (quarter chord)	
        
        real, allocatable :: xBC(:)             ! X location for each blade segment center (quarter chord)         
        real, allocatable :: yBC(:)             ! Y location for each blade segment center (quarter chord)         
        real, allocatable :: zBC(:)             ! Z location for each blade segment center (quarter chord) 

	real, allocatable :: xSE(:)		! X location for each strut segment end
	real, allocatable :: ySE(:)		! Y location for each strut segment end
	real, allocatable :: zSE(:)		! Z location for each strut segment end	

      	real :: dSGeom                          ! Geometry discretization level used in vortex core calculation
        real :: CrRef                           ! Ref chord to radius ratio 	
	real, allocatable :: nxBC(:)		! Normal X for each blade segment
	real, allocatable :: nyBC(:)		! Normal Y for each blade segment
	real, allocatable :: nzBC(:)		! Normal Z for each blade segment 
	real, allocatable :: txBC(:)		! Tangential X for each blade segment 	
	real, allocatable :: tyBC(:)		! Tangential Y for each blade segment 	
	real, allocatable :: tzBC(:)		! Tangential Z for each blade segment 	
        real, allocatable :: sxBC(:)		! Spanwise X for each blade segment 	
	real, allocatable :: syBC(:)		! Spanwise Y for each blade segment 	
	real, allocatable :: szBC(:)		! Spanwise Z for each blade segment 	
	real, allocatable :: eArea(:)		! Element area to radius ratio for each element
	real, allocatable :: eChord(:)		! Element chord to radius ratio for each element
        integer, allocatable :: iSect(:)        ! Array of indicies of the section table to apply to each blade element
        
        ! Current sum of output over all blades (nomalized by machine scale parameters)
        real :: CP_B  ! Power coefficient due to all blades
        real :: CTR_B ! Torque coefficient due to all blades
        real :: CFx_B ! Fx coefficient due to all blades
        real :: CFy_B ! Fy coefficient due to all blades
        real :: CFz_B ! Fz coefficient due to all blades

	CONTAINS


        SUBROUTINE blade_geom_cns(BInd,NElem)

                ! Constructor for the arrays in this module

                integer :: BInd, NElem
        
                Blades(BInd)%NElem=NElem
                allocate(Blades(BInd)%QCx(NElem+1))    
                allocate(Blades(BInd)%QCy(NElem+1)) 
                allocate(Blades(BInd)%QCz(NElem+1)) 
                allocate(Blades(BInd)%nx(NElem+1)) 
                allocate(Blades(BInd)%ny(NElem+1)) 
                allocate(Blades(BInd)%nz(NElem+1)) 
                allocate(Blades(BInd)%tx(NElem+1)) 
                allocate(Blades(BInd)%ty(NElem+1)) 
                allocate(Blades(BInd)%tz(NElem+1)) 
                allocate(Blades(BInd)%CtoR(NElem+1)) 
                allocate(Blades(BInd)%AreaR(NElem)) 
                allocate(Blades(BInd)%iSect(NElem))         
                
        End SUBROUTINE


	SUBROUTINE element_cns(MaxSegEnds,MaxSegEndPerBlade)

		! Constructor for the arrays in this module

		integer :: MaxSegEnds,MaxSegEndPerBlade
		       
		allocate(xBE(MaxSegEnds))
		allocate(yBE(MaxSegEnds))
		allocate(zBE(MaxSegEnds))
                allocate(xBC(MaxSegEnds))
                allocate(yBC(MaxSegEnds))
                allocate(zBC(MaxSegEnds))              
		allocate(xSE(MaxSegEnds))
		allocate(ySE(MaxSegEnds))
		allocate(zSE(MaxSegEnds))
		allocate(nxBC(MaxSegEnds))		
		allocate(nyBC(MaxSegEnds))		
		allocate(nzBC(MaxSegEnds))		
		allocate(txBC(MaxSegEnds))		
		allocate(tyBC(MaxSegEnds))		
		allocate(tzBC(MaxSegEnds))		
		allocate(sxBC(MaxSegEnds))		
		allocate(syBC(MaxSegEnds))		
		allocate(szBC(MaxSegEnds))	        
		allocate(eArea(MaxSegEnds))	
		allocate(eChord(MaxSegEnds))			
                allocate(iSect(MaxSegEnds))              
		
	End SUBROUTINE
        
        
        SUBROUTINE RotateBlade(bind,delt,nrx,nry,nrz,px,py,pz)

                integer :: bind
                real :: delt,nrx,nry,nrz,px,py,pz 
                integer :: nbe
                real :: vrx,vry,vrz,VMag

                ! Rotates data in blade arrays                                                    
                
                ! JCM: Eventually, should just be able to loop through Blades(bind) data structure
                ! While data is still held in arrays concatenated across blades, need to replicate
                ! nbe (stored in configr) from Blades(1).NElem
                nbe=Blades(1)%NElem
                
                nei=1+(bind-1)*(nbe+1) 
                
                ! Blade end locations (quarter chord). xBE(MaxSegEnds)
                do j=0,nbe   
                        nej=nei+j ! element index 

                        Call QuatRot(xBE(nej),yBE(nej),zBE(nej),delt,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
                        xBE(nej)=vrx                                       
                        yBE(nej)=vry                                                                                              
                        zBE(nej)=vrz 
                end do 
                
                                                                                                               
                do j=1,nbe 
                        nej=nei+j                                                                                                      
                        
                        ! Blade center locations (quarter chord). xBC(MaxSegEnds)
                        Call QuatRot(xBC(nej),yBC(nej),zBC(nej),delt,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
                        xBC(nej)=vrx                                      
                        yBC(nej)=vry                                                                                             
                        zBC(nej)=vrz 
                        
                        ! Set element tangent vectors, txBC(MaxSegEnd)  
                        Call QuatRot(txBC(nej),tyBC(nej),tzBC(nej),delt,nrx,nry,nrz,0.0,0.0,0.0,vrx,vry,vrz)
                        txBC(nej)=vrx
                        tyBC(nej)=vry   
                        tzBC(nej)=vrz
                        ! Force normalize
                        VMag=sqrt(txBC(nej)**2+tyBC(nej)**2+tzBC(nej)**2) 
                        txBC(nej)=txBC(nej)/VMag 
                        tyBC(nej)=tyBC(nej)/VMag    
                        tzBC(nej)=tzBC(nej)/VMag 
                        
                        ! Set element normal vectors, nxBC(MaxSegEnd) 
                        Call QuatRot(nxBC(nej),nyBC(nej),nzBC(nej),delt,nrx,nry,nrz,0.0,0.0,0.0,vrx,vry,vrz)
                        nxBC(nej)=vrx
                        nyBC(nej)=vry   
                        nzBC(nej)=vrz
                        ! Force normal to t
                        dp=txBC(nej)*nxBC(nej)+tyBC(nej)*nyBC(nej)+tzBC(nej)*nzBC(nej)
                        nxBC(nej)=nxBC(nej)-dp*txBC(nej)
                        nyBC(nej)=nyBC(nej)-dp*tyBC(nej)
                        nzBC(nej)=nzBC(nej)-dp*tzBC(nej)
                        ! Force normalize
                        VMag=sqrt(nxBC(nej)**2+nyBC(nej)**2+nzBC(nej)**2) 
                        nxBC(nej)=nxBC(nej)/VMag 
                        nyBC(nej)=nyBC(nej)/VMag    
                        nzBC(nej)=nzBC(nej)/VMag 
                    
                        ! Calculate spanwise vector (s = t x n)
                        CALL cross(txBC(nej),tyBC(nej),tzBC(nej),nxBC(nej),nyBC(nej),nzBC(nej),sxBC(nej),syBC(nej),szBC(nej))     
                        ! Force normalize
                        VMag=sqrt(sxBC(nej)**2+syBC(nej)**2+szBC(nej)**2) 
                        sxBC(nej)=sxBC(nej)/VMag 
                        syBC(nej)=syBC(nej)/VMag    
                        szBC(nej)=szBC(nej)/VMag 
                        
                end do 
                           
                
        End SUBROUTINE
	
End
