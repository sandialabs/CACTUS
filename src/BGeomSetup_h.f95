SUBROUTINE BGeomSetup_h(deltr,delt,deltb) 

	use parameters
	
        use configr       
	use element
	use pidef
	
	real :: deltr, delt, deltb
	
	real BladeQC(MaxSegEndPerBlade,3), tB(MaxSegEndPerBlade,3), normB(MaxSegEndPerBlade,3)
	real rB, cSweep, sSweep, dy, dz, deltac
	real RGRo(3,3), RRoP(3,3), RPC(3,3), RCB(3,3), RGB(3,3)	! Rotation matricies
	real BladeQCG(3,1), tG(3,1), nG(3,1)

	! Sets up blade geometry arrays for each blade at each theta position in one revolution, starting from the root of the first blade, and continuing for each blade.        
	
	BladeQC(:,:)=0.0
	tB(:,:)=0.0
	normB(:,:)=0.0
	
	! Setup quarter chord positions and vectors for each element end in blade (planform) axes.                                   
	! Blade (planform) axes have chordwise direction (toward LE) in z, span direction in y, and normal in x 
	! Pos. twist is defined LE toward negative normal direction (-x) 
	
	! Quarter chord positions
	deltac=(eta-.25)*cr(1)    ! define mount point w.r.t. root quarter chord point 
	do j=1,(nbe+1) 
		rB=(j-1)*deltr+hubrr  ! Blade r/R at segment ends
		BladeQC(j,1:3)=[0.0,rB,deltac] ! blade quarter chord locations at segment ends w.r.t. blade mount point in blade (planform) axes
	end do
	
	! Define chord vector (rearward, toward TE) from sweep and twist (around quarter chord line) in blade axes
	! Define normal vector (machine rearward at zero incidence) in blade axes 
	! Sweep used for root section is the same as the first outboard section
	do j=1,(nbe+1)
		if (j>1) then
			dy=BladeQC(j,2)-BladeQC(j-1,2)
			dz=BladeQC(j,3)-BladeQC(j-1,3)	
		else
			dy=BladeQC(2,2)-BladeQC(1,2)
			dz=BladeQC(2,3)-BladeQC(1,3)
		end if
		cSweep=dy/sqrt(dy**2+dz**2)
		sSweep=-dz/sqrt(dy**2+dz**2)
		sTwist=sin(btw(j)*conrad)
		cTwist=cos(btw(j)*conrad)
		tB(j,1:3)=[sTwist,-sSweep*cTwist,-cSweep*cTwist]  
		normB(j,1:3)=[cTwist,sSweep*sTwist,cSweep*sTwist]                                                              
	end do 
	
	
	! Normalized frontal area	
	at=pi     
	     
	! Static rotation matricies
	sinti=sin(Tilt*conrad)    
	costi=cos(Tilt*conrad)
	RGRo=reshape([costi,-sinti,0.0,sinti,costi,0.0,0.0,0.0,1.0],[3,3])	! to Global from Rotor (rotor tilt rotation)     
	
	sinc=sin(bCone*conrad)    
	cosc=cos(bCone*conrad)
	RPC=reshape([cosc,sinc,0.0,-sinc,cosc,0.0,0.0,0.0,1.0],[3,3])	! to Prop from Case (coning angle rotation) 
	
	sini=sin(bi*conrad)    
	cosi=cos(bi*conrad)
	RCB=reshape([cosi,0.0,sini,0.0,1.0,0.0,-sini,0.0,cosi],[3,3])	! to Case to Blade (incidence rotation)
	     
	! Create the turbine blade geometry for each blade                                                             
	do i=1,nb                                                      
		
		nei=1+(i-1)*(nbe+1) 
		
		! JCM: currently, although these values are flor each element, they are held in arrays sized for element ends, where the first value for each blade
		! is simply ignored in bsload (where these values are used). This is kind of stupid and should be changed...            
			
		! JCM: These zeros are ignored...
		eSpan(nei)=0.0   ! Element length  
		eChord(nei)=0.0   ! Element chord                                                                                
		do j=1,nbe 
			nej=nei+j                                                                                                      
			eSpan(nej)=deltr
			eChord(nej)=0.5*(cr(j)+cr(j+1))
		end do 
		
		! For each theta position
		do k=1,nti                              
			
			thetaB1=(k-1)*delt 
			thetaB=thetaB1+(i-1)*deltb     ! blade theta
			sint=sin(thetaB)    
			cost=cos(thetaB)
			RRoP=reshape([1.0,0.0,0.0,0.0,cost,sint,0.0,-sint,cost],[3,3])	! to Rotor from Prop (theta rotation, positive around +x)
			Theta(k)=thetaB1*condeg                 
					
			! Rotation to global from blade 
			RGB=matmul(RGRo,matmul(RRoP,matmul(RPC,RCB)))	
						
			! Blade end locations (quarter chord) in global axes. xBE(MaxTimeStepPerRev,MaxSegEnds)
			! The origin point should be on the axis of rotation.
			do j=0,nbe   
				nej=nei+j ! element index 
				
				BladeQCG=matmul(RGB,reshape(BladeQC(j+1,1:3),[3,1]))  ! rotated to global axes
				xBE(k,nej)=BladeQCG(1,1)                                       
				yBE(k,nej)=BladeQCG(2,1)                                                                                              
				zBE(k,nej)=BladeQCG(3,1) 
				
			end do                                         
                          
			        
			! JCM: currently, although these values are for each element, they are held in arrays sized for element ends, where the first value for each blade
			! is simply ignored in bsload (where these values are used). This is kind of stupid and should be changed...            
			
			! JCM: These zeros are ignored...
			nx(k,nei)=0.0                                                 
			ny(k,nei)=0.0                                                 
			nz(k,nei)=0.0   
			tx(k,nei)=0.0                                                 
			ty(k,nei)=0.0                                                 
			tz(k,nei)=0.0                                              
			! Calculate the normal (machine rearward) vectors for each element. nx(MaxTimeStepPerRev,MaxSegEnd)  
			! Calculate the tangential vectors (rearward chord line) for each element. tx(MaxTimeStepPerRev,MaxSegEnd)                      
			do j=1,nbe                                                     
				nej=nei+j                                                         

				! Chord vector rotated to global
				tG=matmul(RGB,reshape(0.5*(tB(j,1:3)+tB(j+1,1:3)),[3,1]))
 				tG=tG/sqrt(sum(tG**2))
				tx(k,nej)=tG(1,1)
				ty(k,nej)=tG(2,1)                                                                             
				tz(k,nej)=tG(3,1) 
				
				! Normal vector rotated to global
				nG=matmul(RGB,reshape(0.5*(normB(j,1:3)+normB(j+1,1:3)),[3,1]))
				nG=nG/sqrt(sum(nG**2))
				nx(k,nej)=nG(1,1)
				ny(k,nej)=nG(2,1)                                                                             
				nz(k,nej)=nG(3,1) 
 
			end do 
			
		end do 
	end do                                                         
	         
return                                                            
end 
