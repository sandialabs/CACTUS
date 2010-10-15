SUBROUTINE BGeomSetup_v(delty,delt,deltb,at) 

	use parameters
        
        use configr 	
	use element
	use pidef
	
	integer :: i, j, k, nei, nej, nej1
	real :: delty, delt, deltb, deltac, at
	real :: dx, dy, dz, NMag
	real :: rr(MaxSegEndPerBlade)     ! Blade r/R at segment ends 
	real :: yB(MaxSegEndPerBlade)     ! Blade y/R at segment ends
	real :: rSE
	
	! Sets up blade geometry arrays for each blade at each theta position in one revolution, starting from the bottom of the first blade, and continuing for each blade.        
	                                           
	
	! Evaluate the VAWT radius ratio function of height ratio for a parabolic or straight blade shape at blade segment ends.
	! JCM: Note that this could be input from a file if an arbitrary blade shape is desired...      
	deltac=(eta-.25)*cr(1)    ! define mount point w.r.t. root quarter chord point                          

	If (Istraight .EQ. 1) Then	! Straight blades
		do j=1,(nbe+1)  
			yB(j)=(j-1)*delty 
			rr(j)=1.0
		end do 		
	Else				! Parabolic blades
		do j=1,(nbe+1)  
			yB(j)=(j-1)*delty 
			rr(j)=1.0-4.0*(real(j-1)/nbe-0.5)**2                                                                        
		end do 
	End If
	
	! Frontal area normalized by Rmax^2 
	at=0.0                                                                                                                                                                          
	do j=2,(nbe+1)                                                                                               
		at=at+(rr(j)+rr(j-1))/2.0*(yB(j)-yB(j-1))                                                    
	end do                                                         
	at=2.0*at  
	     
	! Create the turbine blade geometry for each blade                                                             
	do i=1,nb                                                      
		
		nei=1+(i-1)*(nbe+1) 
		
		! JCM: currently, although these values are for each element, they are held in arrays sized for element ends, where the first value for each blade
		! is simply ignored in bsload (where these values are used). This is kind of stupid and should be changed...            
			
		! JCM: These zeros are ignored...
		eSpan(nei)=0.0   ! Element length
		eChord(nei)=0.0   ! Element chord                                                                                    
		do j=1,nbe 
			nej=nei+j                                                                                                      
			eSpan(nej)=sqrt((rr(j+1)-rr(j))**2+(yB(j+1)-yB(j))**2)
			eChord(nej)=0.5*(cr(j)+cr(j+1))
		end do 
		
		! For each theta position
		do k=1,nti                              
			
			! Arrays of sin(theta) and cos(theta) for each theta and each blade
			thetaB1=(k-1)*delt 
			thetaB=thetaB1+(i-1)*deltb     ! blade theta
			sint=sin(thetaB)    
			cost=cos(thetaB)
			Theta(k)=thetaB1*condeg                 
						
			! Blade end locations (quarter chord). xBE(MaxTimeStepPerRev,MaxSegEnds)
			do j=0,nbe   
				nej=nei+j ! element index                                          
				yBE(k,nej)=yB(j+1)                                                    
				xBE(k,nej)=-rr(j+1)*sint-deltac*cost                                            
				zBE(k,nej)=-rr(j+1)*cost+deltac*sint    
			end do
                        If (Istrut.EQ.1) Then
                           ! Blade strut element end locations (MFB)
                           Do j=0,nbe
                              nej=nei+j ! element index
                              rSE = rr(j+1) * REAL(j) / REAL(nbe)
                              ySE(k,nej)= 0.5 * hr ! struts are located at blade mid-span location
                              xSE(k,nej) = -rSE*sint
                              zSE(k,nej) = -rSE*cost
                           End Do
                        End If
			        
			! JCM: currently, although these values are for each element, they are held in arrays sized for element ends, where the first value for each blade
			! is simply ignored in bsload (where these values are used). This is kind of stupid and should be changed...            
			
			! JCM: These zeros are ignored..l.
			nx(k,nei)=0.0                                                 
			ny(k,nei)=0.0                                                 
			nz(k,nei)=0.0    
			tx(k,nei)=0.0                                                 
			ty(k,nei)=0.0                                                 
			tz(k,nei)=0.0                                              
			! Calculate the normal (machine inward) and tangential vectors (rearward chord line) for each element. nx(MaxTimeStepPerRev,MaxSegEnd)                      
			do j=1,nbe                                                     
				nej=nei+j                                                         
				nej1=nej-1 
				tx(k,nej)=cost                                                      
				ty(k,nej)=0.0                                                                            
				tz(k,nej)=-sint 
				dx=(xBE(k,nej)-xBE(k,nej1))/eSpan(nej)                                                       
				dy=(yBE(k,nej)-yBE(k,nej1))/eSpan(nej)                                                                              
				dz=(zBE(k,nej)-zBE(k,nej1))/eSpan(nej) 
				CALL cross(tx(k,nej),ty(k,nej),tz(k,nej),dx,dy,dz,nx(k,nej),ny(k,nej),nz(k,nej))     
				! Force normalize the normal vector (in case t and dChord aren't perp.)
				NMag=sqrt(nx(k,nej)**2+ny(k,nej)**2+nz(k,nej)**2) 
				nx(k,nej)=nx(k,nej)/NMag 
				ny(k,nej)=ny(k,nej)/NMag    
				nz(k,nej)=nz(k,nej)/NMag                           
			end do 
			
		end do 
	end do                                                         
	         
return                                                            
end 
