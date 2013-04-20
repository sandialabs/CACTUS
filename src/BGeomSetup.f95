SUBROUTINE BGeomSetup() 

	use parameters
	
        use configr       
	use element      
	use pidef
        use util       
	
        integer i, j       
	real dx, dy, dz, vMag

	! Sets up blade geometry arrays for each blade, starting from the root of the first blade, and continuing for each blade.        
	
        ! Set reference cr
        CrRef=0.0
        do i=1,nb
            do j=1,(nbe+1) 
                CrRef=max(CrRef,Blades(i)%CtoR(j))
            end do  
        end do     
        
        ! Set representative geometry discretization level (for vortex core calculation)
        dSGeom=0.0
        do i=1,nb
            do j=1,nbe
                dx=Blades(i)%QCx(j+1)-Blades(i)%QCx(j)
                dy=Blades(i)%QCy(j+1)-Blades(i)%QCy(j)
                dz=Blades(i)%QCz(j+1)-Blades(i)%QCz(j)
                ds=sqrt(dx**2+dy**2+dz**2);
                dsGeom=max(dsGeom,ds)
            end do  
        end do  
           
	! Set initial turbine geometry (zero theta)                                                            
	do i=1,nb                                                      
		
		nei=1+(i-1)*(nbe+1) 
		
                ! Blade end locations (quarter chord). xBE(MaxSegEnds)
                do j=0,nbe   
                        nej=nei+j ! element index 

                        xBE(nej)=Blades(i)%QCx(j+1)                                       
                        yBE(nej)=Blades(i)%QCy(j+1)                                                                                              
                        zBE(nej)=Blades(i)%QCz(j+1) 
                end do 
                
                
		! JCM: currently, although these values are for each element, they are held in arrays sized for element ends, where the first value for each blade
		! is simply ignored in bsload (where these values are used).           
			
		! JCM: These zeros are ignored...
		eArea(nei)=0.0    
		eChord(nei)=0.0  
                iSect(nei)=0.0                 
                xBC(nei)=0.0                                                 
                yBC(nei)=0.0                                                 
                zBC(nei)=0.0               
                nxBC(nei)=0.0                                                 
                nyBC(nei)=0.0                                                 
                nzBC(nei)=0.0   
                txBC(nei)=0.0                                                 
                tyBC(nei)=0.0                                                 
                tzBC(nei)=0.0
                sxBC(nei)=0.0                                                 
                syBC(nei)=0.0                                                 
                szBC(nei)=0.0                                                                                                
		do j=1,nbe 
			nej=nei+j                                                                                                      
			eArea(nej)=Blades(i)%AreaR(j)
			eChord(nej)=0.5*(Blades(i)%CtoR(j)+Blades(i)%CtoR(j+1))
                        iSect(nej)=Blades(i)%iSect(j)                   
                        
                        ! Blade center locations (quarter chord). xBC(MaxSegEnds)
                        xBC(nej)=0.5*(Blades(i)%QCx(j)+Blades(i)%QCx(j+1))                                      
                        yBC(nej)=0.5*(Blades(i)%QCy(j)+Blades(i)%QCy(j+1))                                                                                              
                        zBC(nej)=0.5*(Blades(i)%QCz(j)+Blades(i)%QCz(j+1)) 
                        
                        ! QC line 
                        dx=Blades(i)%QCx(j+1)-Blades(i)%QCx(j)
                        dy=Blades(i)%QCy(j+1)-Blades(i)%QCy(j)
                        dz=Blades(i)%QCz(j+1)-Blades(i)%QCz(j)
                        
                        ! Set element tangent vectors, txBC(MaxSegEnd)  
                        txBC(nej)=0.5*(Blades(i)%tx(j)+Blades(i)%tx(j+1))
                        tyBC(nej)=0.5*(Blades(i)%ty(j)+Blades(i)%ty(j+1))                                                                             
                        tzBC(nej)=0.5*(Blades(i)%tz(j)+Blades(i)%tz(j+1))
                        ! Force normalize
                        VMag=sqrt(txBC(nej)**2+tyBC(nej)**2+tzBC(nej)**2) 
                        txBC(nej)=txBC(nej)/VMag 
                        tyBC(nej)=tyBC(nej)/VMag    
                        tzBC(nej)=tzBC(nej)/VMag 
                        
                        ! Set element normal vectors, nxBC(MaxSegEnd) 
                        nxBC(nej)=0.5*(Blades(i)%nx(j)+Blades(i)%nx(j+1))
                        nyBC(nej)=0.5*(Blades(i)%ny(j)+Blades(i)%ny(j+1))                                                                             
                        nzBC(nej)=0.5*(Blades(i)%nz(j)+Blades(i)%nz(j+1))
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
                    
                        ! Calculate element spanwise vector (s = t x n)
                        CALL cross(txBC(nej),tyBC(nej),tzBC(nej),nxBC(nej),nyBC(nej),nzBC(nej),sxBC(nej),syBC(nej),szBC(nej))     
                        ! Force normalize
                        VMag=sqrt(sxBC(nej)**2+syBC(nej)**2+szBC(nej)**2) 
                        sxBC(nej)=sxBC(nej)/VMag 
                        syBC(nej)=syBC(nej)/VMag    
                        szBC(nej)=szBC(nej)/VMag 
                        ! Calc projection direction of element spanwise vector on QC line.
                        ! Used to correctly orient circulation in wake grid.
                        ! Note pos. lift creates circulation in neg. element spanwise direction
                        dp=sxBC(nej)*dx+syBC(nej)*dy+szBC(nej)*dz
                        if (dp>0) then
                            CircSign(nej)=-1.0
                        else
                            CircSign(nej)=1.0
                        end if
                        
                        
		end do 
                
        end do
        
return                                                            
end 
