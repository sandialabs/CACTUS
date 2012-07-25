subroutine WGeomSetup() 

        use wallsoln 
        use pidef        
        
        real :: PlaneExtent, dsMinFS, dsMaxFS, dsMinW, dsMaxW, dsC, dsOut, GCS, yPanFS, yPan
        real :: RC, C1, C2, dB, RU, C1U, C2U, dBU, RD, C1D, C2D, dBD, RW, C1W, C2W, dBW, B
        real :: P1(3), P2(3), P3(3), P4(3)
        real, allocatable :: xPanGP(:), zPanGP(:)
        real, allocatable :: xPanFS(:), zPanFS(:)       
        integer :: Ind, NU, ND, NW      

        ! Plane extent (over radius). To be applied in every direction around the turbine location...
        PlaneExtent=10.0 
        
        ! Min and max panel size limits (Note: user can override these with scale factor input...)
        dsMinFS=.25 ! min grid size for free surface grid (to keep problem from getting too large)
        dsMaxFS=1.0 ! max free surface grid size allowed
        dsMinW=.05  ! min grid size for wall grid (to keep problem from getting too large)
        dsMaxW=.1   ! max wall grid size allowed

        if (GPFlag==1) then
        
            ! Setup ground plane grid, clustered in the center of the plane
        
            ! Set near-field grid scale (panel dimension over radius in the center). If h>R use h, else use R. User can modify 
            ! the default scale level with the GPGridSF parameter.
            if (abs(GPy)<1) then
                dsN=.1      
            else
                dsN=.1*abs(GPy) 
            end if  
            dsC=max(dsMinW,min(dsMaxW,dsN))
            ! Apply user specified grid scale factor
            dsC=dsC*GPGridSF
            
            RC=30.0 ! Clustering ratio, must be >= 1 (Ex: 10 for 10x density at center of plane)
            C1=1.0/(1.0+3.0/RC)
            C2=3.0/RC*C1
            dB=dsC/C2
            NumWPx=ceiling(2.0*PlaneExtent/dB)
            NumWPz=NumWPx
            NumWP=NumWPx*NumWPz
            
            ! Allocate local and global arrays
            allocate(xPanGP(NumWPx+1),zPanGP(NumWPz+1))
            Call wallsoln_gp_cns()
            
            dB=2.0/real(NumWPx)
            do i=1,NumWPx+1
                    B=-1.0+dB*(i-1)
                    xPanGP(i)=PlaneExtent*(C1*B**3+C2*B)
            end do
            zPanGP=xPanGP
            yPan=GPy
            
            ! Ordered list of panels going through x, then stepping y, repeat
            Ind=1
            do i=1,NumWPz
                    do j=1,NumWPx
                            P1=[xPanGP(j),yPan,zPanGP(i)]   	! lower panel x, lower panel y corner point
                            P2=[xPanGP(j+1),yPan,zPanGP(i)] 	! pos panel x neighbor point
                            P3=[xPanGP(j),yPan,zPanGP(i+1)]	! pos panel y neighbor point
                            P4=[xPanGP(j+1),yPan,zPanGP(i+1)]  	! pos panel x, pos panel y neighbor point
                            WCPoints(Ind,1:3)=0.25*(P1+P2+P3+P4) 	! panel center
                            WXVec(Ind,1:3)=P2-P1		! panel x tangential vector
                            WPL(Ind)=sqrt(sum(WXVec(Ind,1:3)**2))		! panel x length
                            WXVec(Ind,1:3)=WXVec(Ind,1:3)/WPL(Ind)		! normalize
                            WYVec(Ind,1:3)=P1-P3		! panel y tangential vector, set so that panel normal will be in the domain inward direction
                            WPW(Ind)=sqrt(sum(WYVec(Ind,1:3)**2))		! panel y length
                            WYVec(Ind,1:3)=WYVec(Ind,1:3)/WPW(Ind)		! normalize
                            Call cross(WXVec(Ind,1),WXVec(Ind,2),WXVec(Ind,3),WYVec(Ind,1),WYVec(Ind,2),WYVec(Ind,3),WZVec(Ind,1),WZVec(Ind,2),WZVec(Ind,3))	! panel normal vector	
                            WZVec(Ind,1:3)=WZVec(Ind,1:3)/sqrt(sum(WZVec(Ind,1:3)**2))		! normalize	
                            Ind=Ind+1
                    end do
            end do
    
            ! Panel edge tolerance (needs to be less than 1/2 of the min panel dimension)
            WEdgeTol=min(minval(WPL),minval(WPW))/10.0
        
        end if
        
        
        if (FSFlag==1) then
        
            ! Setup free surface grid. 
            ! Set near field grid scale based on depth when d/R > 1 or based on radius when d/R < 1.
            ! Set far-field grid scale based on far-field deep water transverse wavelength for given inflow conditions (Froude number).
            ! At low Froude number, the far-field grid scale will drop below the near field scale, at which point, the free surface
            ! boundary condition is de-activated and replaced with a wall boundary condition. Note that as Froude -> 0, the free 
            ! surface boundary condition becomes equivalent to the wall boundary condition...
            
            ! Set near-field grid scale (panel dimension over radius in the center). If d>R use d, else use R. User can modify 
            ! the default scale level with the FSGridSF parameter.
            if (abs(FSy)<1) then
                dsN=.1    
            else
                dsN=.1*abs(FSy)  
            end if 
            dsC=max(dsMinFS,min(dsMaxFS,dsN))
            ! Apply user specified grid scale factor
            dsC=dsC*FSGridSF 
            
            ! Calc far field length scale (scaled by deep water transverse wavelength)
            GCS=1.0/6.0
            dsOut=GCS*2*pi*FnR**2  ! outer ds required to capture wave train 
            dsOut=min(dsMaxFS,dsOut)
            
            ! If dsOut is less than dsC, use wall solution instead of free surface. The
            ! wall solution is the low Froude limit of the free surface solution and
            ! does not require a small grid size to capture the low wavelength modes of
            ! the free surface at low Froude number.
            UseFSWall=.FALSE.
            if (dsOut<dsC) then
                UseFSWall=.TRUE.
                ! Reset grid level
                dsC=max(dsMinW,min(dsMaxW,dsN)) 
            end if
            
            ! Setup geometry
            if (UseFSWall) then
            
                ! Use wall BC
            
                RC=30.0 ! Clustering ratio, must be >= 1 (Ex: 10 for 10x density at center of plane)
                C1=1.0/(1.0+3.0/RC)
                C2=3.0/RC*C1
                dB=dsC/C2
                NumFSPx=ceiling(2.0*PlaneExtent/dB)
                NumFSPz=NumFSPz
                NumFSCPx=NumFSPx
                NumFSCPz=NumFSPz
                NumFSP=NumFSPx*NumFSPz
                NumFSCP=NumFSP ! colocation on panel
                
                ! Allocate local and global arrays
                allocate(xPanFS(NumFSPx+1),zPanFS(NumFSPx+1))
                Call wallsoln_fs_cns()
                
                dB=2.0/real(NumFSPx)
                do i=1,NumFSPx+1
                        B=-1.0+dB*(i-1)
                        xPanFS(i)=PlaneExtent*(C1*B**3+C2*B)
                end do
                zPanFS=xPanFS
                yPanFS=FSy
                
                ! Ordered list of panels going through x, then stepping y, repeat
                Ind=1
                do i=1,NumFSPz
                        do j=1,NumFSPx
                                P1=[xPanFS(j),yPanFS,zPanFS(i)]           ! lower panel x, lower panel y corner point
                                P2=[xPanFS(j+1),yPanFS,zPanFS(i)]         ! pos panel x neighbor point
                                P3=[xPanFS(j),yPanFS,zPanFS(i+1)] ! pos panel y neighbor point
                                P4=[xPanFS(j+1),yPanFS,zPanFS(i+1)]       ! pos panel x, pos panel y neighbor point
                                FSCPoints(Ind,1:3)=0.25*(P1+P2+P3+P4)        ! panel center
                                FSCPPoints(Ind,1:3)=FSCPoints(Ind,1:3)        ! colocation points (on panel)
                                FSXVec(Ind,1:3)=P2-P1                ! panel x tangential vector
                                FSPL(Ind)=sqrt(sum(FSXVec(Ind,1:3)**2))               ! panel x length
                                FSXVec(Ind,1:3)=FSXVec(Ind,1:3)/FSPL(Ind)              ! normalize
                                FSCXVec(Ind,1:3)=FSXVec(Ind,1:3)   ! colocation on panel
                                FSYVec(Ind,1:3)=P3-P1                ! panel y tangential vector, set so that panel normal will be in the domain inward direction
                                FSPW(Ind)=sqrt(sum(FSYVec(Ind,1:3)**2))               ! panel y length
                                FSYVec(Ind,1:3)=FSYVec(Ind,1:3)/FSPW(Ind)              ! normalize
                                FSCYVec(Ind,1:3)=FSYVec(Ind,1:3)   ! colocation on panel
                                Call cross(FSXVec(Ind,1),FSXVec(Ind,2),FSXVec(Ind,3),FSYVec(Ind,1),FSYVec(Ind,2),FSYVec(Ind,3),FSZVec(Ind,1),FSZVec(Ind,2),FSZVec(Ind,3))    ! panel normal vector   
                                FSZVec(Ind,1:3)=FSZVec(Ind,1:3)/sqrt(sum(FSZVec(Ind,1:3)**2))          ! normalize    
                                FSCZVec(Ind,1:3)=FSZVec(Ind,1:3)   ! colocation on panel    
                                Ind=Ind+1
                        end do
                end do
    
            else
            
                ! Use free surface BC
            
                ! Clustered panels
                RU=dsOut/dsC ! clustering factor at center, must be >= 1 (Ex: 10 for 10x density ratio at center)
                C1U=1.0/(1.0+3.0/RU)
                C2U=3.0/RU*C1U
                dBU=dsC/C2U
                NU=ceiling(PlaneExtent/dBU)
                
                RD=1.0 ! dont expand grid downstream of turbine as this occasionally causes instability...
                C1D=1.0/(1.0+3.0/RD)
                C2D=3.0/RD*C1D;
                dBD=dsC/C2D
                ND=ceiling(PlaneExtent/dBD)
                
                RW=dsOut/dsC
                C1W=1.0/(1.0+3.0/RW)
                C2W=3.0/RW*C1W
                dBW=dsC/C2W
                NW=ceiling(2.0*PlaneExtent/dBW)
            
                ! Number of panels and colocation points
                NumFSCPx=NU+ND
                NumFSCPz=NW
                NumFSPx=NumFSCPx+1 ! add extra panel for BC application
                NumFSPz=NumFSCPz
                NumFSP=NumFSPz*NumFSPx
                NumFSCP=NumFSCPz*NumFSCPx         
            
                ! Allocate local and global arrays
                allocate(xPanFS(NumFSPx+1),zPanFS(NumFSPz+1))
                Call wallsoln_fs_cns()
                
                dB=1.0/real(NU)
                do i=1,NU+1
                        B=-1.0+dB*(i-1)
                        xPanFS(i)=PlaneExtent*(C1U*B**3+C2U*B)
                end do
                
                dB=1.0/real(ND)
                do i=1,ND
                        B=0.0+dB*i
                        xPanFS(i+NU+1)=PlaneExtent*(C1D*B**3+C2D*B)
                end do
                
                dB=2.0/real(NW)
                do i=1,NW+1
                        B=-1.0+dB*(i-1)
                        zPanFS(i)=PlaneExtent*(C1W*B**3+C2W*B)
                end do
                
                yPanFS=FSy
   
                ! Offset actual panels a distance above domain avoid evaluating
                ! the tangential velocity close to the panel, especially near the edges
                ! of panels where singularities occur.
                yPanFSOff=FSy+2.0*dsC
                
                ! Add extra panel in x direction for BC application
                xPanFS(NumFSPx+1)=xPanFS(NumFSPx)+(xPanFS(NumFSPx)-xPanFS(NumFSPx-1))
                
                ! Ordered list of panels going through x, then stepping y, repeat
                Ind=1
                do i=1,NumFSPz
                        do j=1,NumFSPx
                                P1=[xPanFS(j),yPanFSOff,zPanFS(i)]           ! lower panel x, lower panel y corner point
                                P2=[xPanFS(j+1),yPanFSOff,zPanFS(i)]         ! pos panel x neighbor point
                                P3=[xPanFS(j),yPanFSOff,zPanFS(i+1)] ! pos panel y neighbor point
                                P4=[xPanFS(j+1),yPanFSOff,zPanFS(i+1)]       ! pos panel x, pos panel y neighbor point
                                FSCPoints(Ind,1:3)=0.25*(P1+P2+P3+P4)        ! panel center
                                FSXVec(Ind,1:3)=P2-P1                ! panel x tangential vector
                                FSPL(Ind)=sqrt(sum(FSXVec(Ind,1:3)**2))               ! panel x length
                                FSXVec(Ind,1:3)=FSXVec(Ind,1:3)/FSPL(Ind)              ! normalize
                                FSYVec(Ind,1:3)=P3-P1                ! panel y tangential vector, set so that panel normal will be in the domain inward direction
                                FSPW(Ind)=sqrt(sum(FSYVec(Ind,1:3)**2))               ! panel y length
                                FSYVec(Ind,1:3)=FSYVec(Ind,1:3)/FSPW(Ind)              ! normalize
                                Call cross(FSXVec(Ind,1),FSXVec(Ind,2),FSXVec(Ind,3),FSYVec(Ind,1),FSYVec(Ind,2),FSYVec(Ind,3),FSZVec(Ind,1),FSZVec(Ind,2),FSZVec(Ind,3))    ! panel normal vector   
                                FSZVec(Ind,1:3)=FSZVec(Ind,1:3)/sqrt(sum(FSZVec(Ind,1:3)**2))          ! normalize    
                                Ind=Ind+1
                        end do
                end do     
 
                ! Ordered list of colocation points going through x, then stepping y, repeat
                Ind=1
                do i=1,NumFSCPz
                        do j=1,NumFSCPx
                                P1=[xPanFS(j),yPanFS,zPanFS(i)]           ! lower panel x, lower panel y corner point
                                P2=[xPanFS(j+1),yPanFS,zPanFS(i)]         ! pos panel x neighbor point
                                P3=[xPanFS(j),yPanFS,zPanFS(i+1)] ! pos panel y neighbor point
                                P4=[xPanFS(j+1),yPanFS,zPanFS(i+1)]       ! pos panel x, pos panel y neighbor point
                                FSCPPoints(Ind,1:3)=0.25*(P1+P2+P3+P4)        ! colocation points (on panel)
                                FSCXVec(Ind,1:3)=(P2-P1)/sqrt(sum((P2-P1)**2))                 ! panel x tangential vector
                                FSCYVec(Ind,1:3)=(P3-P1)/sqrt(sum((P1-P3)**2))                ! panel y tangential vector, set so that panel normal will be in the domain inward direction
                                Call cross(FSCXVec(Ind,1),FSCXVec(Ind,2),FSCXVec(Ind,3),FSCYVec(Ind,1),FSCYVec(Ind,2),FSCYVec(Ind,3),FSCZVec(Ind,1),FSCZVec(Ind,2),FSCZVec(Ind,3))    ! panel normal vector   
                                FSCZVec(Ind,1:3)=FSCZVec(Ind,1:3)/sqrt(sum(FSCZVec(Ind,1:3)**2))          ! normalize    
                                Ind=Ind+1
                        end do
                end do                             
    
            end if
    
            ! Panel edge tolerance (needs to be less than 1/2 of the min panel dimension)
            FSEdgeTol=min(minval(FSPL),minval(FSPW))/10.0
        
        end if
	
return
end
