subroutine WGeomSetup() 

	use wallsoln  
	
	real :: PlaneExtent
	real :: P1(3,1), P2(3,1), P3(3,1), P4(3,1)
	real, allocatable :: xPan(:), zPan(:), B(:)
	integer :: Ind

	! Plane extent (over radius). To be applied in every direction around the turbine location...
	PlaneExtent=10.0 

	! Setup grid, clustered in the center of the plane
	RC=30.0 ! Clustering ratio, must be >= 1 (Ex: 10 for 10x density at center of plane)
	C1=1.0/(1.0+3.0/RC)
	C2=3.0/RC*C1
	dsC=0.1	! Panel dimension over radius in the center
	dB=dsC/C2
	NumWPx=ceiling(2.0*PlaneExtent/dB)
	NumWP=NumWPx**2
	
	! Allocate local and global arrays
	allocate(xPan(NumWPx+1),zPan(NumWPx+1),B(NumWPx+1))
	Call wallsoln_cns()
	
	dB=2.0/real(NumWPx)
	do i=0,NumWPx
		B(i+1)=-1.0+dB*i
	end do
	xPan=PlaneExtent*(C1*B**3+C2*B)
	zPan=xPan
	yPan=GPy
	
	! Ordered list of panels going through x, then stepping y, repeat
	Ind=1
	do i=1,NumWPx
		do j=1,NumWPx
			P1=reshape([xPan(j),yPan,zPan(i)],[3,1])   	! lower panel x, lower panel y corner point
			P2=reshape([xPan(j+1),yPan,zPan(i)],[3,1]) 	! pos panel x neighbor point
			P3=reshape([xPan(j),yPan,zPan(i+1)],[3,1]) 	! pos panel y neighbor point
			P4=reshape([xPan(j+1),yPan,zPan(i+1)],[3,1])  	! pos panel x, pos panel y neighbor point
			WCPoints(1:3,Ind)=0.25*(P1(1:3,1)+P2(1:3,1)+P3(1:3,1)+P4(1:3,1)) 	! panel center
			WXVec(1:3,Ind)=P2(1:3,1)-P1(1:3,1)		! panel x tangential vector
			WPL(Ind)=sqrt(sum(WXVec(1:3,Ind)**2))		! panel x length
			WXVec(1:3,Ind)=WXVec(1:3,Ind)/WPL(Ind)		! normalize
			WYVec(1:3,Ind)=P1(1:3,1)-P3(1:3,1)		! panel y tangential vector, set so that panel normal will be in the domain inward direction
			WPW(Ind)=sqrt(sum(WYVec(1:3,Ind)**2))		! panel y length
			WYVec(1:3,Ind)=WYVec(1:3,Ind)/WPW(Ind)		! normalize
			Call cross(WXVec(1,Ind),WXVec(2,Ind),WXVec(3,Ind),WYVec(1,Ind),WYVec(2,Ind),WYVec(3,Ind),WZVec(1,Ind),WZVec(2,Ind),WZVec(3,Ind))	! panel normal vector	
			WZVec(1:3,Ind)=WZVec(1:3,Ind)/sqrt(sum(WZVec(1:3,Ind)**2))		! normalize	
			Ind=Ind+1
		end do
	end do

	! Panel edge tolerance (needs to be less than 1/2 of the min panel dimension)
	WEdgeTol=min(minval(WPL),minval(WPW))/10.0
	
return
end
