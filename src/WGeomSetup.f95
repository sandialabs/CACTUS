subroutine WGeomSetup() 

	use wallsoln  
	
	real :: PlaneExtent
	real :: P1(3), P2(3), P3(3), P4(3)
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
			P1=[xPan(j),yPan,zPan(i)]   	! lower panel x, lower panel y corner point
			P2=[xPan(j+1),yPan,zPan(i)] 	! pos panel x neighbor point
			P3=[xPan(j),yPan,zPan(i+1)]	! pos panel y neighbor point
			P4=[xPan(j+1),yPan,zPan(i+1)]  	! pos panel x, pos panel y neighbor point
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
	
return
end
