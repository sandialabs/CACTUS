MODULE wallsoln
	
	! Wall panel geometry and solution arrays
	! Panel functions included...

	integer :: GPFlag			! Set to 1 to do ground plane calculation, otherwise 0

	real, allocatable :: WCPoints(:,:)	! Panel center points (over radius)
	real, allocatable :: WXVec(:,:)		! Panel tangential vectors in the length direction
	real, allocatable :: WYVec(:,:)		! Panel tangential vectors in the width
	real, allocatable :: WZVec(:,:)		! Panel normal vectors
	real, allocatable :: WPL(:)		! Panel lengths (over radius)
	real, allocatable :: WPW(:)		! Panel widths (over radius)
	real :: GPy				! y location of ground plane (over radius)
	real :: WEdgeTol			! Tolerance around panel edge in which to evaluate influence in special way (to avoid inf...)
	integer :: NumWPx			! Number of wall panels in the x direction
	integer :: NumWP			! Total number of wall panels
	
	real, allocatable :: WInCoeffN(:,:)	! Wall normal velocity self influence matrix
	real, allocatable :: WSource(:,:)       ! Wall source density values (column vector)
	real, allocatable :: WSMat(:,:) 	! Wall solution matrix
	real, allocatable :: WSMatI(:,:) 	! Inverse of the wall solution matrix
	real, allocatable :: WRHS(:,:) 		! Right hand side vector for the wall solution
	real, allocatable :: WWakeNVel(:) 	! Storage for wake induced wall normal velocity
	
	CONTAINS
	
	SUBROUTINE wallsoln_cns()
		
		! Constructor for the arrays in this module
		
		allocate(WCPoints(3,NumWP))
		allocate(WXVec(3,NumWP))
		allocate(WYVec(3,NumWP))
		allocate(WZVec(3,NumWP))
		allocate(WPL(NumWP))
		allocate(WPW(NumWP))
		
		allocate(WInCoeffN(NumWP,NumWP))
		allocate(WSource(NumWP,1))
		allocate(WRHS(NumWP,1))
		allocate(WWakeNVel(NumWP))
		allocate(WSMat(NumWP,NumWP))
		allocate(WSMatI(NumWP,NumWP))
	
	End SUBROUTINE
	
	SUBROUTINE RectSourceVel(PointI,L,Wth,Source,SelfInfluence,EdgeTol,CalcDer,Vel,dudx)
	
		! Calculate velocity induced by a rectangular source panel
		
		! Point is location in panel coord from center of panel
		! L is panel length, Wth is width
		! Source is panel source strength density
		! SelfInfluence is 1 when looking for influence of a panel on its own midpoint (at top of panel).
		! EdgeTol is the tolerance around edge in which to apply limiting conditions
		
		! Vel is the velocity in panel coordinates
		! dudx is the derivative of the panel x velocity in the panel x direction
		! (for application of the linear free surface method), only calculated if CalcDer is 1
		
		real :: PointI(3,1), L, Wth, Source, EdgeTol
		integer :: SelfInfluence, CalcDer 
		
		real :: Vel(3,1), dudx
		
		real :: pi, R2, R2P, dP1(3), dP2(3), dP3(3), dP4(3), Point(3)
		real :: u, v, w, sZ, Rp1, Rp2, Rp3, Rp4, h1, h2, h3, h4, R, A
		integer :: Flag
		
		! Define pi
		pi = 4.0*atan(1.0)
		
		! Make point rank 1 internally for code simplicity...
		Point(1)=PointI(1,1)
		Point(2)=PointI(2,1)
		Point(3)=PointI(3,1)
		
		R2=sum(Point*Point)
		R2P=((L/2.0+EdgeTol)**2+(Wth/2.0+EdgeTol)**2+EdgeTol**2)
		
		if (SelfInfluence==1) then
			! Self
			Vel=reshape([0.0,0.0,Source/2.0],[3,1])
			if (CalcDer==1) then
				dudx=2.0*sqrt(1.0/(1.0+(L/Wth)**2))*Source/(pi*L)
			else
				dudx=0.0
			end if
		else if (R2 < R2P) then
			! Near-field (check edge conditions)
			dP1=Point+[L/2.0,Wth/2.0,0.0]
			dP2=Point+[L/2.0,-Wth/2.0,0.0]
			dP3=Point+[-L/2.0,-Wth/2.0,0.0]
			dP4=Point+[-L/2.0,Wth/2.0,0.0]
			
			! Check edges
			Flag=0
			if (dP1(1)>(-EdgeTol) .AND. dP1(1)<(L+EdgeTol) .AND. abs(dP1(2))<EdgeTol .AND. abs(dP1(3))<EdgeTol) then
				Flag=Flag+1
			end if
			if (dP2(2)>(-Wth-EdgeTol) .AND. dP2(2)<(EdgeTol) .AND. abs(dP2(1))<EdgeTol .AND. abs(dP2(3))<EdgeTol) then
				Flag=Flag+1
			end if
			if (dP3(1)>(-L-EdgeTol) .AND. dP3(1)<(EdgeTol) .AND. abs(dP3(2))<EdgeTol .AND. abs(dP3(3))<EdgeTol) then
				Flag=Flag+1
			end if
			if (dP4(2)>(-EdgeTol) .AND. dP4(2)<(Wth+EdgeTol) .AND. abs(dP4(1))<EdgeTol .AND. abs(dP4(3))<EdgeTol) then
				Flag=Flag+1
			end if
			
			sZ=sign(1.0,Point(3))
			if (Flag==2) then
				! Tolerance may overlap with three other panels. Set to average
				! panel normal velocity, sum(Source/2)/4 (average over all 4 panels)
				Vel=sZ*reshape([0.0,0.0,Source/8.0],[3,1])
				dudx=0.0
			else if (Flag==1) then
				! Tolerance may overlap with one other panel. Set to average
				! panel normal velocity, sum(Source/2)/2 (average over both panels)
				Vel=sZ*reshape([0.0,0.0,Source/4.0],[3,1]) 
				dudx=0.0
			else
				! Full panel influence
				Rp1=sqrt(sum(dP1**2))
				Rp2=sqrt(sum(dP2**2))
				Rp3=sqrt(sum(dP3**2))
				Rp4=sqrt(sum(dP4**2))
				
				h1=dP1(1)*dP1(2)
				h2=dP2(1)*dP2(2)
				h3=dP3(1)*dP3(2)
				h4=dP4(1)*dP4(2)
				
				u=Source/(4.0*pi)*log(((Rp1+Rp2-Wth)*(Rp3+Rp4+Wth))/((Rp1+Rp2+Wth)*(Rp3+Rp4-Wth)))
				v=Source/(4.0*pi)*log(((Rp4+Rp1-L)*(Rp2+Rp3+L))/((Rp4+Rp1+L)*(Rp2+Rp3-L)))
				w=Source/(4.0*pi)*(atan(h1/(Point(3)*Rp1))+atan(h3/(Point(3)*Rp3))-atan(h2/(Point(3)*Rp2))-atan(h4/(Point(3)*Rp4)))
				
				Vel=reshape([u,v,w],[3,1]) 
				if (CalcDer==1) then
					dudx=Source/(2.0*pi)*Wth*((dP1(1)/Rp1+dP2(1)/Rp2)/((Rp1+Rp2-Wth)*(Rp1+Rp2+Wth)) - (dP3(1)/Rp3+dP4(1)/Rp4)/((Rp3+Rp4-Wth)*(Rp3+Rp4+Wth)))
				else 
					dudx=0.0
				end if
			end if
		
		else if (R2 < (6.0**2)*R2P) then
			! Mid-field (full panel influence)
			dP1=Point+[L/2.0,Wth/2.0,0.0]
			dP2=Point+[L/2.0,-Wth/2.0,0.0]
			dP3=Point+[-L/2.0,-Wth/2.0,0.0]
			dP4=Point+[-L/2.0,Wth/2.0,0.0]
			
			Rp1=sqrt(sum(dP1**2))
			Rp2=sqrt(sum(dP2**2))
			Rp3=sqrt(sum(dP3**2))
			Rp4=sqrt(sum(dP4**2))
			
			h1=dP1(1)*dP1(2)
			h2=dP2(1)*dP2(2)
			h3=dP3(1)*dP3(2)
			h4=dP4(1)*dP4(2)
			
			u=Source/(4.0*pi)*log(((Rp1+Rp2-Wth)*(Rp3+Rp4+Wth))/((Rp1+Rp2+Wth)*(Rp3+Rp4-Wth)))
			v=Source/(4.0*pi)*log(((Rp4+Rp1-L)*(Rp2+Rp3+L))/((Rp4+Rp1+L)*(Rp2+Rp3-L)))
			w=Source/(4.0*pi)*(atan(h1/(Point(3)*Rp1))+atan(h3/(Point(3)*Rp3))-atan(h2/(Point(3)*Rp2))-atan(h4/(Point(3)*Rp4)))
			
			Vel=reshape([u,v,w],[3,1]) 
			if (CalcDer==1) then
				dudx=Source/(2.0*pi)*Wth*((dP1(1)/Rp1+dP2(1)/Rp2)/((Rp1+Rp2-Wth)*(Rp1+Rp2+Wth)) - (dP3(1)/Rp3+dP4(1)/Rp4)/((Rp3+Rp4-Wth)*(Rp3+Rp4+Wth)))
			else
				dudx=0.0
			end if 
		else
			! Far-field (point source influence), at greater than 6 panel radii
			R=sqrt(R2)
			A=L*Wth
			Vel=reshape(Source*A/(4.0*pi)*Point/R**3,[3,1])
			if (CalcDer==1) then
				dudx=Source*A/(4.0*pi*R**3)*(1-3*Point(1)**2/R2)
			else
				dudx=0.0
			end if
 		end if
	
	End SUBROUTINE
	
	SUBROUTINE WallIndVel(PointG,Vel)
	
		real :: PointG(3,1), Vel(3,1)
		
		integer :: i
		real :: R(3,3), Point(3,1), dVel(3,1), dudx
	
		! Calculate velocity induced by all wall panels being used in the calculation
		Vel(:,:)=0.0
		if (GPFlag == 1) then
		
			do i=1,NumWP
			
				! Rotation from panel i to global
				R(1:3,1)=WXVec(1:3,i)
				R(1:3,2)=WYVec(1:3,i)
				R(1:3,3)=WZVec(1:3,i)
				
				! Calc influence in panel frame
				Point(1:3,1)=matmul(transpose(R),(PointG(1:3,1)-WCPoints(1:3,i)))
				Call RectSourceVel(Point,WPL(i),WPW(i),WSource(i,1),0,WEdgeTol,0,dVel,dudx)
				
				! Rotate to global frame
				dVel=matmul(R,dVel)
				Vel=Vel+dVel
			end do
			
		end if
	
	End SUBROUTINE
	
End