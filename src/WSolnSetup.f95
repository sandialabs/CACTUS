subroutine WSolnSetup()

	use wallsoln 
	
	integer :: i, j, Self
	integer :: INFO
	integer :: IPIV(NumWP)
	real :: R(3,3), Point(3), dPG(3), dVel(3), dVelG(3), dudx
       
	
	! Setup wall self influence matrix 
	do i=1,NumWP
		do j=1,NumWP
			if (j==i) then
				Self=1
			else
				Self=0
			end if
			
                        ! Rotation from global to panel i
                        R(1,1:3)=WXVec(i,1:3)
                        R(2,1:3)=WYVec(i,1:3)
                        R(3,1:3)=WZVec(i,1:3)
			
			! Calc influence in panel frame
                        dPG=WCPoints(j,1:3)-WCPoints(i,1:3)                         
                        Call CalcRotation3(R,dPG,Point,0)                     
			Call RectSourceVel(Point,WPL(i),WPW(i),1.0,Self,WEdgeTol,0,dVel,dudx)
                        
			! Rotate to global frame
                        Call CalcRotation3(R,dVel,dVelG,1)                      
			WInCoeffN(j,i)=sum(dVelG*WZVec(j,1:3))
		end do
	end do
	
	! Store wall solution matrix and inverse
	WSMat=WInCoeffN
	! LAPACK => SGESV: Linear equation solution A*X=B where A(N,N) X(N,NRHS) B(N,NRHS)
	! Note that if NRHS = N, and B is the identity, X is the inverse of A...
	! Initialize inverse to the identity
	WSMatI(:,:)=0.0
	do i=1,NumWP
		do j=1,NumWP
			if (j==i) then
				WSMatI(i,j)=1.0	
			end if
		end do
	end do
	
	Call SGESV(NumWP,NumWP,WSMat,NumWP,IPIV,WSMatI,NumWP,INFO)
	
	! Initialize source strengths to zero
	WSource(:,:)=0.0
	! Initialize wake induced normal velocity to zero
	WWakeNVel(:)=0.0

return
end

