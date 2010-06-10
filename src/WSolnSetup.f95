subroutine WSolnSetup()

	use wallsoln 
	
	integer :: i, j, Self
	integer :: INFO
	integer :: IPIV(NumWP)
	real :: R(3,3), Point(3,1), dVel(3,1), dudx
	
	! Setup wall self influence matrix 
	do i=1,NumWP
		do j=1,NumWP
			if (j==i) then
				Self=1
			else
				Self=0
			end if
			
			! Rotation from panel i to global
			R(1:3,1)=WXVec(1:3,i)
			R(1:3,2)=WYVec(1:3,i)
			R(1:3,3)=WZVec(1:3,i)
			
			! Calc influence in panel frame
			Point(1:3,1)=matmul(transpose(R),(WCPoints(1:3,j)-WCPoints(1:3,i)))
			Call RectSourceVel(Point,WPL(i),WPW(i),1.0,Self,WEdgeTol,0,dVel,dudx)
			
			! Rotate to global frame
			dVel=matmul(R,dVel)
			WInCoeffN(j,i)=sum(dVel(1:3,1)*WZVec(1:3,j))
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

