subroutine WSolnSetup()
    use wallsoln 

    integer :: i, j, Self, IBCInd, BCRow
    integer :: INFO
    integer, allocatable :: IPIV(:)
    real :: R(3,3), Point(3), dPG(3), dVel(3), dVelG(3), dudx


    if (GPFlag==1) then

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

        ! LAPACK => DGESV: Linear equation solution A*X=B where A(N,N) X(N,NRHS) B(N,NRHS)
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

        allocate(IPIV(NumWP)) ! allocation storage for pivot array
        Call DGESV(NumWP,NumWP,WSMat,NumWP,IPIV,WSMatI,NumWP,INFO)
        deallocate(IPIV)
        if (INFO>0) then
            write(6,*) 'Matrix inversion failed in WSolnSetup. Exiting...'
            stop
        end if

        ! Initialize source strengths and RHS to zero
        WSource(:,:)=0.0
        WRHS(:,:)=0.0

    end if


    if (FSFlag==1) then

        ! Setup free surface self influence matrix 
        do i=1,NumFSP
            do j=1,NumFSCP
                ! Rotation from global to panel i
                R(1,1:3)=FSXVec(i,1:3)
                R(2,1:3)=FSYVec(i,1:3)
                R(3,1:3)=FSZVec(i,1:3)

                ! Calc influence in panel frame
                dPG=FSCPPoints(j,1:3)-FSCPoints(i,1:3)                         
                Call CalcRotation3(R,dPG,Point,0)                     
                Call RectSourceVel(Point,FSPL(i),FSPW(i),1.0,0,FSEdgeTol,1,dVel,dudx)

                ! Rotate to global frame
                Call CalcRotation3(R,dVel,dVelG,1)                      
                FSInCoeffN(j,i)=sum(dVelG*FSCZVec(j,1:3))
                FSInCoeffT(j,i)=sum(dVelG*FSCXVec(j,1:3))
                FSInCoeffdUdX(j,i)=dudx
            end do
        end do

        ! Create free surface solution matrix and inverse
        if (UseFSWall) then 
            ! Wall solution matrix
            FSSMat=FSInCoeffN 
        else
            ! Free surface solution matrix
            FSSMat(1:NumFSCP,1:NumFSP)=FSInCoeffN-(FnR**2)*FSInCoeffdUdX
            ! Add inflow BCs
            FSBCRow(:)=0
            do i=1,NumFSCPz
                IBCInd=1+(i-1)*NumFSCPx
                BCRow=NumFSCP+i
                FSSMat(BCRow,1:NumFSP)=FSInCoeffdUdX(IBCInd,1:NumFSP)
                ! Set BC row index
                FSBCRow(IBCInd)=BCRow
            end do
        end if

        ! LAPACK => DGESV: Linear equation solution A*X=B where A(N,N) X(N,NRHS) B(N,NRHS)
        ! Note that if NRHS = N, and B is the identity, X is the inverse of A...
        ! Initialize inverse to the identity
        FSSMatI(:,:)=0.0
        do i=1,NumFSP
            do j=1,NumFSP
                if (j==i) then
                    FSSMatI(i,j)=1.0	
                end if
            end do
        end do

        allocate(IPIV(NumFSP)) ! allocation storage for pivot array
        Call DGESV(NumFSP,NumFSP,FSSMat,NumFSP,IPIV,FSSMatI,NumFSP,INFO)
        deallocate(IPIV)
        if (INFO>0) then
            write(6,*) 'Matrix inversion failed in WSolnSetup. Exiting...'
            stop
        end if

        ! Initialize source strengths and RHS to zero
        FSSource(:,:)=0.0
        FSRHS(:,:)=0.0
        FSRHSInd=1

    end if

    return
end subroutine WSolnSetup

