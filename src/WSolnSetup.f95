subroutine WSolnSetup()

    use util
    use wallsoln
    use wallsystem
!$  use omp_lib

    integer :: i, j, Self, IBCInd, BCRow
    integer :: INFO
    integer, allocatable :: IPIV(:)
    real :: R(3,3), Point(3), dPG(3), dVel(3), dVelG(3), dudx

    real :: t0,t1 ! cpu time variables for matrix inversion

    ! Set up ground plane system
    if (GPFlag==1) then

        !! This assumes WGeomSetup() has already been called, and that the ground plane has been configured
        !  as a wall system of one wall.

        ! Setup wall self influence matrix 
        write(*,*) 'Generating wall influence matrix...'
        call cpu_time(t0)
!$      t0 = omp_get_wtime()
        call gen_influence_matrix()
        call cpu_time(t1)
!$      t1 = omp_get_wtime()
        print '("Time to generate influence matrix = ",f15.3," seconds.")',t1-t0

        ! Store wall solution matrix and inverse
        write(*,*) 'Inverting wall influence matrix...'
        call cpu_time(t0)
!$      t0 = omp_get_wtime()
        call invert_influence_matrix()
        call cpu_time(t1)
!$      t1 = omp_get_wtime()
        print '("Time to generate influence matrix = ",f15.3," seconds.")',t1-t0

    end if


    ! Set up generic wall system
    if (WPFlag==1) then

        !! This assumes that the wall geometry has already been loaded as a wall system of one wall.

        ! Setup wall self influence matrix 
        write(*,*) 'Generating wall influence matrix...'
        call cpu_time(t0)
!$      t0 = omp_get_wtime()
        call gen_influence_matrix()
        call cpu_time(t1)
!$      t1 = omp_get_wtime()
        print '("Time to generate influence matrix = ",f15.3," seconds.")',t1-t0

        ! Store wall solution matrix and inverse
        write(*,*) 'Inverting wall influence matrix...'
        call cpu_time(t0)
!$      t0 = omp_get_wtime()
        call invert_influence_matrix()
        call cpu_time(t1)
!$      t1 = omp_get_wtime()
        print '("Time to generate influence matrix = ",f15.3," seconds.")',t1-t0

    end if


    ! Set up free surface system
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
        if (INFO>0) then
            write(6,'(A)') 'Matrix inversion failed in WSolnSetup. Exiting...'
            stop
        end if

        ! Initialize source strengths and RHS to zero
        FSSource(:,:)=0.0
        FSRHS(:,:)=0.0
        FSRHSInd=1

    end if

    return
end subroutine WSolnSetup

