MODULE wallsoln

    ! Wall panel geometry and solution arrays
    ! Panel functions included...

    use wallsystem

    ! Ground plane (rectangular source panel implementation)

    integer :: GPFlag          ! Set to 1 to do ground plane calculation, otherwise 0

    real    :: GPy             ! y location of ground plane (over radius)
    real    :: GPGridSF        ! Grid scale factor (factor on default grid discretization level)       
    real    :: GPGridExtent    ! Distance ground plane extends in both directions from the origin
    real    :: WEdgeTol        ! Tolerance around panel edge in which to evaluate influence in special way (to avoid inf...)
    integer :: NumWP           ! Total number of wall panels    


    ! Wall System
    integer :: WallFlag        ! Set to 1 to do wall calculation, otherwise 0


    ! Free surface (rectangular source panel implementation)
    ! Note: this method assumes that both the nominal freestream is oriented in the positive x-direction and
    ! gravity is oriented in the negative y-direction. The free surface and ground plane panel systems must have
    ! their x tangent vectors aligned with the nominal freestream...

    integer :: FSFlag                        ! Set to 1 to do free surface calculation, otherwise 0

    real, allocatable :: FSCPoints(:,:)      ! Panel center points (over radius)
    real, allocatable :: FSCPPoints(:,:)     ! Colocation points (over radius)
    real, allocatable :: FSXVec(:,:)         ! Panel tangential vectors in the length direction
    real, allocatable :: FSYVec(:,:)         ! Panel tangential vectors in the width
    real, allocatable :: FSZVec(:,:)         ! Panel normal vectors
    real, allocatable :: FSPL(:)             ! Panel lengths (over radius)
    real, allocatable :: FSPW(:)             ! Panel widths (over radius)
    real, allocatable :: FSCXVec(:,:)        ! Colocation tangential vectors in the length direction
    real, allocatable :: FSCYVec(:,:)        ! Colocation tangential vectors in the width
    real, allocatable :: FSCZVec(:,:)        ! Colocation normal vectors
    real :: FSy                              ! y location of undeflected free surface (over radius)
    real :: FSGridSF                         ! Near field grid scale factor (factor on default near field grid discretization level)
    real :: FnR                              ! Froude number based on turbine radius, FnR = Uinf/sqrt(g*R)
    real :: FSEdgeTol                        ! Tolerance around panel edge in which to evaluate influence in special way (to avoid inf...)
    integer :: NumFSPx                       ! Number of panels in the x direction
    integer :: NumFSPz                       ! Number of panels in the z direction
    integer :: NumFSCPx                      ! Number of colocation points in the x direction
    integer :: NumFSCPz                      ! Number of colocation points in the z direction
    integer :: NumFSP                        ! Total number of wall panels
    integer :: NumFSCP                       ! Total number of colocation points
    integer :: NFSRHSAve                     ! Number of RHS evaluations to use in running average (should cover approx 1 revolution...)
    integer :: FSRHSInd                      ! RHS averaging array index
    logical :: UseFSWall                     ! Use wall BC on FS or not

    real, allocatable :: FSInCoeffN(:,:)     ! FS normal velocity self influence matrix
    real, allocatable :: FSInCoeffT(:,:)     ! FS tangent (along length) velocity self influence matrix
    real, allocatable :: FSInCoeffdUdX(:,:)  ! FS tangent (along length) velocity gradient self influence matrix
    real, allocatable :: FSSource(:,:)       ! FS source density values (column vector) (non-dimensional, normalized by freestream velocity)
    real, allocatable :: FSSMat(:,:)         ! FS solution matrix
    real, allocatable :: FSSMatI(:,:)        ! Inverse of the FS solution matrix
    real, allocatable :: FSRHS(:,:)          ! Right hand side vector for the FS solution (list for running average)
    real, allocatable :: FSRHSAve(:,:)       ! Average right hand side vector for the FS solution (over approx 1 revolution)
    integer, allocatable :: FSBCRow(:)       ! Row index of BC equation in solution matrix for each colocation point (if CP has BC, 0 otherwise)

    ! Free surface data output
    character(1000) :: FSOutHead = 'X/R (-),Y/R (-),Z/R (-),U/Uinf (-),dH/R (-)'    
    real, allocatable :: FSVT(:,:)          ! Right hand side vector for the FS solution (list for running average)
    real, allocatable :: FSVTAve(:,:)       ! Average right hand side vector for the FS solution (over approx 1 revolution)         

    ! Wall update interval
    integer :: iWall                         ! Wall update interval (number of timesteps)  

    ! Output function params       
    integer :: WallOutFlag                   ! Output wall panel data
    integer :: OutCount                      ! global counter for wall output function
    integer :: OutIterFS                     ! timestep iteration on final revolution on which to write FS output


CONTAINS

    SUBROUTINE wallsoln_fs_cns()

        ! Free surface                
        allocate(FSCPoints(NumFSP,3))
        allocate(FSXVec(NumFSP,3))
        allocate(FSYVec(NumFSP,3))
        allocate(FSZVec(NumFSP,3))
        allocate(FSPL(NumFSP))
        allocate(FSPW(NumFSP))
        allocate(FSCPPoints(NumFSCP,3))
        allocate(FSCXVec(NumFSCP,3))
        allocate(FSCYVec(NumFSCP,3))
        allocate(FSCZVec(NumFSCP,3))
        allocate(FSInCoeffN(NumFSCP,NumFSP))
        allocate(FSInCoeffT(NumFSCP,NumFSP))
        allocate(FSInCoeffdUdX(NumFSCP,NumFSP))
        allocate(FSSource(NumFSP,1))
        allocate(FSRHS(NumFSP,NFSRHSAve))
        allocate(FSRHSAve(NumFSP,1))
        allocate(FSSMat(NumFSP,NumFSP))
        allocate(FSSMatI(NumFSP,NumFSP))
        allocate(FSBCRow(NumFSCP))
        ! Output
        allocate(FSVT(NumFSCP,NFSRHSAve))
        allocate(FSVTAve(NumFSCP,1))

    End SUBROUTINE wallsoln_fs_cns


    SUBROUTINE RectSourceVel(Point,L,Wth,Source,SelfInfluence,EdgeTol,CalcDer,Vel,dudx)

        ! Calculate velocity induced by a rectangular source panel

        ! Point is location in panel coord from center of panel (row vector)
        ! L is panel length, Wth is width
        ! Source is panel source strength density
        ! SelfInfluence is 1 when looking for influence of a panel on its own midpoint (at top of panel).
        ! EdgeTol is the tolerance around edge in which to apply limiting conditions

        ! Vel is the velocity in panel coordinates (row vector)
        ! dudx is the derivative of the panel x velocity in the panel x direction
        ! (for application of the linear free surface method), only calculated if CalcDer is 1

        ! Note: the use of fortran 95 array math intrinsic functions (reshape, matmul) has been avoided to speed things up...

        real :: Point(3), L, Wth, Source, EdgeTol
        integer :: SelfInfluence, CalcDer 

        real :: Vel(3), dudx

        real :: pi, R2, R2P, dP1(3), dP2(3), dP3(3), dP4(3)
        real :: u, v, w, sZ, Rp1, Rp2, Rp3, Rp4, h1, h2, h3, h4, R, A
        integer :: Flag

        ! Define pi
        pi = 4.0*atan(1.0)

        R2=sum(Point**2)           
        R2P=((L/2.0+EdgeTol)**2+(Wth/2.0+EdgeTol)**2+EdgeTol**2)

        if (SelfInfluence==1) then
            ! Self
            Vel=[0.0,0.0,Source/2.0]
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
                Vel=sZ*[0.0,0.0,Source/8.0]
                dudx=0.0
            else if (Flag==1) then
                ! Tolerance may overlap with one other panel. Set to average
                ! panel normal velocity, sum(Source/2)/2 (average over both panels)
                Vel=sZ*[0.0,0.0,Source/4.0]
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

                Vel=[u,v,w]
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

            Vel=[u,v,w]
            if (CalcDer==1) then
                dudx=Source/(2.0*pi)*Wth*((dP1(1)/Rp1+dP2(1)/Rp2)/((Rp1+Rp2-Wth)*(Rp1+Rp2+Wth)) - (dP3(1)/Rp3+dP4(1)/Rp4)/((Rp3+Rp4-Wth)*(Rp3+Rp4+Wth)))
            else
                dudx=0.0
            end if
        else
            ! Far-field (point source influence), at greater than 6 panel radii
            R=sqrt(R2)
            A=L*Wth
            Vel=Source*A/(4.0*pi)*Point/R**3
            if (CalcDer==1) then
                dudx=Source*A/(4.0*pi*R**3)*(1-3*Point(1)**2/R2)
            else
                dudx=0.0
            end if
        end if

    End SUBROUTINE RectSourceVel


    SUBROUTINE GPIndVel(PointG,CalcDer,Vel,dudx)

        real :: PointG(3), Vel(3), dudx
        integer :: CalcDer

        ! Calculate velocity induced by ground plane panels.
        ! Use CalcDer=1 to calc dudx

        ! initialize to zero
        Vel(:)=0.0

        ! compute the velocity from all panels in wallsystem
        call wall_ind_vel(PointG,CalcDer,Vel,dudx)


    End SUBROUTINE GPIndVel


    SUBROUTINE WPIndVel(PointG,CalcDer,Vel,dudx)

        real :: PointG(3), Vel(3), dudx
        integer :: CalcDer

        ! Calculate velocity induced by ground plane panels.
        ! Use CalcDer=1 to calc dudx

        ! initialize to zero
        Vel(:)=0.0

        ! compute the velocity from all panels in wallsystem
        call wall_ind_vel(PointG,CalcDer,Vel,dudx)


    End SUBROUTINE WPIndVel


    SUBROUTINE FSIndVel(PointG,CalcDer,Vel,dudx)

        real :: PointG(3), Vel(3), dudx
        integer :: CalcDer

        integer :: i
        real :: R(3,3), Point(3), dPG(3), dVel(3), dVelG(3)

        ! Calculate velocity induced by free surface panels.
        ! Use CalcDer=1 to calc dudx
        ! Note: the use of fortran 95 array math intrinsic functions (reshape, matmul) has been avoided to speed things up...
        Vel(:)=0.0

        do i=1,NumFSP

            ! Rotation from global to panel i
            R(1,1:3)=FSXVec(i,1:3)
            R(2,1:3)=FSYVec(i,1:3)
            R(3,1:3)=FSZVec(i,1:3)

            ! Calc influence in panel frame
            dPG=PointG-FSCPoints(i,1:3)                            
            Call CalcRotation3(R,dPG,Point,0)                        
            Call RectSourceVel(Point,FSPL(i),FSPW(i),FSSource(i,1),0,FSEdgeTol,CalcDer,dVel,dudx)

            ! Rotate to global frame
            Call CalcRotation3(R,dVel,dVelG,1)                             

            Vel=Vel+dVelG

        end do

    End SUBROUTINE FSIndVel


    SUBROUTINE WallIndVel(PointG,Vel)

        ! Calculate velocity induced by all wall panels being used in the calculation
        ! Note: the use of fortran 95 array math intrinsic functions (reshape, matmul) has been avoided to speed things up...

        real, intent(in) :: PointG(3)
        real, intent(out) :: Vel(3)

        real :: dVel(3), dudx

        Vel(:)=0.0

        if (GPFlag == 1) then
            Call GPIndVel(PointG,0,dVel,dudx)
            Vel=Vel+dVel
        end if

        if (WallFlag == 1) then
            Call WPIndVel(PointG,0,dVel,dudx)
            Vel=Vel+dVel
        end if

        if (FSFlag == 1) then
            Call FSIndVel(PointG,0,dVel,dudx)
            Vel=Vel+dVel
        end if

    End SUBROUTINE WallIndVel

End MODULE wallsoln
