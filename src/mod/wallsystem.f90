module wallsystem

    ! wallsystem module for managing multiple walls and wall solution (strengths)

    use util, only: cross
    use vecutils, only: cross3, mag3, calcrotation3
    use plot3d, only: read_p3d_multiblock
    use wallgeom
    use quadsourcepanel
    use ieee_arithmetic, only: ieee_is_nan

    implicit none

    ! wallsystem geometry
    integer           :: Nwalls              ! number of walls
    integer           :: NumWP_total         ! total number of wall panels
    integer           :: NumWPNodes_total    ! total number of nodes

    ! concatenated source panel geometry
    real, allocatable :: W1Vec(:,:)
    real, allocatable :: W2Vec(:,:)
    real, allocatable :: W3Vec(:,:)
    real, allocatable :: WCPoints(:,:)

    ! wallsystem solution
    real, allocatable :: WInf(:,:)                     ! global influence matrix
    real, allocatable :: WInfI(:,:)                    ! inverse of the global influence matrix
    real, allocatable :: WRHS(:,:)                     ! global right hand side vector for wall solution
    real, allocatable :: WSource(:,:)                  ! global wall strengths

    ! array to hold Walls
    type(WallType), allocatable :: Walls(:)

contains

    subroutine wallsystem_cns()

        ! wallsystem_cns() : allocate storage for the wallsystem

        integer :: iw, vp

        ! set the number of walls
        Nwalls = size(Walls)

        ! count the number of panels and nodes
        NumWP_total = 0
        NumWPNodes_total = 0
        do iw=1,Nwalls
            NumWP_total = NumWP_total + Walls(iw)%NumWP
            NumWPNodes_total = NumWPNodes_total + Walls(iw)%NumWPNodes
        end do

        ! allocate storage for the concatenated geometry arrays
        allocate(W1Vec(NumWP_total,3))
        allocate(W2Vec(NumWP_total,3))
        allocate(W3Vec(NumWP_total,3))
        allocate(WCPoints(NumWP_total,3))

        ! allocate storage for the solution arrays
        allocate(WInf(NumWP_total,NumWP_total))
        allocate(WInfI(NumWP_total,NumWP_total))
        allocate(WRHS(NumWP_total,1))
        allocate(WSource(NumWP_total,1))

        ! concatenate wall geometry into single long arrays
        vp = 1
        do iw=1,NWalls
            W1Vec(   vp:vp+Walls(iw)%NumWP-1,:) = Walls(iw)%W1Vec
            W2Vec(   vp:vp+Walls(iw)%NumWP-1,:) = Walls(iw)%W2Vec
            W3Vec(   vp:vp+Walls(iw)%NumWP-1,:) = Walls(iw)%W3Vec
            WCPoints(vp:vp+Walls(iw)%NumWP-1,:) = Walls(iw)%WCPoints

            vp = vp + Walls(iw)%NumWP
        end do

        ! initialize source panel strengths and RHS to zero
        WSource(:,:) = 0.0
        WRHS(:,:)    = 0.0

    end subroutine wallsystem_cns


    subroutine gen_influence_matrix()

        ! gen_influence_matrix() : Generates the influence matrix WInf(:,:) for the wall system.

        real    :: p(3), p_center(3)
        real    :: p1(3), p2(3), p3(3), p4(3)
        real    :: p_plane(3), p1_plane(3), p2_plane(3), p3_plane(3), p4_plane(3)

        real    :: vel(3), temp
        real    :: vel_global(3)
        real    :: R(3,3)

        ! temporary variables cont'd: counters for looping through the contatenation
        integer :: ip_global,jp_global      ! inner/outer wall panel index (in entire wallsystem)
        integer :: info
        integer :: selfinfluence

        ! compute the influence matrix
!$omp parallel do private(jp_global,selfinfluence,R,p1,p2,p3,p4,p,p_center,p_plane,p1_plane,p2_plane,p3_plane,p4_plane,vel,temp,info,vel_global)
        do ip_global=1,NumWP_total
            do jp_global=1,NumWP_total

                ! check if we are computing induced velocity of panel on itself
                if (jp_global==ip_global) then
                    selfinfluence = 1
                else
                    selfinfluence = 0
                end if

                ! rotation matrix from global to panel jp
                R(1,1:3)=W1Vec(jp_global,1:3)
                R(2,1:3)=W2Vec(jp_global,1:3)
                R(3,1:3)=W3Vec(jp_global,1:3)

                ! get the coordinates of the nodes which make up the jp-th panel
                call ip_global_to_nodepoints(jp_global,p1,p2,p3,p4)

                ! get the relative position of the panel centers (outer-inner)
                p        = WCPoints(ip_global,1:3) ! the point of interest
                p_center = WCPoints(jp_global,1:3)

                ! rotate the point locations to the local panel coordinate system
                Call CalcRotation3(R,p- p_center,p_plane ,0)
                Call CalcRotation3(R,p1-p_center,p1_plane,0)
                Call CalcRotation3(R,p2-p_center,p2_plane,0)
                Call CalcRotation3(R,p3-p_center,p3_plane,0)
                Call CalcRotation3(R,p4-p_center,p4_plane,0)

                Call quadsourcevel(p_plane,p1_plane,p2_plane,p3_plane,p4_plane,1.0,selfinfluence,0,vel,temp,info)

                ! Rotate back to global frame
                Call CalcRotation3(R,vel,vel_global,1)

                ! Calc dot product of induced velocity with the point of interest panel normal
                ! WInf(jp_global,ip_global) = sum(vel_global*W3Vec(ip_global,1:3))
                WInf(ip_global,jp_global) = sum(vel_global*W3Vec(ip_global,1:3))

            end do
        end do
!$omp end parallel do

    end subroutine gen_influence_matrix


    subroutine invert_influence_matrix()

        ! invert_influence_matrix() : inverts the influence matrix using LAPACK DGESV
        !   WInf^-1 is stored in WInfI
        real, allocatable    :: IPIV(:)
        integer :: info

        integer :: i, j

        ! allocate memory
        allocate(IPIV(NumWP_total))    ! pivot matrix

        ! LAPACK => DGESV: Linear equation solution A*X=B where A(N,N) X(N,NRHS) B(N,NRHS)
        ! Note that if NRHS = N, and B is the identity, X is the inverse of A...
        ! Initialize inverse to the identity
        WInfI(:,:) = 0.0
        do i=1,NumWP_total
            do j=1,NumWP_total
                if (j==i) then
                    WInfI(i,j) = 1.0
                end if
            end do
        end do

        ! solve A*X=I for X. The solution is written to WInfI
        call DGESV(NumWP_total,NumWP_total,WInf,NumWP_total,IPIV,WInfI,NumWP_total,info)

        if (info > 0) then
            write(6,'(A)') 'Matrix inversion failed in wallsystem.f95. Exiting.'
            stop
        end if

    end subroutine invert_influence_matrix


    subroutine wall_ind_vel(p,calcder,vel,dudx)

        ! wall_ind_vel() : calculates the velocity induced by all the wall panels at a point

        real, intent(in)    :: p(3)
        integer, intent(in) :: calcder

        real, intent(out)   :: vel(3), dudx

        ! temporary variables
        real                :: R(3,3)
        integer             :: i
        real                :: p_center(3)
        real                :: p1(3), p2(3), p3(3), p4(3)
        real                :: p_plane(3), p1_plane(3), p2_plane(3), p3_plane(3), p4_plane(3)
        real                :: dvel(3), dvel_global(3)
        real                :: sigma
        integer             :: info
        real                :: quad(4,3)

        ! initialize velocity to 0.0
        vel(:) = 0.0

        do i=1,NumWP_total
            ! get the node coordinates
            call ip_global_to_nodepoints(i,p1,p2,p3,p4)

            ! rotation matrix from global to panel coordinates
            R(1,1:3) = W1Vec(i,1:3)
            R(2,1:3) = W2Vec(i,1:3)
            R(3,1:3) = W3Vec(i,1:3)

            ! get the center of the ith panel
            p_center = WCPoints(i,1:3)

            ! rotate the point locations to the local panel plane
            ! (relative to the panel center)
            Call calcrotation3(R,p -p_center,p_plane ,0)
            Call calcrotation3(R,p1-p_center,p1_plane,0)
            Call calcrotation3(R,p2-p_center,p2_plane,0)
            Call calcrotation3(R,p3-p_center,p3_plane,0)
            Call calcrotation3(R,p4-p_center,p4_plane,0)

            ! compute the induced velocity (panel coordinates)
            sigma = WSource(i,1)

            ! compute the induced velocity dvel at the the point
            ! (selfinfluence=0, calcder=0)
            Call quadsourcevel(p_plane,p1_plane,p2_plane,p3_plane,p4_plane,sigma,0,0,dvel,dudx,info)

            ! Rotate back to global frame
            Call calcrotation3(R,dvel,dvel_global,1)

            ! add the velocity
            vel = vel + dvel_global

            ! check for problems in dvel
            if (ieee_is_nan(dvel(1)) .or. ieee_is_nan(dvel(2)) .or. ieee_is_nan(dvel(3))) then
                quad(1,:) = p1_plane(1:3)
                quad(2,:) = p2_plane(1:3)
                quad(3,:) = p3_plane(1:3)
                quad(4,:) = p4_plane(1:3)

                write(*,*) "Warning: NaN encountered in wall_ind_vel()"
                ! stop
            end if
        end do

    end subroutine wall_ind_vel


    subroutine ip_global_to_iwip_local(ip_global,iw,ip_local)
        ! ip_global_to_iwip_local() : returns the wall index and local panel index of a given "global" wall panel

        integer, intent(in) :: ip_global
        integer, intent(out) :: iw, ip_local

        iw = 1
        ip_local = ip_global
        do while (ip_local - Walls(iw)%NumWP > 0)
            ip_local = ip_local - Walls(iw)%NumWP
            iw = iw + 1
        end do

    end subroutine ip_global_to_iwip_local


    subroutine ip_global_to_nodepoints(ip_global,p1,p2,p3,p4)

        ! ip_global_to_nodepoints() : returns the four node points p1,p2,p3,p4 corresponding to a given "global" wall panel

        integer, intent(in) :: ip_global
        real, intent(out)   :: p1(3), p2(3), p3(3), p4(3)

        integer             :: iw, ip_local
        integer             :: i_panel, j_panel

        ! get the wall number, iw
        call ip_global_to_iwip_local(ip_global,iw,ip_local)

        ! get the local panel indices, i_panel & j_panel
        call ip_local_to_ij_panel(Walls(iw),ip_local,i_panel,j_panel)

        ! get the four nodes which make up the quadrilateral
        call ij_local_to_nodepoints(Walls(iw),i_panel,j_panel,p1,p2,p3,p4)

    end subroutine ip_global_to_nodepoints


    subroutine read_p3d_walls(WallMeshPath)
        ! read_p3d_walls() : Read in a 2-D wall mesh from a Plot3d grid file.
        !   Grid data from each block is loaded into a type(Wall), an array of which is stored
        !   as Walls(:).
        !
        ! WallMeshPath (input) : Path to file containing a multi-block structured mesh.

        character(len=*), intent(in) :: WallMeshPath
        integer                      :: nblocks
        integer, allocatable         :: ni(:), nj(:), nk(:)
        real, allocatable            :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)
        integer                      :: m,i,j
        integer, parameter           :: k=1
        integer                      :: nnodes,npanels
        real                         :: p1(3), p2(3), p3(3), p4(3)
        integer                      :: idx

        ! read the plot3d mesh
        Call read_p3d_multiblock(WallMeshPath,nblocks,ni,nj,nk,x,y,z)

        ! set the module variable
        Nwalls = nblocks

        ! allocate for Walls
        allocate(Walls(nblocks))

        ! set initial value fo NumWP_total
        NumWP_total = 0
        NumWPNodes_total = 0

        ! for each block (or each "wall panel")
        do m = 1,nblocks
            if (nk(m) > 1) then
                write(*,*) "Error in read_p3d_walls(): Wall dimension in 3rd direction is > 1 -- mesh is a volume, not a wall!"
                stop
            end if

            ! compute the number of panels
            nnodes  = ni(m)*nj(m)
            npanels = (ni(m)-1)*(nj(m)-1)

            ! add to the total number of wall panels
            NumWP_total = NumWP_total + npanels
            NumWPNodes_total = NumWPNodes_total + nnodes

            ! set some WallType variables
            Walls(m)%NumWP1     = ni(m)-1
            Walls(m)%NumWP2     = nj(m)-1
            Walls(m)%NumWP      = npanels
            Walls(m)%NumWPNodes = nnodes

            ! call the wall constructor subroutine (allocator for internal wall variables)
            Call wall_cns(Walls(m))

            ! loop through grid coordinates
            idx = 1

            do j=1,nj(m)-1
                do i=1,ni(m)-1
                    ! get cell corners (ordered clockwise from lowest i,j when viewed from behind panel)
                    p1 = [x(i  ,j  ,k,m), y(i  ,j  ,k,m), z(i  ,j  ,k,m)]
                    p2 = [x(i+1,j  ,k,m), y(i+1,j  ,k,m), z(i+1,j  ,k,m)]
                    p3 = [x(i+1,j+1,k,m), y(i+1,j+1,k,m), z(i+1,j+1,k,m)]
                    p4 = [x(i  ,j+1,k,m), y(i  ,j+1,k,m), z(i  ,j+1,k,m)]

                    Walls(m)%WCPoints(idx,1:3) = 0.25*(p1+p2+p3+p4)                                ! panel center
                    Walls(m)%W1Vec(idx,1:3)    = (p2-p1)/mag3(p2-p1)                     ! panel 1-tangential vector
                    Walls(m)%W2Vec(idx,1:3)    = (p4-p1)/mag3(p4-p1)                     ! panel 2-tangential vector

                    ! compute panel normal vector
                    Call cross(Walls(m)%W1Vec(idx,1),Walls(m)%W1Vec(idx,2),Walls(m)%W1Vec(idx,3),&
                               Walls(m)%W2Vec(idx,1),Walls(m)%W2Vec(idx,2),Walls(m)%W2Vec(idx,3),&
                               Walls(m)%W3Vec(idx,1),Walls(m)%W3Vec(idx,2),Walls(m)%W3Vec(idx,3))

                    idx = idx + 1
                end do
            end do

            ! store node locations in an array (flattening)
            idx = 1
            do j=1,nj(m)
                do i=1,ni(m)
                    ! store node locations in array
                    Walls(m)%pnodes(idx,1) = x(i,j,k,m)
                    Walls(m)%pnodes(idx,2) = y(i,j,k,m)
                    Walls(m)%pnodes(idx,3) = z(i,j,k,m)
                    idx = idx + 1
                end do
            end do
        end do

    end subroutine read_p3d_walls

end module wallsystem