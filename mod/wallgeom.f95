module wallgeom

    ! Wall geometry and outputs module

    use util
    use plot3d

    implicit none

    type WallType

        ! wall geometry
        integer           :: NumWP1                   ! Number of wall panels in the 1-direction
        integer           :: NumWP2                   ! Number of wall panels in the 2-direction
        integer           :: NumWP                    ! Number of wall panels (on this particular wall)
        integer           :: NumNodes                 ! Number of nodes (on this particula wall)

        real, allocatable :: WCPoints(:,:)            ! Panel center points (over radius)
        real, allocatable :: W1Vec(:,:)               ! Panel tangential vectors in the length direction
        real, allocatable :: W2Vec(:,:)               ! Panel tangential vectors in the width direction
        real, allocatable :: W3Vec(:,:)               ! Panel normal vectors

        real, allocatable :: pnodes(:,:)              ! Panel node locations (x,y,z)
        
        ! wall strength
        real, allocatable :: WSource(:,:)             ! Wall source density values (column vector) (non-dimensional, normalized by freestream velocity)
        
        ! wall velocities (for output)
        real, allocatable :: u(:)                     ! u velocity over Uinf
        real, allocatable :: v(:)                     ! v velocity over Uinf
        real, allocatable :: w(:)                     ! w velocity over Uinf
        real, allocatable :: ur(:)                    ! velocity mag over Uinf

    end type WallType

    integer                     :: Nwalls             ! number of walls
    integer                     :: NumWP_total        ! total number of wall panels
    integer                     :: NumNodes_total     ! total number of nodes
    type(WallType), allocatable :: Walls(:)

contains

    subroutine wall_cns(WallInd,NumWP,NumNodes)

        ! Constructor for the arrays in this module

        integer :: WallInd,NumWP,NumNodes

        allocate(Walls(WallInd)%WCPoints(NumWP,3))
        allocate(Walls(WallInd)%W1Vec(NumWP,3))
        allocate(Walls(WallInd)%W2Vec(NumWP,3))
        allocate(Walls(WallInd)%W3Vec(NumWP,3))

        allocate(Walls(WallInd)%pnodes(NumNodes,3))

        allocate(Walls(WallInd)%WSource(NumWP,1))

        allocate(Walls(WallInd)%u(NumWP))
        allocate(Walls(WallInd)%v(NumWP))
        allocate(Walls(WallInd)%w(NumWP))
        allocate(Walls(WallInd)%ur(NumWP))

    end subroutine wall_cns


    subroutine nodepoints_from_ij(Wall,i,j,p1,p2,p3,p4)
        ! nodepoints_from_ij() : returns the coordinates of the four node points of the
        !   quadrilateral source panel whose center is indexed by i,j

        type(WallType), intent(in) :: Wall
        integer, intent(in)        :: i, j
        real, intent(out)          :: p1(3), p2(3), p3(3), p4(3)

        p1 = Wall%pnodes(i  ,j  )
        p2 = Wall%pnodes(i+1,j  )
        p3 = Wall%pnodes(i+1,j+1)
        p4 = Wall%pnodes(i  ,j+1)

    end subroutine nodepoints_from_ij


    subroutine read_p3d_walls(WallMeshPath)
        ! read_p3d_walls() : Read in a wall mesh from a Plot3d grid file.
        !   Grid data from each block is loaded into a Wall type.
        !
        ! WallMeshPath (input) : Path to file containing a multi-block structured mesh.
    
        character(len=*)     :: WallMeshPath
        integer              :: nblocks
        integer, allocatable :: ni(:), nj(:), nk(:)
        real, allocatable    :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)
        integer              :: m,i,j
        integer, parameter   :: k=1
        integer              :: nnodes,npanels
        real                 :: p1(3), p2(3), p3(3), p4(3)
        integer              :: idx

        ! read the plot3d mesh
        Call read_p3d_multiblock(WallMeshPath,nblocks,ni,nj,nk,x,y,z)

        ! set the module variable
        Nwalls = nblocks

        ! allocate for Walls
        allocate(Walls(nblocks))

        ! set initial value fo NumWP_total
        NumWP_total = 0
        NumNodes_total = 0

        ! for each block (or each "wall panel")
        do m = 1,nblocks
            if (nk(m) > 1) then
                write(*,*) "Error: Wall dimension in 3rd direction is > 1 -- wall is a volume, not a wall!"
            end if

            ! compute the number of panels
            nnodes  = ni(m)*nj(m)
            npanels = (ni(m)-1)*(nj(m)-1)

            ! add to the total number of wall panels
            NumWP_total = NumWP_total + npanels
            NumNodes_total = NumNodes_total + nnodes

            ! set some WallType variables
            Walls(m)%NumWP1   = ni(m)-1
            Walls(m)%NumWP2   = nj(m)-1
            Walls(m)%NumWP    = npanels
            Walls(m)%NumNodes = nnodes

            ! call the wall constructor subroutine (allocator for internal wall variables)
            Call wall_cns(m,npanels,nnodes)

            ! loop through grid coordinates
            idx = 1
            
            do i=1,ni(m)-1
                do j=1,nj(m)-1
                    ! get cell corners (ordered clockwise from lowest i,j)
                    p1 = [x(i  ,j  ,k,m), y(i  ,j  ,k,m), z(i  ,j  ,k,m)]
                    p2 = [x(i+1,j  ,k,m), y(i+1,j  ,k,m), z(i+1,j  ,k,m)]
                    p3 = [x(i+1,j+1,k,m), y(i+1,j+1,k,m), z(i+1,j+1,k,m)]
                    p4 = [x(i  ,j+1,k,m), y(i  ,j+1,k,m), z(i  ,j+1,k,m)]

                    Walls(m)%WCPoints(idx,1:3) = 0.25*(p1+p2+p3+p4)                                ! panel center
                    Walls(m)%W1Vec(idx,1:3)    = (p4-p1)/sqrt(sum((p4-p1)**2))                     ! panel 1-tangential vector
                    Walls(m)%W2Vec(idx,1:3)    = (p2-p1)/sqrt(sum((p2-p1)**2))                     ! panel 2-tangential vector

                    ! compute ! panel normal vector
                    Call cross(Walls(m)%W1Vec(idx,1),Walls(m)%W1Vec(idx,2),Walls(m)%W1Vec(idx,3),&
                               Walls(m)%W2Vec(idx,1),Walls(m)%W2Vec(idx,2),Walls(m)%W2Vec(idx,3),&
                               Walls(m)%W3Vec(idx,1),Walls(m)%W3Vec(idx,2),Walls(m)%W3Vec(idx,3))    

                    idx = idx + 1
                end do
            end do

            ! store node locations in an array (flattening)
            idx = 1
            do i=1,ni(m)
                do j=1,nj(m)
                    ! store node locations in array
                    Walls(m)%pnodes(idx,1) = x(i,j,k,m)
                    Walls(m)%pnodes(idx,2) = y(i,j,k,m)
                    Walls(m)%pnodes(idx,3) = z(i,j,k,m)
                    idx = idx + 1
                end do
            end do
        end do

        !deallocate x,y,k,ni,nj,nk?

    end subroutine read_p3d_walls


    ! subroutine compute_wall_velocity()
    ! ! compute_wall_velocity() : Compute the wall velocity at all locations along the wall panels
    ! end subroutine compute_wall_velocity

end module wallgeom
