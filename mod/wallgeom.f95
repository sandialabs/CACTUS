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
        
        real, allocatable :: WCPoints(:,:)            ! Panel center points (over radius)
        real, allocatable :: W1Vec(:,:)               ! Panel tangential vectors in the length direction
        real, allocatable :: W2Vec(:,:)               ! Panel tangential vectors in the width direction
        real, allocatable :: W3Vec(:,:)               ! Panel normal vectors
        real, allocatable :: WP1(:)                   ! Panel lengths (over radius)
        real, allocatable :: WP2(:)                   ! Panel widths (over radius)
        
        ! wall strength
        real, allocatable :: WSource(:,:)             ! Wall source density values (column vector) (non-dimensional, normalized by freestream velocity)
        
        ! wall velocities (for output)
        real, allocatable :: u(:)                     ! u velocity over Uinf
        real, allocatable :: v(:)                     ! v velocity over Uinf
        real, allocatable :: w(:)                     ! w velocity over Uinf
        real, allocatable :: ur(:)                    ! velocity mag over Uinf

    end type WallType

    integer                     :: Nwalls             ! number of walls
    type(WallType), allocatable :: Walls(:)

contains

    subroutine wall_cns(WallInd,NumWP)

        ! Constructor for the arrays in this module

        integer :: WallInd,NumWP

        allocate(Walls(WallInd)%WCPoints(NumWP,3))
        allocate(Walls(WallInd)%W1Vec(NumWP,3))
        allocate(Walls(WallInd)%W2Vec(NumWP,3))
        allocate(Walls(WallInd)%W3Vec(NumWP,3))
        allocate(Walls(WallInd)%WP1(NumWP))
        allocate(Walls(WallInd)%WP2(NumWP))

        allocate(Walls(WallInd)%WSource(NumWP,1))

        allocate(Walls(WallInd)%u(NumWP))
        allocate(Walls(WallInd)%v(NumWP))
        allocate(Walls(WallInd)%w(NumWP))
        allocate(Walls(WallInd)%ur(NumWP))

    end subroutine wall_cns


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

        ! for each block (or each "wall panel")
        do m = 1,nblocks
            if (nk(m) > 1) then
                write(*,*) "Error: Wall dimension in 3rd direction is > 1 -- wall is a volume, not a wall!"
            end if

            ! compute the number of panels
            nnodes  = ni(m)*nj(m)
            npanels = (ni(m)-1)*(nj(m)-1)

            ! set some WallType variables
            Walls(m)%NumWP1 = ni(m)-1
            Walls(m)%NumWP2 = nj(m)-1
            Walls(m)%NumWP  = npanels

            ! call the wall constructor subroutine (allocator for internal wall variables)
            Call wall_cns(m,npanels)

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

                    Walls(m)%W1Vec(idx,1:3)    = p2-p1                                             ! panel x tangential vector
                    Walls(m)%WP1(idx)          = sqrt(sum((p2-p1)**2))                             ! panel x length
                    Walls(m)%W1Vec(idx,1:3)    = Walls(m)%W1Vec(idx,1:3)/Walls(m)%WP1(idx)         ! normalize

                    Walls(m)%W2Vec(idx,1:3)    = p4-p1                                             ! panel y tangential vector, set so that panel normal will be in the domain inward direction
                    Walls(m)%WP2(idx)          = sqrt(sum((p4-p1)**2))                             ! panel y length
                    Walls(m)%W2Vec(idx,1:3)    = Walls(m)%W2Vec(idx,1:3)/Walls(m)%WP2(idx)         ! normalize

                    ! compute ! panel normal vector
                    Call cross(Walls(m)%W1Vec(idx,1),Walls(m)%W1Vec(idx,2),Walls(m)%W1Vec(idx,3),&
                               Walls(m)%W2Vec(idx,1),Walls(m)%W2Vec(idx,2),Walls(m)%W2Vec(idx,3),&
                               Walls(m)%W3Vec(idx,1),Walls(m)%W3Vec(idx,2),Walls(m)%W3Vec(idx,3))    

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
