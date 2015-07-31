module wallgeom

    ! wall geometry and outputs module

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

contains

    subroutine wall_cns(Wall,NumWP,NumNodes)

        ! Constructor for the arrays in this module
        type(WallType) :: Wall
        integer        :: WallInd,NumWP,NumNodes

        allocate(Wall%WCPoints(NumWP,3))
        allocate(Wall%W1Vec(NumWP,3))
        allocate(Wall%W2Vec(NumWP,3))
        allocate(Wall%W3Vec(NumWP,3))
        allocate(Wall%pnodes(NumNodes,3))

        allocate(Wall%WSource(NumWP,1))

        allocate(Wall%u(NumWP))
        allocate(Wall%v(NumWP))
        allocate(Wall%w(NumWP))
        allocate(Wall%ur(NumWP))

    end subroutine wall_cns


    subroutine ij_from_ip_local(Wall,ip_local,i,j)
        ! ij_from_ip_local() : returns the i,j index of the panel with "flat" index given by ip_local
        
        type(WallType), intent(in) :: Wall
        integer, intent(in)        :: ip_local
        integer, intent(out)       :: i,j
        
        i = mod(ip_local-1,Wall%NumWP1) + 1
        j = (ip_local-1)/Wall%NumWP1    + 1 ! forced to integer

    end subroutine ij_from_ip_local


    subroutine nodepoints_from_ij_local(Wall,i,j,p1,p2,p3,p4)
        ! nodepoints_from_ij_local() : returns the coordinates of the four node points of the
        !   quadrilateral source panel whose center is indexed by i,j

        type(WallType), intent(in) :: Wall
        integer, intent(in)        :: i, j
        real, intent(out)          :: p1(3), p2(3), p3(3), p4(3)
        integer                    :: NumWP1

        NumWP1 = Wall%NumWP1

        p1 = Wall%pnodes(i   + (j-1)*(NumWP1+1),1:3)
        p2 = Wall%pnodes(i+1 + (j-1)*(NumWP1+1),1:3)
        p3 = Wall%pnodes(i+1 + (j  )*(NumWP1+1),1:3)
        p4 = Wall%pnodes(i   + (j  )*(NumWP1+1),1:3)

    end subroutine nodepoints_from_ij_local


end module wallgeom
