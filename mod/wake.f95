MODULE wake

	! Wake data

    ! Wake circulation data
    real, allocatable :: GT(:,:)            ! Trailing wake (streamwise) vorticity
    real, allocatable :: GS(:,:)            ! Shed wake (spanwise) vorticity

    ! Wake location data 
    real, allocatable :: X(:,:)             ! X position
    real, allocatable :: Y(:,:)             ! Y position
    real, allocatable :: Z(:,:)             ! Z position

    ! Induced velocity on wake lattice points
    real, allocatable :: U(:,:)             ! Lattice point x velocity for each wake point
    real, allocatable :: V(:,:)             ! Lattice point y velocity for each wake point
    real, allocatable :: W(:,:)             ! Lattice point z velocity for each wake point

    ! Freestream velocity at wake lattice points
    real, allocatable :: UFS(:,:)		! Lattice point freestream x velocity
    real, allocatable :: VFS(:,:)		! Lattice point freestream y velocity
    real, allocatable :: WFS(:,:)		! Lattice point freestream z velocity

    ! Old wake lattice point velocities
    real, allocatable :: UO(:,:)            ! Last lattice point x velocity for each wake point
    real, allocatable :: VO(:,:)            ! Last lattice point y velocity for each wake point
    real, allocatable :: WO(:,:)            ! Last lattice point z velocity for each wake point

CONTAINS

    SUBROUTINE wake_cns(MaxWakeNodes, MaxSegEnds)

		! Constructor - allocates memory for arrays
		
        integer :: MaxWakeNodes, MaxSegEnds

        allocate(GT(MaxWakeNodes,MaxSegEnds))
        allocate(GS(MaxWakeNodes+1,MaxSegEnds))         ! needs extra spanwise station for shedvor
        allocate(X(MaxWakeNodes,MaxSegEnds))
        allocate(Y(MaxWakeNodes,MaxSegEnds))
        allocate(Z(MaxWakeNodes,MaxSegEnds))
		allocate(U(MaxWakeNodes,MaxSegEnds))
        allocate(V(MaxWakeNodes,MaxSegEnds))
        allocate(W(MaxWakeNodes,MaxSegEnds))
        allocate(UFS(MaxWakeNodes,MaxSegEnds))
        allocate(VFS(MaxWakeNodes,MaxSegEnds))
        allocate(WFS(MaxWakeNodes,MaxSegEnds))
        allocate(UO(MaxWakeNodes,MaxSegEnds))
        allocate(VO(MaxWakeNodes,MaxSegEnds))
        allocate(WO(MaxWakeNodes,MaxSegEnds))    

    End SUBROUTINE wake_cns

End MODULE wake
