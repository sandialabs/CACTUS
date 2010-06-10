MODULE vel

	! Wake lattice point velocity data

	real, allocatable :: U(:,:)		! Lattice point x velocity for each wake point
	real, allocatable :: V(:,:)		! Lattice point y velocity for each wake point
	real, allocatable :: W(:,:)		! Lattice point z velocity for each wake point

	real, allocatable :: UIWake(:)		! Velocity induced at blade from wake
	real, allocatable :: VIWake(:)		! Velocity induced at blade from wake
	real, allocatable :: WIWake(:)		! Velocity induced at blade from wake

	CONTAINS

	SUBROUTINE vel_cns(MaxWakeNodes, MaxSegEnds)

		! Constructor for the arrays in this module

		integer :: MaxWakeNodes, MaxSegEnds
		
		allocate(U(MaxWakeNodes,MaxSegEnds))
		allocate(V(MaxWakeNodes,MaxSegEnds))
		allocate(W(MaxWakeNodes,MaxSegEnds))

		allocate(UIWake(MaxSegEnds))
		allocate(VIWake(MaxSegEnds))
		allocate(WIWake(MaxSegEnds))
		
	End SUBROUTINE

End
