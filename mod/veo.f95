MODULE veo
	
	! Old wake lattice point velocities
	
	real, allocatable :: UO(:,:)		! Last lattice point x velocity for each wake point
	real, allocatable :: VO(:,:)		! Last lattice point y velocity for each wake point
	real, allocatable :: WO(:,:)		! Last lattice point z velocity for each wake point

	CONTAINS

	SUBROUTINE veo_cns(MaxWakeNodes, MaxSegEnds)

		! Constructor for the arrays in this module
		
		integer :: MaxWakeNodes, MaxSegEnds
		
		allocate(UO(MaxWakeNodes,MaxSegEnds))
		allocate(VO(MaxWakeNodes,MaxSegEnds))
		allocate(WO(MaxWakeNodes,MaxSegEnds))
		
	End SUBROUTINE
	
End