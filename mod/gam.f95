MODULE gam
	
	! Circulation data
	
	real, allocatable :: GT(:,:)		! Trailing wake (streamwise) vorticity
	real, allocatable :: GS(:,:)		! Shed wake (spanwise) vorticity
	real, allocatable :: GB(:)		! Bound vorticity
	real, allocatable :: OGB(:)		! Old bound vorticity (previous time step)
		
	CONTAINS

	SUBROUTINE gam_cns(MaxWakeNodes, MaxSegEnds)

		! Constructor for the arrays in this module

		integer :: MaxWakeNodes, MaxSegEnds
	
		allocate(GT(MaxWakeNodes,MaxSegEnds))
		allocate(GS(MaxWakeNodes,MaxSegEnds))
		allocate(GB(MaxSegEnds))
		allocate(OGB(MaxSegEnds))		
		
	End SUBROUTINE

End