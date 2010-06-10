MODULE wakeloc
	
	! Blade element and wake location data 

	real, allocatable :: X(:,:)		! X position history for each blade element
	real, allocatable :: Y(:,:)		! Y position history for each blade element
	real, allocatable :: Z(:,:)		! Z position history for each blade element
	
	CONTAINS
	
	SUBROUTINE wakeloc_cns(MaxWakeNodes, MaxSegEnds)
		
		! Constructor for the arrays in this module
		
		integer :: MaxWakeNodes, MaxSegEnds
		
		allocate(X(MaxWakeNodes,MaxSegEnds))
		allocate(Y(MaxWakeNodes,MaxSegEnds))
		allocate(Z(MaxWakeNodes,MaxSegEnds))
	
	End SUBROUTINE
	
End
