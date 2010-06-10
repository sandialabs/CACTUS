MODULE freestream
	
	! Freestream velocity data

	real, allocatable :: UFS(:,:)		! Lattice point freestream x velocity 
	real, allocatable :: VFS(:,:)		! Lattice point freestream y velocity 
	real, allocatable :: WFS(:,:)		! Lattice point freestream z velocity 

	CONTAINS

	SUBROUTINE freestream_cns(MaxWakeNodes, MaxSegEnds)

		! Constructor for the arrays in this module

		integer :: MaxWakeNodes, MaxSegEnds
		
		allocate(UFS(MaxWakeNodes,MaxSegEnds))
		allocate(VFS(MaxWakeNodes,MaxSegEnds))
		allocate(WFS(MaxWakeNodes,MaxSegEnds))
		
	End SUBROUTINE

End