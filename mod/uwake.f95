MODULE uwake

	! Fixed wake velocity data
	
	real, allocatable :: UFW(:,:,:)	! Fixed wake x velocities
	real, allocatable :: VFW(:,:,:)	! Fixed wake y velocities
	real, allocatable :: WFW(:,:,:)	! Fixed wake z velocities
		

	CONTAINS

	SUBROUTINE uwake_cns(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ)

		! Constructor for the arrays in this module

		integer :: MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ
		
		allocate(UFW(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ))
		allocate(VFW(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ))
		allocate(WFW(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ))
		
	End SUBROUTINE
	
End