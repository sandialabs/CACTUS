MODULE wakedata

	! Wake visualization data for WriteWakeData
	
	real, allocatable :: WakeDefAve(:,:)
	integer :: WakeOut
	
	CONTAINS

	SUBROUTINE wakedata_cns(MaxFixWakeY,MaxFixWakeZ)

		! Constructor for the arrays in this module

		integer :: MaxFixWakeY,MaxFixWakeZ
		
		allocate(WakeDefAve(MaxFixWakeY,MaxFixWakeZ))
		
	End SUBROUTINE
	
End