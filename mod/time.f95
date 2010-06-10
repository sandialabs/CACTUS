MODULE time

	! Info about program real time usage
	
	character :: DMY*9, HMS*8		! Time strings for program run date and time
	real :: t0				! Time in seconds marker for the beginning of the performance iteration
	real :: Time1				! Delta time from beginning of the last iteration
	real :: TimeF				! Delta time at the end of the performance loop
	real, allocatable :: DTime(:)		! Seconds used for this revolution
	real, allocatable :: ETime(:)		! Elapsed time in seconds since beginning of code after each revolution

	CONTAINS

	SUBROUTINE time_cns(MaxRevs)

		! Constructor for the arrays in this module

		integer :: MaxRevs
		
		allocate(DTime(MaxRevs))
		allocate(ETime(MaxRevs))
		
	End SUBROUTINE
	
End