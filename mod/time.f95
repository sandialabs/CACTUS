MODULE time

	! Info about program real time usage

	real :: t0				! Time in seconds marker for the beginning of the performance iteration
	real :: Time1			! Time at end of last rev
	real :: Time2			! Time at end of this rev
	real :: dtime           ! Delta time over this rev
	real :: etime           ! Elapsed time from begining at end of this rev

End MODULE time
