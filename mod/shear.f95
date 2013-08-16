MODULE shear

	! Ground shear layer model data

	real :: yGC		! Ground clearance height to radius ratio
	real :: yRef		! Freestream reference height to radius ratio (99% location maybe)
	real :: slex		! Exponent for shear layer model (Ex. 1/2 for parabolic BL model, 0 for constant freestream...)	
	real :: Tempr		! Temperature in degF

End MODULE shear
