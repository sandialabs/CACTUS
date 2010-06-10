MODULE xwake

	! Fixed wake grid data
	
	real, allocatable :: XFW(:)		! Fixed wake grid x stations
	real, allocatable :: YFW(:)		! Fixed wake grid x stations
	real, allocatable :: ZFW(:)		! Fixed wake grid x stations
	
	integer iFW			! Number of fixed wake grid points in X
	integer jFW			! Number of fixed wake grid points in Y
	integer kFW			! Number of fixed wake grid points in Z	
	

	CONTAINS

	SUBROUTINE xwake_cns(MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ)

		! Constructor for the arrays in this module

		integer :: MaxFixWakeX,MaxFixWakeY,MaxFixWakeZ
		
		allocate(XFW(MaxFixWakeX))
		allocate(YFW(MaxFixWakeY))
		allocate(ZFW(MaxFixWakeZ))
		
	End SUBROUTINE
	
End