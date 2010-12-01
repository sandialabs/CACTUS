MODULE airfoil
	
	! Airfoil section designation data
	
	character*80, allocatable :: aftitle(:)  	! Title for each airfoil section
	character*4, allocatable ::  camber(:)		! Camber designation string ('IN','OUT', or 'NONE') for each section
	integer, allocatable :: camb(:)			! Camber flag (derived from camber string) for each section
	real, allocatable :: tc(:)			! Thickness to chord ratio for each section      
	

	CONTAINS

	SUBROUTINE airfoil_cns(MaxAirfoilSect)

		! Constructor for the arrays in this module

		integer :: MaxAirfoilSect
		
		allocate(aftitle(MaxAirfoilSect))
		allocate(camber(MaxAirfoilSect))
		allocate(camb(MaxAirfoilSect))
		allocate(tc(MaxAirfoilSect))		

		
	End SUBROUTINE

End