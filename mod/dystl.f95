MODULE dystl
	
	! Dynamic stall table data

	real, allocatable :: GammaXL(:)			! Variables apparently related to the transonic effect of 
	real, allocatable :: GammaXM(:)			! section thickness/chord on the dynamic stall behaviour with AOA...
	real, allocatable :: dGammaL(:)			! Need description of Boeing-Vertol dynamic stall model...
	real, allocatable :: dGammaM(:)			! (do these tip speeds get that fast?)...
	real, allocatable :: SMachL(:)			! ...
	real, allocatable :: SMachM(:)			! ...
	real :: K1Pos					! ...
	real :: K1Neg					! ...
	
	real, allocatable :: alzer(:)			! Zero lift AOA for each section
	real, allocatable :: restl(:,:)			! Re numbers in the stall AOA table for each section
	real, allocatable :: alstlp(:,:)		! Stall AOA (positive) at all Re numbers for each section
	real, allocatable :: alstln(:,:)		! Stall AOA (negative) at all Re numbers for each section
	real, allocatable :: dapdre(:,:)		! (Delta positive stall angle) / (Delta Re number)
	real, allocatable :: dandre(:,:)		! (Delta negative stall angle) / (Delta Re number)
	
	real, allocatable :: alfold(:)			! last time value of alpha for each segment
	real :: MInf 					! Freestream mach number
	
	integer, allocatable :: nstl(:)			! Number of Re num values in the stall AOA table for each section
	
	CONTAINS
	
	SUBROUTINE dystl_cns(MaxAirfoilSect, MaxReVals, MaxSegEnds)
		
		! Constructor for the arrays in this module
		
		integer :: MaxAirfoilSect, MaxReVals, MaxSegEnds
		
		allocate(GammaXL(MaxAirfoilSect))
		allocate(GammaXM(MaxAirfoilSect))
		allocate(dGammaL(MaxAirfoilSect))
		allocate(dGammaM(MaxAirfoilSect))
		allocate(SMachL(MaxAirfoilSect))
		allocate(SMachM(MaxAirfoilSect))
		
		allocate(alzer(MaxAirfoilSect))
		allocate(restl(MaxReVals,MaxAirfoilSect))
		allocate(alstlp(MaxReVals,MaxAirfoilSect))
		allocate(alstln(MaxReVals,MaxAirfoilSect))
		allocate(dapdre(MaxReVals,MaxAirfoilSect))
		allocate(dandre(MaxReVals,MaxAirfoilSect))
		
		allocate(alfold(MaxSegEnds))
		allocate(nstl(MaxAirfoilSect))
	
	End SUBROUTINE
	
End 