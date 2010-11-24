MODULE cltab
	
	! Airfoil section coefficient data
	
	character*80, allocatable :: dftitle(:)	! Airfoil section data title
	real, allocatable :: TA(:,:,:)		! Table AOA values
	real, allocatable :: TCL(:,:,:)		! Table CL values
	real, allocatable :: TCD(:,:,:)		! Table CD values
	real, allocatable :: TRE(:,:)		! Table Re values
	integer, allocatable :: nTBL(:,:)	! Number of AOA values for each Re number, in each section data table
        integer, allocatable :: nRET(:)         ! Number of Re number values in each section data table       
	integer, allocatable :: iSect(:)	! Array of indicies of the section table to apply to each blade element
	

	CONTAINS

	SUBROUTINE cltab_cns(MaxAOAVals,MaxReVals,MaxAirfoilSect,MaxSegEnds)

		! Constructor for the arrays in this module

		integer :: MaxAOAVals,MaxReVals,MaxAirfoilSect,MaxSegEnds
		
		allocate(dftitle(MaxAirfoilSect))
		allocate(TA(MaxAOAVals,MaxReVals,MaxAirfoilSect))
		allocate(TCL(MaxAOAVals,MaxReVals,MaxAirfoilSect))
		allocate(TCD(MaxAOAVals,MaxReVals,MaxAirfoilSect))
		allocate(TRE(MaxReVals,MaxAirfoilSect))		
		allocate(nTBL(MaxReVals,MaxAirfoilSect))
                allocate(nRET(MaxAirfoilSect))              		
		allocate(iSect(MaxSegEnds))		
		
	End SUBROUTINE

End