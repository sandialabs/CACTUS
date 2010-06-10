MODULE element

	! Blade element geometry data
	
	real, allocatable :: xBE(:,:)		! X location for each blade segment end (quarter chord) at each theta in a revolution	
	real, allocatable :: yBE(:,:)		! Y location for each blade segment end (quarter chord) at each theta in a revolution		
	real, allocatable :: zBE(:,:)		! Z location for each blade segment end (quarter chord) at each theta in a revolution	
	
	real :: hr				! Height to radius ratio. (VAWT only)
	real :: hubrr				! Hub radius ratio (Currently only used for HAWT calculation)
	real :: bi				! Blade incidence (Currently only used for HAWT calculation)
	real :: bCone				! Blade coning angle (Currently only use for HAWT calculation)
	real :: Tilt				! Rotor Tilt angle (Currently only use for HAWT calculation)
	real :: eta				! Blade mount point ratio (distance behind leading edge of the blade mount point / root chord)
      	real, allocatable :: btw(:)		! Blade twist for each element end in a blade (root to tip). (Currently only used for HAWT calculation)	
      	real, allocatable :: cr(:)		! Chord to radius ratio for each element end in a blade (root to tip)	
	real, allocatable :: nx(:,:)		! Normal X for each blade segment at each theta in a revolution	
	real, allocatable :: ny(:,:)		! Normal Y for each blade segment at each theta in a revolution	
	real, allocatable :: nz(:,:)		! Normal Z for each blade segment at each theta in a revolution	
	real, allocatable :: tx(:,:)		! Tangential X for each blade segment at each theta in a revolution	
	real, allocatable :: ty(:,:)		! Tangential Y for each blade segment at each theta in a revolution	
	real, allocatable :: tz(:,:)		! Tangential Z for each blade segment at each theta in a revolution	
	real, allocatable :: eSpan(:)		! Element span to equitorial radius ratio for each element
	real, allocatable :: eChord(:)		! Element chord to equitorial radius ratio for each element
	real, allocatable :: Theta(:)		! Theta at each point in a revolution (for the first blade)
	

	CONTAINS

	SUBROUTINE element_cns(MaxTimeStepPerRev,MaxSegEnds,MaxSegEndPerBlade)

		! Constructor for the arrays in this module

		integer :: MaxTimeStepPerRev,MaxSegEnds,MaxSegEndPerBlade
		
		allocate(xBE(MaxTimeStepPerRev,MaxSegEnds))
		allocate(yBE(MaxTimeStepPerRev,MaxSegEnds))
		allocate(zBE(MaxTimeStepPerRev,MaxSegEnds))
		allocate(btw(MaxSegEndPerBlade))
		allocate(cr(MaxSegEndPerBlade))
		allocate(nx(MaxTimeStepPerRev,MaxSegEnds))		
		allocate(ny(MaxTimeStepPerRev,MaxSegEnds))		
		allocate(nz(MaxTimeStepPerRev,MaxSegEnds))		
		allocate(tx(MaxTimeStepPerRev,MaxSegEnds))		
		allocate(ty(MaxTimeStepPerRev,MaxSegEnds))		
		allocate(tz(MaxTimeStepPerRev,MaxSegEnds))		
		allocate(eSpan(MaxSegEnds))	
		allocate(eChord(MaxSegEnds))			
		allocate(Theta(MaxTimeStepPerRev))
		
	End SUBROUTINE
	
End
