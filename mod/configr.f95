MODULE configr
	
	! Configuration data
	     
	integer :: GeomFlag		! Set to 1 for VAWT calculation, 0 for HAWT calculation
	integer :: nb			! Number of blades
	integer :: nbe			! Number of blade segments/elements in a blade
	integer :: ne			! Number of blade segment ends (total)
	integer :: iut			! Number of time steps over which the wake velocities are left constant
	integer :: nr			! Number of revolutions to perform
	integer :: iRev			! Revolution counter
	integer :: nt			! Time step counter
	integer :: nti			! Number of time steps per revolution
	integer :: ntTerm		! Time step level beyond which wake is ignored (if iXTerm = 1)
	integer :: nSect		! Number of airfoil section data tables used
	integer :: FitStartRev		! Initialization rev for the cpf fit
	integer :: npw			! Number of points in the wake
	integer :: nsw			! Next iteration at which wake velocities will be calculated
	
	real :: CrRef				! Reference chord to radius ratio
	real :: ut				! Tip speed ratio
	real :: wRotX				! Normalized machine angular velocity X
	real :: wRotY				! Normalized machine angular velocity Y
	real :: wRotZ				! Normalized machine angular velocity Z
	real :: RotX				! Machine rotation axis X
	real :: RotY				! Machine rotation axis Y
	real :: RotZ				! Machine rotation axis Z
	real :: convrg				! Convergence level
	real :: convrgf				! Convergence level to be used for final convergence step (if ifc = 1)	
	real, allocatable :: Cp(:)		! Average power coeff for each revolution
	real, allocatable :: Kp(:)		! Average power coeff (based on tip speed rather than Uinf?) for each revolution
	real, allocatable :: Cpf(:)		! Estimate of infinite rev Cp for each revolution based on moving fit of last three Cp values 
	real, allocatable :: Kpf(:)		! Estimate of infinite rev Kp for each revolution based on moving fit of last three Cp values 
	real :: AT				! Normalized frontal area (frontal area / (equitorial radius)^2)
	real :: AreaT				! Frontal area
	real :: PowerC				! Output power normalization
	real :: uMPH				! Uinf in mph
	real :: Power				! Output power for the current revolution
	real :: ReM				! Reynolds number based on blade chord

	CONTAINS

	SUBROUTINE configr_cns(MaxRevs)

		! Constructor for the arrays in this module

		integer :: MaxRevs
		
		allocate(Cp(MaxRevs))
		allocate(Kp(MaxRevs))
		allocate(Cpf(MaxRevs))
		allocate(Kpf(MaxRevs))		
		
	End SUBROUTINE
	
End
