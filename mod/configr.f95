MODULE configr
	
	! Configuration data
	     
        integer :: DiagOutFlag          ! Set to 1 to print iteration info to stdout      
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
	integer :: Istraight		! Set to 1 for straight-bladed VAWT, 0 for parabolic blade shape
	integer :: Istrut 		! Set to 1 for equatorial ! cross-flow turbine blade struts, 0 for no struts
        real :: sthick                  ! Strut sectional thickness-to-chord ratio
        real :: Cdpar                   ! Additional strut parasitic drag coefficient based on equatorial blade "chord area" (chord squared)
        real :: CTExcrM                 ! Additional machine level excrescence torque based on tip speed and Rmax
        real :: VCRFB                   ! Vortex core radius factor (on max blade chord) for bound vortex
        real :: VCRFT                   ! Vortex core radius factor (on blade discretization level) for trailing wake vorticies
        real :: VCRFS                   ! Vortex core radius factor (on temporal discretization level) for spanwise wake vorticies
        integer :: PRFlag               ! 0 for no pitch rate aero effects, 1 to include these effects
	
	real :: ut				! Tip speed ratio
        real :: dt                              ! Normalized timestep       
	real :: wRotX				! Normalized machine angular velocity X
	real :: wRotY				! Normalized machine angular velocity Y
	real :: wRotZ				! Normalized machine angular velocity Z
	real :: RotX				! Machine rotation axis X
	real :: RotY				! Machine rotation axis Y
	real :: RotZ				! Machine rotation axis Z
	real :: convrg				! Convergence level
	real :: convrgf				! Convergence level to be used for final convergence step (if ifc = 1)	
	real :: AT				! Normalized frontal area (frontal area / (equitorial radius)^2)
	real :: AreaT				! Frontal area
        real :: TorqueC                         ! Output torque normalization
        real :: Torque                          ! Output torque for the current revolution       
	real :: PowerC				! Output power normalization
	real :: Power				! Output power for the current revolution
	real :: ReM				! Machine Reynolds number based on Rmax
	
End
