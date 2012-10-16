MODULE configr
	
        ! Configuration data
	     
        integer :: DiagOutFlag          ! Set to 1 to print iteration info to stdout      
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
        integer :: nsw			! Next iteration at which wake velocities will be calculated
        integer :: nsWall               ! Next iteration at which the wall models will be updated     
        integer :: Istraight		! Set to 1 for straight-bladed VAWT, 0 for parabolic blade shape
        real :: convrg                  ! Convergence level
        real :: convrgf                 ! Convergence level to be used for final convergence step (if ifc = 1)         
        real :: CTExcrM                 ! Additional machine level excrescence torque based on tip speed and Rmax
        real :: VCRFB                   ! Vortex core radius factor (on max blade chord) for bound vortex
        real :: VCRFT                   ! Vortex core radius factor (on blade discretization level) for trailing wake vorticies
        real :: VCRFS                   ! Vortex core radius factor (on temporal discretization level) for spanwise wake vorticies
        integer :: PRFlag               ! 0 for no pitch rate aero effects, 1 to include these effects
	
        real :: ut				! Tip speed ratio       
        real :: dt                              ! Normalized timestep
        real :: delt                            ! Phase angle step       
        real :: wRotX				! Normalized machine angular velocity X
        real :: wRotY				! Normalized machine angular velocity Y
        real :: wRotZ				! Normalized machine angular velocity Z
        real :: RotX				! Machine rotation axis X
        real :: RotY				! Machine rotation axis Y
        real :: RotZ				! Machine rotation axis Z
        real :: RotPX                           ! Rotation origin X
        real :: RotPY                           ! Rotation origin Y
        real :: RotPZ                           ! Rotation origin Z  
        real :: ReM                             ! Machine Reynolds number based on Rmax
        real :: MInf                            ! Freestream mach number       

        real :: ForceC                          ! Output force normalization
        real :: TorqueC                         ! Output torque normalization
        real :: PowerC				! Output power normalization
        real :: AT                              ! Normalized frontal area (frontal area / (equitorial radius)^2)
        real :: AreaT                           ! Frontal area
        
        ! Current normalized time and phase angle (Theta)
        real :: Theta                           ! Current turbine rotation phase angle
        real :: TimeN                           ! Current normalized time
        
        ! Rev average quantities
        real :: CPAve                           ! Rev average power coeff
        real :: KPAve                           ! Rev average tip power coeff
        real :: CTRAve                          ! Rev average torque coeff
        real :: CFxAve                          ! Rev average Fx coeff
        real :: CFyAve                          ! Rev average Fy coeff
        real :: CFzAve                          ! Rev average Fz coeff
        real :: PowerAve                        ! Rev average power
        real :: TorqueAve                       ! Rev average torque
        real :: CPSum
        real :: CTRSum
        real :: CFxSum
        real :: CFySum
        real :: CFzSum
	
End
