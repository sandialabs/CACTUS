MODULE configr

    ! Configuration data

    ! Input flags
    character*80 :: jbtitle         ! job title
    integer :: ifc                  ! Flag to use final convergence step
    integer :: iXTerm               ! Flag to ignore wake points beyond xstop 
    integer :: DiagOutFlag          ! Set to 1 to print iteration info to stdout   
    integer :: PRFlag               ! 0 for no pitch rate aero effects, 1 to include these effects

    ! Params  
    integer :: nb           ! Number of blades
    integer :: nbe          ! Number of blade segments/elements in a blade
    integer :: ne           ! Number of blade segment ends (total)
    integer :: iut          ! Number of time steps over which the wake velocities are left constant (Note: values of 0, -1, and -2 have defined meanings)
    integer :: nr           ! Number of revolutions to perform
    integer :: iRev         ! Revolution counter
    integer :: nt           ! Time step counter
    integer :: nti          ! Number of time steps per revolution
    integer :: nric                 ! Intermediate rev at which to switch to final convergence (if ifc = 1)
    integer :: ntif                 ! final number of time steps per revolution (used instead of nti during final convergence step if ifc = 1)
    integer :: iutf                 ! final number of time steps between updating the wake velocities (used instead of iut during final convergence step if ifc = 1)
    integer :: ntTerm       ! Time step level beyond which wake is ignored (if iXTerm = 1)
    integer :: nSect        ! Number of airfoil section data tables used
    integer :: nsw          ! Next iteration at which wake velocities will be calculated
    integer :: nsWall               ! Next iteration at which the wall models will be updated
    real :: convrg                  ! Convergence level (Note: for no convergence check, input -1)
    real :: convrgf                 ! Convergence level to be used for final convergence step (if ifc = 1) (Note: for no convergence check, input -1)
    real :: XStop                   ! If iXTerm = 1, ignore wake beyond x = xstop      
    real :: CTExcrM                 ! Additional machine level excrescence torque based on tip speed and Rmax
    real :: VCRFB                   ! Vortex core radius factor (on max blade chord) for bound vortex
    real :: VCRFT                   ! Vortex core radius factor (on blade discretization level) for trailing wake vorticies
    real :: VCRFS                   ! Vortex core radius factor (on temporal discretization level) for spanwise wake vorticies
    
    integer :: WakeElementOutIntervalTimesteps             ! Number of revolutions between writing wake data
    integer :: WakeElementOutStartTimestep                 ! Revolution number at which to start writing wake data
    integer :: WakeElementOutEndTimestep                   ! Revolution number at which to stop writing wake data

    integer :: WakeGridOutIntervalTimesteps                ! Number of revolutions between writing wake data
    integer :: WakeGridOutStartTimestep                    ! Revolution number at which to start writing wake data
    integer :: WakeGridOutEndTimestep                      ! Revolution number at which to stop writing wake data
    
    integer :: WallOutIntervalTimesteps                    ! Number of revolutions between writing wall data
    integer :: WallOutStartTimestep                        ! Revolution number at which to start writing wall data
    integer :: WallOutEndTimestep                          ! Revolution number at which to stop writing wall data

    real :: ut              ! Tip speed ratio       
    real :: dt                              ! Normalized timestep
    real :: delt                            ! Phase angle step       
    real :: wRotX               ! Normalized machine angular velocity X
    real :: wRotY               ! Normalized machine angular velocity Y
    real :: wRotZ               ! Normalized machine angular velocity Z
    real :: RotX                ! Machine rotation axis X
    real :: RotY                ! Machine rotation axis Y
    real :: RotZ                ! Machine rotation axis Z
    real :: RotPX                           ! Rotation origin X
    real :: RotPY                           ! Rotation origin Y
    real :: RotPZ                           ! Rotation origin Z  
    real :: ReM                             ! Machine Reynolds number based on Rmax
    real :: MInf                            ! Freestream mach number       

    real :: ForceC                          ! Output force normalization
    real :: TorqueC                         ! Output torque normalization
    real :: PowerC              ! Output power normalization
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

End MODULE configr
