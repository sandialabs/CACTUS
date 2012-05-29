MODULE ioption

	! Input options
	
	integer :: ivtxcor	! switches whether to use finite vortex core model...
	integer :: ifc		! Flag to use final convergence step
	integer :: nric		! Intermediate rev at which to switch to final convergence (optional)
	integer :: ntif		! final number of time steps per revolution (used instead of nti during final convergence step)
	integer :: iutf		! final number of time steps between updating the wake velocities (used instead of iut during final convergence step)
	integer :: iXTerm	! Flag to ignore wake points beyond xstop
	integer :: XStop	! see above
        integer :: BladeFileFlag ! Flag to read blade geometry from file 'blade.inp'
	
	character*80 :: jbtitle  ! job title
	character*80 :: BladeFileName = 'blade.inp'  ! Blade geometry file name
End
