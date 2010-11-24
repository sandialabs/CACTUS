MODULE wakedata

	! Wake visualization data for WriteWakeData
	
        integer :: WakeOutFlag
	real, allocatable :: WakeDefAve(:,:)
        integer, allocatable :: WakeLineInd(:) 
        integer :: NWakeInd     
        character(1000) :: WakeOutHead = 'Timestep,Element,X/R,Y/R,Z/R,U/Uinf,V/Uinf,W/Uinf'
                 
	
	CONTAINS

	SUBROUTINE wakedata_cns()

		! Constructor for the arrays in this module

                allocate(WakeLineInd(NWakeInd))              
		
	End SUBROUTINE
	
End