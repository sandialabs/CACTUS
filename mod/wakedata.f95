MODULE wakedata

	! Wake visualization data for WriteWakeData
	
        integer :: WakeOutFlag
        integer, allocatable :: WakeLineInd(:) 
        integer :: NWakeInd     
        character(1000) :: WakeOutHead = 'Timestep,Element,X/R,Y/R,Z/R,U/Uinf,V/Uinf,W/Uinf'
        
        ! Wake deficit calculation performed if WakeOutFlag=2
        ! JCM test: wake deficit output plane is currently hardcoded...
        integer :: nxgrid = 141
        integer :: nzgrid = 81
        real :: ygrid = 1.2325
        real :: xgridL = -2
        real :: xgridU = 5
        real :: zgridL = -2
        real :: zgridU = 2
        integer :: xcount, zcount, ntcount
        real :: dxgrid, dzgrid
        real :: PointGrid(3), IndVelGrid(3)
        real, allocatable :: XGrid(:,:) 
        real, allocatable :: ZGrid(:,:) 
        real, allocatable :: SVDef(:,:) 
                 
	
	CONTAINS

	SUBROUTINE wakedata_cns()

		! Constructor for the arrays in this module

                allocate(WakeLineInd(NWakeInd))     
                         
                ! JCM test: wake deficit output
                allocate(XGrid(nxgrid,nzgrid))
                allocate(ZGrid(nxgrid,nzgrid))
                allocate(SVDef(nxgrid,nzgrid))         
		
	End SUBROUTINE
	
End