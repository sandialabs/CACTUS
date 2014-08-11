MODULE wakedata

    ! Wake visualization data for WriteWakeData

    integer :: WakeElementOutFlag, WakeGridOutFlag
    integer, allocatable :: WakeLineInd(:) 
    integer :: NWakeInd     
    character(1000) :: WakeOutHead = 'Normalized Time (-),Element,X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'

    ! Wake deficit calculation on a grid
    character(1000) :: GridVelOutHead = 'Normalized Time (-),X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'

    ! number of grid elements in each direction    
    integer :: nxgrid
    integer :: nygrid
    integer :: nzgrid

    ! grid extents
    real :: xgridL
    real :: xgridU
    real :: ygridL
    real :: ygridU
    real :: zgridL
    real :: zgridU

    ! arrays holding grid locations
    real, allocatable :: XGrid(:,:,:) 
    real, allocatable :: YGrid(:,:,:) 
    real, allocatable :: ZGrid(:,:,:) 
    real, allocatable :: VXInd(:,:,:) 
    real, allocatable :: VYInd(:,:,:) 
    real, allocatable :: VZInd(:,:,:) 

    ! global counter
    integer :: ntcount


CONTAINS

    SUBROUTINE wakedata_cns()

     ! Constructor for the arrays in this module

        allocate(WakeLineInd(NWakeInd))     

        ! Wake deficit output, horizontal plane
        allocate(XGrid(nxgrid,nygrid,nzgrid))
        allocate(YGrid(nxgrid,nygrid,nzgrid))
        allocate(ZGrid(nxgrid,nygrid,nzgrid))
        allocate(VXInd(nxgrid,nygrid,nzgrid))    
        allocate(VYInd(nxgrid,nygrid,nzgrid)) 
        allocate(VZInd(nxgrid,nygrid,nzgrid))    

    End SUBROUTINE wakedata_cns

End MODULE wakedata
