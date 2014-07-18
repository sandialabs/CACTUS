MODULE wakedata

	! Wake visualization data for WriteWakeData

    integer :: WakeElementOutFlag, WakePlaneOutFlag
    integer, allocatable :: WakeLineInd(:) 
    integer :: NWakeInd     
    character(1000) :: WakeOutHead = 'Timestep,Element,X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'

    ! Wake deficit calculation performed on a horizontal plane
    character(1000) :: HGridVelOutHead = 'Timestep,X/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'
    integer :: nxhgrid
    integer :: nzhgrid
    real :: yhgrid
    real :: xhgridL
    real :: xhgridU
    real :: zhgridL
    real :: zhgridU
    real, allocatable :: XHGrid(:,:) 
    real, allocatable :: ZHGrid(:,:) 
    real, allocatable :: VXIndH(:,:) 
    real, allocatable :: VYIndH(:,:) 
    real, allocatable :: VZIndH(:,:) 

    ! Wake deficit calculation performed on a vertical plane
    character(1000) :: VGridVelOutHead = 'Timestep,X/R (-),Y/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'
    integer :: nxvgrid
    integer :: nyvgrid
    real :: zvgrid
    real :: xvgridL
    real :: xvgridU
    real :: yvgridL
    real :: yvgridU
    real, allocatable :: XVGrid(:,:) 
    real, allocatable :: YVGrid(:,:) 
    real, allocatable :: VXIndV(:,:) 
    real, allocatable :: VYIndV(:,:) 
    real, allocatable :: VZIndV(:,:) 

    ! Wake deficit calculation performed on a cross-section plane
    character(1000) :: CGridVelOutHead = 'Timestep,Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'
    integer :: nycgrid
    integer :: nzcgrid
    real :: xcgrid
    real :: ycgridL
    real :: ycgridU
    real :: zcgridL
    real :: zcgridU
    real, allocatable :: YCGrid(:,:)
    real, allocatable :: ZCGrid(:,:)
    real, allocatable :: VXIndC(:,:)
    real, allocatable :: VYIndC(:,:)
    real, allocatable :: VZIndC(:,:)

    ! global counter
    integer :: ntcount


CONTAINS

	SUBROUTINE wakedata_cns()

     ! Constructor for the arrays in this module

        allocate(WakeLineInd(NWakeInd))     

        ! Wake deficit output, horizontal plane
        allocate(XHGrid(nxhgrid,nzhgrid))
        allocate(ZHGrid(nxhgrid,nzhgrid))
        allocate(VXIndH(nxhgrid,nzhgrid))    
        allocate(VYIndH(nxhgrid,nzhgrid)) 
        allocate(VZIndH(nxhgrid,nzhgrid))    

        ! Wake deficit output, vertical plane
        allocate(XVGrid(nxvgrid,nyvgrid))
        allocate(YVGrid(nxvgrid,nyvgrid))
        allocate(VXIndV(nxvgrid,nyvgrid))    
        allocate(VYIndV(nxvgrid,nyvgrid)) 
        allocate(VZIndV(nxvgrid,nyvgrid)) 

        ! Wake deficit output, cross-section plane
        allocate(YCGrid(nycgrid,nzcgrid))
        allocate(ZCGrid(nycgrid,nzcgrid))
        allocate(VXIndC(nycgrid,nzcgrid))
        allocate(VYIndC(nycgrid,nzcgrid))
        allocate(VZIndC(nycgrid,nzcgrid))

	End SUBROUTINE wakedata_cns

End MODULE wakedata
