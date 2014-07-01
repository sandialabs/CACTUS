MODULE wakedata

	! Wake visualization data for WriteWakeData

    integer :: WakeOutFlag
    integer, allocatable :: WakeLineInd(:) 
    integer :: NWakeInd     
    character(1000) :: WakeOutHead = 'Timestep,Element,X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'

    ! Wake deficit calculation performed on a horizontal plane if WakeOutFlag=2
    ! JCM test: wake deficit output plane is currently hardcoded...
    character(1000) :: HGridVelOutHead = 'X/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'
    integer :: nxhgrid = 141
    integer :: nzhgrid = 81
    real :: yhgrid = 1.2325
    real :: xhgridL = -2
    real :: xhgridU = 5
    real :: zhgridL = -2
    real :: zhgridU = 2
    real, allocatable :: XHGrid(:,:) 
    real, allocatable :: ZHGrid(:,:) 
    real, allocatable :: VXIndH(:,:) 
    real, allocatable :: VYIndH(:,:) 
    real, allocatable :: VZIndH(:,:) 

    ! Wake deficit calculation performed on a vertical plane if WakeOutFlag=3
    ! JCM test: wake deficit output plane is currently hardcoded...
    character(1000) :: VGridVelOutHead = 'X/R (-),Y/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'
    integer :: nxvgrid = 141
    integer :: nyvgrid = 101
    real :: zvgrid = 0.0
    real :: xvgridL = -2
    real :: xvgridU = 5
    real :: yvgridL = -2
    real :: yvgridU = 3
    real, allocatable :: XVGrid(:,:) 
    real, allocatable :: YVGrid(:,:) 
    real, allocatable :: VXIndV(:,:) 
    real, allocatable :: VYIndV(:,:) 
    real, allocatable :: VZIndV(:,:) 

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

	End SUBROUTINE wakedata_cns

End MODULE wakedata
