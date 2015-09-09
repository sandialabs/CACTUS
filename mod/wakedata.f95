module wakedata

    implicit none

    ! Wake visualization data for WriteWakeData

    integer :: WakeElementOutFlag, WakeGridOutFlag
    integer, allocatable :: WakeLineInd(:)
    integer :: NWakeInd
    character(1000) :: WakeOutHead = 'Normalized Time (-),Node ID,Origin Node,X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'

    ! Wake deficit calculation on a grid
    character(1000) :: GridVelOutHead = 'Normalized Time (-),X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-),Ufs/Uinf (-),Vfs/Uinf (-),Wfs/Uinf (-)'

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

    ! grid spacing
    real :: dxgrid, dygrid, dzgrid

    ! arrays holding grid locations
    real, allocatable :: XGrid(:,:,:)
    real, allocatable :: YGrid(:,:,:)
    real, allocatable :: ZGrid(:,:,:)
    real, allocatable :: VXInd(:,:,:)
    real, allocatable :: VYInd(:,:,:)
    real, allocatable :: VZInd(:,:,:)
    real, allocatable :: UfsGrid(:,:,:)
    real, allocatable :: VfsGrid(:,:,:)
    real, allocatable :: WfsGrid(:,:,:)

    ! array index counters
    integer :: xcount, ycount, zcount

    ! global counter
    integer :: ntcount


contains

    subroutine wakedata_cns()

     ! Constructor for the arrays in this module

        allocate(WakeLineInd(NWakeInd))

        ! Wake deficit output, horizontal plane
        allocate(XGrid(nxgrid,nygrid,nzgrid))
        allocate(YGrid(nxgrid,nygrid,nzgrid))
        allocate(ZGrid(nxgrid,nygrid,nzgrid))
        allocate(VXInd(nxgrid,nygrid,nzgrid))
        allocate(VYInd(nxgrid,nygrid,nzgrid))
        allocate(VZInd(nxgrid,nygrid,nzgrid))
        allocate(UfsGrid(nxgrid,nygrid,nzgrid))
        allocate(VfsGrid(nxgrid,nygrid,nzgrid))
        allocate(WfsGrid(nxgrid,nygrid,nzgrid))

        !! Set up grid for induced velocity output

        ! Compute the appropriate grid spacing
        ! (if number of grids is 1 in any direction, catch a divide-by-zero error)
        if (nxgrid==1) then
            dxgrid=0.0
        else
            dxgrid=(xgridU-xgridL)/(nxgrid-1)
        end if

        if (nygrid==1) then
            dygrid=0.0
        else
            dygrid=(ygridU-ygridL)/(nygrid-1)
        end if

        if (nzgrid==1) then
            dzgrid=0.0
        else
            dzgrid=(zgridU-zgridL)/(nzgrid-1)
        end if

        ! Set up grid node locations and initialize induced velocities
        do zcount=1,nzgrid
            do ycount=1,nygrid
                do xcount=1,nxgrid
                    ! Setup grid node locations
                    XGrid(xcount,ycount,zcount)=xgridL+(xcount-1)*dxgrid
                    YGrid(xcount,ycount,zcount)=ygridL+(ycount-1)*dygrid
                    ZGrid(xcount,ycount,zcount)=zgridL+(zcount-1)*dzgrid

                    ! Initialize induced velocities to zero
                    VXInd(xcount,ycount,zcount)=0.0
                    VYInd(xcount,ycount,zcount)=0.0
                    VZInd(xcount,ycount,zcount)=0.0

                    ! Initialize freestream velocities to zero
                    UfsGrid(xcount,ycount,zcount)=0.0
                    VfsGrid(xcount,ycount,zcount)=0.0
                    WfsGrid(xcount,ycount,zcount)=0.0
                end do
            end do
        end do

    end subroutine wakedata_cns

end module wakedata
