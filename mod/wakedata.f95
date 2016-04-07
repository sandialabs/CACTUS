module wakedata

    implicit none

    ! Wake visualization data for WriteWakeData

    integer :: WakeElemOutFlag
    integer, allocatable :: WakeLineInd(:)
    integer :: NWakeInd
    character(1000) :: WakeOutHead = 'Normalized Time (-),Node ID,Origin Node,X/R (-),Y/R (-),Z/R (-),U/Uinf (-),V/Uinf (-),W/Uinf (-)'

contains

    subroutine wakedata_cns()

     ! Constructor for the arrays in this module

        allocate(WakeLineInd(NWakeInd))

    end subroutine wakedata_cns

end module wakedata
