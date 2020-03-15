module output

    ! Arrays to be output to csv data files

    ! Output flags
    integer :: BladeElemOutFlag                ! Set to 1 to output loads data at each timestep for each element, 0 to omit this output
    integer :: DynStallOutFlag                ! Set to 1 to output dynamic stall diagnostic output

    ! Scale parameters and flow properties
    character(10000) :: ParamsOutHead = 'R (ft),Frontal Area (ft^2),RPM,U (ft/s),rho (slugs/ft^3),tempr (degF),vis (slugs/(ft*s)),q (lb/ft^2),TSR (-),ReM (-),FnR (-)'
    real :: ParamsOutData(1,11)                             ! Data

    ! Revolution average data
    character(10000) :: RevOutHead = 'Rev,Power Coeff. (-),Tip Power Coeff. (-),Torque Coeff. (-),Fx Coeff. (-),Fy Coeff. (-),Fz Coeff. (-),Power (kW),Torque (ft-lbs)'
    real :: RevOutData(1,9)                ! Revolution average data for each revolution

    ! Timestep data
    character(10000) :: TSOutHead = 'Normalized Time (-),Theta (rad),Rev,Torque Coeff. (-),Power Coeff. (-),Fx Coeff. (-),Fy Coeff. (-),Fz Coeff. (-)'
    character(1000) :: TSOutHeadBlade = 'Blade Fx Coeff. (-),Blade Fy Coeff. (-),Blade Fz Coeff. (-),Blade Torque Coeff. (-)'
    character(1000) :: TSOutHeadStrut = 'Strut Fx Coeff. (-),Strut Fy Coeff. (-),Strut Fz Coeff. (-),Strut Torque Coeff. (-)'
    integer :: TSOutNCols                                ! Number of timestep outputs
    real, allocatable :: TSOutData(:,:)                 ! Timestep data for each timestep

    ! Element loads data
    character(10000) :: BladeElemOutHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,x/R (-),y/R (-),z/R (-),AOA25 (deg),AOA50 (deg),AOA75 (deg),AdotNorm (-),Re (-),Mach (-),Ur (-),IndU (-),IndV (-),IndW (-),GB (?),CL (-),CD (-),CM25 (-),&
&CLCirc (-),CN (-),CT (-),Fx (-),Fy (-),Fz (-),te (-)'
    integer :: BladeElemOutNCols = 29                           ! Number of element loads outputs
    integer :: BladeElemOutRow                                 ! Rows of element row output
    real, allocatable :: BladeElemOutData(:,:)                 ! Element loads data for each timestep

    ! Dynamic stall diagnostic output
    integer :: DynStallOutType
    character(10000) :: DynStallOutBVHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,alpha (deg),adotnorm (-),alrefL (deg),alrefD (deg),DynamicFlagL,DynamicFlagD'
    integer :: DynStallOutBVNCols = 11                           ! Number of BV diagnostic outputs
    integer :: DynStallOutBVRow                                 ! Rows of element row output
    real, allocatable :: DynStallOutBVData(:,:)                 ! Diagnostic data for each timestep
    character(10000) :: DynStallOutLBHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,L1,L2,L3,L4,L5,L6,L7,L8,L9,Check'
    integer :: DynStallOutLBNCols = 15                           ! Number of LB diagnostic outputs
    integer :: DynStallOutLBRow                                 ! Rows of element row output
    real, allocatable :: DynStallOutLBData(:,:)                 ! Diagnostic data for each timestep

contains

    subroutine output_cns(MaxSeg, MaxBlades, MaxStruts, DSFlag)

        ! Constructor for the arrays in this module

        integer :: MaxSeg, MaxBlades, MaxStruts, DSFlag

        integer :: i

        ! Calculate number of timestep outputs
        TSOutNCols=8+4*MaxBlades+4*MaxStruts
        allocate(TSOutData(1,TSOutNCols))
        do i=1,MaxBlades
            TSOutHead=trim(TSOutHead)//','//trim(TSOutHeadBlade)
        end do
        do i=1,MaxStruts
            TSOutHead=trim(TSOutHead)//','//trim(TSOutHeadStrut)
        end do

        if (BladeElemOutFlag == 1) then
            allocate(BladeElemOutData(MaxSeg,BladeElemOutNCols))
        end if

        if (DynStallOutFlag == 1 .AND. DSFlag == 1) then
            DynStallOutType = 1
            allocate(DynStallOutBVData(MaxSeg,DynStallOutBVNCols))
        else if (DynStallOutFlag == 1 .AND. DSFlag == 2) then
            DynStallOutType = 2
            allocate(DynStallOutLBData(MaxSeg,DynStallOutLBNCols))
        end if


    end subroutine output_cns

end module output
