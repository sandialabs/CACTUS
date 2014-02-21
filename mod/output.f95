MODULE output

    ! Arrays to be output to csv data files

    ! Output flags
    integer :: Output_ELFlag                ! Set to 1 to output loads data at each timestep for each element, 0 to omit this output
    integer :: Output_DSFlag                ! Set to 1 to output dynamic stall diagnostic output

    ! Scale parameters and flow properties
    character(10000) :: Output_SFHead = 'R (ft),Frontal Area (ft^2),RPM,U (ft/s),rho (slugs/ft^3),tempr (degF),vis (slugs/(ft*s)),q (lb/ft^2),TSR (-),ReM (-),FnR (-)'
    real :: Output_SFData(1,11)                             ! Data

    ! Revolution average data
    character(10000) :: Output_RevHead = 'Rev,Power Coeff. (-),Tip Power Coeff. (-),Torque Coeff. (-),Fx Coeff. (-),Fy Coeff. (-),Fz Coeff. (-),Power (kW),Torque (ft-lbs)'
    real :: Output_RevData(1,9)                ! Revolution average data for each revolution

    ! Timestep data
    character(10000) :: Output_TSHead = 'Normalized Time (-),Theta (rad),Rev,Torque Coeff. (-),Power Coeff. (-),Fx Coeff. (-),Fy Coeff. (-),Fz Coeff. (-)'
    character(1000) :: Output_TSHeadBlade = 'Blade Fx Coeff. (-),Blade Fy Coeff. (-),Blade Fz Coeff. (-),Blade Torque Coeff. (-)'
    character(1000) :: Output_TSHeadStrut = 'Strut Fx Coeff. (-),Strut Fy Coeff. (-),Strut Fz Coeff. (-),Strut Torque Coeff. (-)'
    integer :: Output_NTSOut                                ! Number of timestep outputs
    real, allocatable :: Output_TSData(:,:)                 ! Timestep data for each timestep

    ! Element loads data
    character(10000) :: Output_ELHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,AOA25 (deg),AOA50 (deg),AOA75 (deg),AdotNorm (-),Re (-),Mach (-),Ur (-),CL (-),CD (-),CM25 (-),CLCirc (-),CN (-),CT (-),Fx (-),Fy (-),Fz (-),te (-)'
    integer :: Output_NELOut = 22                           ! Number of element loads outputs
    integer :: Output_ELRow                                 ! Rows of element row output
    real, allocatable :: Output_ELData(:,:)                 ! Element loads data for each timestep

    ! Dynamic stall diagnostic output
    integer :: Output_DSType 
    character(10000) :: Output_BVHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,alpha (deg),adotnorm (-),alrefL (deg),alrefD (deg),DynamicFlagL,DynamicFlagD'
    integer :: Output_NBVOut = 11                           ! Number of BV diagnostic outputs
    integer :: Output_BVRow                                 ! Rows of element row output
    real, allocatable :: Output_BVData(:,:)                 ! Diagnostic data for each timestep   
    character(10000) :: Output_LBHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,L1,L2,L3,L4,L5,L6,L7,L8,L9,Check'
    integer :: Output_NLBOut = 15                           ! Number of LB diagnostic outputs
    integer :: Output_LBRow                                 ! Rows of element row output
    real, allocatable :: Output_LBData(:,:)                 ! Diagnostic data for each timestep

CONTAINS

    SUBROUTINE output_cns(MaxSeg, MaxBlades, MaxStruts, DSFlag)

        ! Constructor for the arrays in this module

        integer :: MaxSeg, MaxBlades, MaxStruts, DSFlag

        integer :: i

        ! Calculate number of timestep outputs
        Output_NTSOut=8+4*MaxBlades+4*MaxStruts
        allocate(Output_TSData(1,Output_NTSOut))
        do i=1,MaxBlades
            Output_TSHead=trim(Output_TSHead)//','//trim(Output_TSHeadBlade)  
        end do
        do i=1,MaxStruts
            Output_TSHead=trim(Output_TSHead)//','//trim(Output_TSHeadStrut)  
        end do

        if (Output_ELFlag == 1) then
            allocate(Output_ELData(MaxSeg,Output_NELOut))
        end if

        if (Output_DSFlag == 1 .AND. DSFlag == 1) then
            Output_DSType = 1
            allocate(Output_BVData(MaxSeg,Output_NBVOut))
        else if (Output_DSFlag == 1 .AND. DSFlag == 2) then
            Output_DSType = 2
            allocate(Output_LBData(MaxSeg,Output_NLBOut))
        end if


    End SUBROUTINE output_cns

End MODULE output
