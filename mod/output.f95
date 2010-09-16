MODULE output

        ! Arrays to be output to csv data files

        ! Output flags
        integer :: Output_ELFlag                ! Set to 1 to output loads data at each timestep for each element, 0 to omit this output
        
        ! Scale parameters and flow properties
        character(1000) :: Output_SFHead = 'R (ft),Frontal Area (ft^2),RPM,U (ft/s),rho (slugs/ft^3),tempr (degF),vis (slugs/(ft*s)),q (lb/ft^2),TSR (-),ReM (-)'
        real :: Output_SFData(1,10)                             ! Data
        
        ! Revolution average data
        character(1000) :: Output_RevHead = 'Rev,Power Coeff. (-),Tip Power Coeff. (-),Power (kW),Torque (ft-lbs),Delta CPU Time (s),Total CPU Time (s)'
        integer :: Output_NRevOut = 7                           ! Number of revolution average outputs
        integer :: Output_RevRow                                ! Rows of revolution average output
        real, allocatable :: Output_RevData(:,:)                ! Revolution average data for each revolution 
  
        ! Timestep data
        character(1000) :: Output_TSHead = 'Normalized Time (-),Rev,Torque Coeff. (-),Power Coeff. (-)'
        character(1000) :: Output_TSHeadBlade = 'Blade Fx Coeff. (-),Blade Fy Coeff. (-),Blade Fz Coeff. (-),Blade Torque Coeff. (-)'
        integer :: Output_NTSOut                                ! Number of timestep outputs
        integer :: Output_TSRow                                 ! Rows of timestep output
        real, allocatable :: Output_TSData(:,:)                 ! Timestep data for each timestep
        
        ! Element loads data
        character(1000) :: Output_ELHead = 'Normalized Time (-),Blade,Element,Rev,DynamicFlagL,DynamicFlagD,AOA (deg),Re (-),Mach (-),Ur (-),CN (-),CT (-),Fx (-),Fy (-),Fz (-),te (-)'
        integer :: Output_NELOut =16                            ! Number of element loads outputs
        integer :: Output_ELRow                                 ! Rows of element row output
        real, allocatable :: Output_ELData(:,:)                 ! Element loads data for each timestep
        
        CONTAINS
        
        SUBROUTINE output_cns(MaxRevs, MaxTimeSteps, MaxSeg, MaxBlades)
                
                ! Constructor for the arrays in this module
                
                integer :: MaxRevs, MaxTimeSteps, MaxSeg, MaxBlades
                
                integer :: i
                
                ! Calculate number of timestep outputs
                Output_NTSOut=4+4*MaxBlades
                do i=1,MaxBlades
                       Output_TSHead=trim(Output_TSHead)//','//trim(Output_TSHeadBlade)  
                end do
                
                allocate(Output_RevData(MaxRevs,Output_NRevOut))
                allocate(Output_TSData(MaxTimeSteps,Output_NTSOut))
                if (Output_ELFlag == 1) then
                        allocate(Output_ELData(MaxTimeSteps * MaxSeg,Output_NELOut))
                end if
        
        End SUBROUTINE       

End