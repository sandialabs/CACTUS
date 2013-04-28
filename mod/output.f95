MODULE output

        ! Arrays to be output to csv data files

        ! Output flags
        integer :: Output_ELFlag                ! Set to 1 to output loads data at each timestep for each element, 0 to omit this output
        
        ! Scale parameters and flow properties
        character(10000) :: Output_SFHead = 'R (ft),Frontal Area (ft^2),RPM,U (ft/s),rho (slugs/ft^3),tempr (degF),vis (slugs/(ft*s)),q (lb/ft^2),TSR (-),ReM (-),FnR (-)'
        real :: Output_SFData(1,11)                             ! Data
        
        ! Revolution average data
        character(10000) :: Output_RevHead = 'Rev,Power Coeff. (-),Tip Power Coeff. (-),Torque Coeff. (-),Fx Coeff. (-),Fy Coeff. (-),Fz Coeff. (-),Power (kW),Torque (ft-lbs),Delta CPU Time (s),Total CPU Time (s)'
        real :: Output_RevData(1,11)                ! Revolution average data for each revolution 
  
        ! Timestep data
        character(10000) :: Output_TSHead = 'Normalized Time (-),Theta (rad),Rev,Torque Coeff. (-),Power Coeff. (-),Fx Coeff. (-),Fy Coeff. (-),Fz Coeff. (-)'
        character(1000) :: Output_TSHeadBlade = 'Blade Fx Coeff. (-),Blade Fy Coeff. (-),Blade Fz Coeff. (-),Blade Torque Coeff. (-)'
        character(1000) :: Output_TSHeadStrut = 'Strut Fx Coeff. (-),Strut Fy Coeff. (-),Strut Fz Coeff. (-),Strut Torque Coeff. (-)'
        integer :: Output_NTSOut                                ! Number of timestep outputs
        real, allocatable :: Output_TSData(:,:)                 ! Timestep data for each timestep
        
        ! Element loads data
        character(10000) :: Output_ELHead = 'Normalized Time (-),Theta (rad),Blade,Element,Rev,AOA (deg),Re (-),Mach (-),Ur (-),CN (-),CT (-),Fx (-),Fy (-),Fz (-),te (-),L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12'
        integer :: Output_NELOut = 27                           ! Number of element loads outputs
        integer :: Output_ELRow                                 ! Rows of element row output
        real, allocatable :: Output_ELData(:,:)                 ! Element loads data for each timestep
        
        CONTAINS
        
        SUBROUTINE output_cns(MaxSeg, MaxBlades, MaxStruts)
                
                ! Constructor for the arrays in this module
                
                integer :: MaxSeg, MaxBlades, MaxStruts
                
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
        
        End SUBROUTINE       

End