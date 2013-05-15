SUBROUTINE BladeLoads(NLTol,iConv)  

        use parameters
        use pidef
	use blade            
	use configr
        use regtest
        use output       
        use element
        use airfoil
        use dystl

        Implicit None
	
        integer iConv
        integer i, j, nei, nej, nej1, IsBE, Loop, LBCheck
        real alpha, Re, umach, ur, CN, CT, te, NLTol, dgb, Fx, Fy, Fz
        real CTExcr
        
        ! Calculates blade performance, bound and new shed vorticity        
        
        
        ! Zero out current blades loads sum
        CP_B=0.0
        CTR_B=0.0
        CFx_B=0.0
        CFy_B=0.0
        CFz_B=0.0 
       
	iConv=0											                                           
	do i=1,nb
                                                                                                             
                ! Zero out current blade loads
                Blades(i)%CP=0.0
                Blades(i)%CTR=0.0
                Blades(i)%CFx=0.0
                Blades(i)%CFy=0.0
                Blades(i)%CFz=0.0                                                                                             
                                                                                                             
		nei=1+(i-1)*(nbe+1)                                               
		do j=1,nbe                                                     
			nej=nei+j                                                         
			nej1=nej-1 
			
			IsBE=0
			if (j==1 .OR. j==nbe) then
				IsBE=1
			end if
			                                                                                                                             
			! Calculate the loads on the blade segment                                                                                              
			CALL bsload(nej,IsBE,alpha,Re,umach,ur,CN,CT,Fx,Fy,Fz,te) 
			                                                                                  
                        if (TSFilFlag == 1) then
                            ! Don't allow NL iteration to update bound vorticity (fixed at filtered value)
                            ! This trivializes the NL iteration...
                            dgb=0.0                                                                       
                        else
                            ! Set GB directly equal to GB_Raw (no filtering)
                            GB(nej1)=GB_Raw(nej1)                                                                      
                                                                                                  
                            ! Calculate the bound vortex strength change (w.r.t reference circulation)                                                                                                   
                            dgb=abs((GB(nej1)-GS(nt,nej1))/(CrRef*ut))                          
                            
                            ! If change outside tolerance for any element, set flag
                            if (dgb .gt. NLTol) iConv=1 
                            
                            ! Reset the bound circulation as the current entry in the spanwise
                            ! vorticity array for velocity calculation (and eventual wake convection)                                  
                            GS(nt,nej1)=GB(nej1)  
                        end if                     
                         
                        
                        ! Update output data 
                        
                        ! Calc LB dynamic stall model logic checksum
                        Call LB_LogicChecksum(nej,LBCheck)                          
                                                   
                        ! Element loads output  
                        ! JCM: when blade/element module is reorg'd, this data will be held for each element as part of the 
                        ! blades structure and sampled for output in EndTS.                                                  
                        if (Output_ELFlag == 1) then
                                Output_ELRow=(i-1)*nbe+j
                                Output_ELData(Output_ELRow,1)=TimeN         ! Normalized simulation time (t*Uinf/Rmax) 
                                Output_ELData(Output_ELRow,2)=Theta         ! Phase angle
                                Output_ELData(Output_ELRow,3)=i      
                                Output_ELData(Output_ELRow,4)=j 
                                Output_ELData(Output_ELRow,5)=irev
                                Output_ELData(Output_ELRow,6)=alpha*condeg              ! Element angle of attack
                                Output_ELData(Output_ELRow,7)=Re                        ! Element Reynolds number based on local chord and flow velocity
                                Output_ELData(Output_ELRow,8)=umach                     ! Element Mach number based on local flow velocity
                                Output_ELData(Output_ELRow,9)=ur                        ! Element velocity ratio with freestream
                                Output_ELData(Output_ELRow,10)=CN                        ! Element normal force coefficient (per span) based on local chord and flow velocity
                                Output_ELData(Output_ELRow,11)=CT                       ! Element tangential force coefficient (per span) based on local chord and flow velocity
                                Output_ELData(Output_ELRow,12)=Fx                       ! Element global x force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,13)=Fy                       ! Element global y force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,14)=Fz                       ! Element global z force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,15)=te                       ! Element torque coefficient contribution based on freestream flow, turbine area, and Rmax       
                                
                                ! BV Logic
                                Output_ELData(Output_ELRow,16)=BVLogicOutputs(1)          ! BV Dynamic stall flag for lift coefficient
                                Output_ELData(Output_ELRow,17)=BVLogicOutputs(2)          ! BV Dynamic stall flag for drag coefficient
                                
                                ! LB Logic 
                                do Loop=1,9
                                        Output_ELData(Output_ELRow,17+Loop)=LBLogicOutputs(nej,Loop)   
                                end do
                                Output_ELData(Output_ELRow,27)=LBCheck           
                                
                        end if

                        
                        ! Add to blade output
                        Blades(i)%CTR=Blades(i)%CTR + te        ! Torque coeff from this blade, based on freestream flow, turbine area, and Rmax
                        Blades(i)%CP=Blades(i)%CP + te*ut       ! Power coeff from this blade, based on freestream flow, turbine area, and Rmax
                        Blades(i)%CFx=Blades(i)%CFx + Fx        ! Force coeff along global x on this blade, based on freestream flow and turbine area
                        Blades(i)%CFy=Blades(i)%CFy + Fy        ! Force coeff along global y on this blade, based on freestream flow and turbine area
                        Blades(i)%CFz=Blades(i)%CFz + Fz        ! Force coeff along global z on this blade, based on freestream flow and turbine area
                        
                          
                        ! Regression test
                        if (RegTFlag == 1) then
                                Reg_ElemNum=nej1 
                                Reg_DFL=BVLogicOutputs(1)
                                Reg_LBC=LBCheck
                                Reg_ElemAOA=alpha*180.0/3.14159 
                                Reg_ElemCirc=GB(nej1)
                                Reg_dElemCirc=dgb   
                                Call WriteRegTOutput(1)                  
                        end if                                                
                                                                  
		end do
                
                ! Add to total blades output
                CTR_B=CTR_B + Blades(i)%CTR
                CP_B=CP_B + Blades(i)%CP
                CFx_B=CFx_B + Blades(i)%CFx
                CFy_B=CFy_B + Blades(i)%CFy
                CFz_B=CFz_B + Blades(i)%CFz


	end do
                                                                 
        ! Apply any user specified machine level excrescence torque. CTExcrM = TorqueExcr / (1/2*rho*Utip^2*Rmax^3)
        CTExcr = CTExcrM*ut**2/at 
        CTR_B=CTR_B + - CTExcr
        CP_B=CP_B - CTExcr*ut                                                   

Return
End                                                               
