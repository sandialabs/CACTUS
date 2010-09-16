SUBROUTINE bvort(nGeom,NLTol,iConv)  

        use parameters
        use pidef
	use blade            
	use test
	use configr
        use regtest
        use output       
	
        integer IsBE, DynamicFlagL, DynamicFlagD, offset
        real alpha, Re, umach, ur, CN, CT, te, NLTol, dgb, cpl, Fx, Fy, Fz
        real BladeLoad(MaxBlades,3), BladeTorque(MaxBlades)
	real cp, ctr
        
        cp=0.0
        ctr=0.0
        BladeLoad(:,:)=0.0
        BladeTorque(:)=0.0
        
	! Calculates blade performance, bound and new shed vorticity
	iConv=0											                                           
	do i=1,nb                                                                                                     
		nei=1+(i-1)*(nbe+1)                                               
		do j=1,nbe                                                     
			nej=nei+j                                                         
			nej1=nej-1 
			
			IsBE=0
			if (j==1 .OR. j==nbe) then
				IsBE=1
			end if
			                                                                                                                             
			! Calculate the loads on the blade segment                                                                                              
			CALL bsload(nej,nGeom,IsBE,DynamicFlagL,DynamicFlagD,alpha,Re,umach,ur,CN,CT,Fx,Fy,Fz,te) 
			                                                                                  
			! Calculate the bound vortex strength change                                                                                                   
			dgb=abs((GB(nej1)-GS(nt,nej1))/GB(nej1))                          
			
			! If change outside tolerance for any element, set flag
			if (dgb .gt. NLTol) iConv=1 
			 
                        ! Set the bound circulation as the current entry in the spanwise
                        ! vorticity array for velocity calculation (and eventual wake convection)                                  
                        GS(nt,nej1)=GB(nej1)                       
                         
                        
                        ! Update output data 
                                                   
                        ! Element loads output                                                     
                        if (Output_ELFlag == 1) then
                                Output_ELRow=(nt-1)*(nb*nbe)+(i-1)*nbe+j
                                Output_ELData(Output_ELRow,1)=(nt-1)*DT         ! Normalized simulation time (t*Uinf/Rmax) 
                                Output_ELData(Output_ELRow,2)=i      
                                Output_ELData(Output_ELRow,3)=j 
                                Output_ELData(Output_ELRow,4)=irev
                                Output_ELData(Output_ELRow,5)=DynamicFlagL      ! Dynamic stall flag for lift coefficient
                                Output_ELData(Output_ELRow,6)=DynamicFlagD      ! Dynamic stall flag for drag coefficient 
                                Output_ELData(Output_ELRow,7)=alpha*condeg      ! Element angle of attack
                                Output_ELData(Output_ELRow,8)=Re                ! Element Reynolds number based on local chord and flow velocity
                                Output_ELData(Output_ELRow,9)=umach             ! Element Mach number based on local flow velocity
                                Output_ELData(Output_ELRow,10)=ur               ! Element velocity ratio with freestream
                                Output_ELData(Output_ELRow,11)=CN               ! Element normal force coefficient (per span) based on local chord and flow velocity
                                Output_ELData(Output_ELRow,12)=CT               ! Element tangential force coefficient (per span) based on local chord and flow velocity
                                Output_ELData(Output_ELRow,13)=Fx               ! Element global x force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,14)=Fy               ! Element global y force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,15)=Fz               ! Element global z force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,16)=te               ! Element torque coefficient contribution based on freestream flow, turbine area, and Rmax              
                        end if

                        ! Machine level output
                        ctr=ctr+te        ! Torque coeff based on freestream speed and Rmax     
                        cp=cp+te*ut       ! Power coeff based on freestream speed and Rmax 
                        ! Blade level output
                        BladeTorque(i)=BladeTorque(i)+te        ! Torque coeff from this blade, based on freestream flow, turbine area, and Rmax
                        BladeLoad(i,1)=BladeLoad(i,1)+Fx        ! Force coeff along global x on this blade, based on freestream flow and turbine area
                        BladeLoad(i,2)=BladeLoad(i,2)+Fy        ! Force coeff along global y on this blade, based on freestream flow and turbine area
                        BladeLoad(i,3)=BladeLoad(i,3)+Fz        ! Force coeff along global z on this blade, based on freestream flow and turbine area     
                                                
                        ! Regression test
                        if (RegTFlag == 1) then
                                Reg_ElemNum=nej1 
                                Reg_DFL=DynamicFlagL
                                Reg_DFD=DynamicFlagD
                                Reg_ElemAOA=alpha*180.0/3.14159 
                                Reg_ElemCirc=GB(nej1)
                                Reg_dElemCirc=dgb   
                                Call WriteRegTOutput(1)                  
                        end if                                                
                                                                  
		end do 
	end do
                                                                 
        ! Machine level output
        Output_TSRow=nt                  
        Output_TSData(Output_TSRow,1)=(nt-1)*DT         ! Normalized simulation time (t*Uinf/Rmax)
        Output_TSData(Output_TSRow,2)=irev
        Output_TSData(Output_TSRow,3)=ctr               ! Torque coeff
        Output_TSData(Output_TSRow,4)=cp                ! Power coeff
        ! Blade level output
        do i=1,nb
                offset=4+(i-1)*4
                Output_TSData(Output_TSRow,offset+1)=BladeLoad(i,1)     ! Blade Fx coeff
                Output_TSData(Output_TSRow,offset+2)=BladeLoad(i,2)     ! Blade Fy coeff
                Output_TSData(Output_TSRow,offset+3)=BladeLoad(i,3)     ! Blade Fz coeff  
                Output_TSData(Output_TSRow,offset+4)=BladeTorque(i)     ! Blade torque coeff
        end do                                                        
	
        if (RegTFlag == 1) then
                Reg_CPOut=cp      
        end if
        
Return
End                                                               
