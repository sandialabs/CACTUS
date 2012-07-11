SUBROUTINE BladeLoads(nGeom,NLTol,iConv)  

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
	
        integer nGeom, iConv
        integer i, j, nei, nej, nej1, IsBE, offset, Loop, LBCheck
        real alpha, Re, umach, ur, CN, CT, te, NLTol, dgb, Fx, Fy, Fz
        real BladeLoad(MaxBlades,3), BladeTorque(MaxBlades)
        real CTExcr
	real cp, ctr
        
	real xe,ye,ze,rade,dr,uFSs,vFSs,wFSs,uBlade,vBlade,wBlade,us,vs,ws
	real uTot,vTot,wTot,Delem,Dtorq,Cd0,Cdj,t_ave,Djunc,restrut,cflam,cfturb
        real cdlam,cdturb,fblend
        real, parameter :: recrit=3.0e5
        integer ygcerr 
        
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
			CALL bsload(nej,nGeom,IsBE,alpha,Re,umach,ur,CN,CT,Fx,Fy,Fz,te) 
			                                                                                  
			! Calculate the bound vortex strength change (w.r.t reference circulation)                                                                                                   
			dgb=abs((GB(nej1)-GS(nt,nej1))/(CrRef*ut))                          
			
			! If change outside tolerance for any element, set flag
			if (dgb .gt. NLTol) iConv=1 
			 
                        ! Set the bound circulation as the current entry in the spanwise
                        ! vorticity array for velocity calculation (and eventual wake convection)                                  
                        GS(nt,nej1)=GB(nej1)                       
                         
                        
                        ! Update output data 
                        
                        ! Calc LB dynamic stall model logic checksum
                        Call LB_LogicChecksum(nej,LBCheck)                          
                                                   
                        ! Element loads output                                                     
                        if (Output_ELFlag == 1) then
                                Output_ELRow=(nt-1)*(nb*nbe)+(i-1)*nbe+j
                                Output_ELData(Output_ELRow,1)=(nt-1)*DT         ! Normalized simulation time (t*Uinf/Rmax) 
                                Output_ELData(Output_ELRow,2)=i      
                                Output_ELData(Output_ELRow,3)=j 
                                Output_ELData(Output_ELRow,4)=irev
                                Output_ELData(Output_ELRow,5)=alpha*condeg              ! Element angle of attack
                                Output_ELData(Output_ELRow,6)=Re                        ! Element Reynolds number based on local chord and flow velocity
                                Output_ELData(Output_ELRow,7)=umach                     ! Element Mach number based on local flow velocity
                                Output_ELData(Output_ELRow,8)=ur                        ! Element velocity ratio with freestream
                                Output_ELData(Output_ELRow,9)=CN                        ! Element normal force coefficient (per span) based on local chord and flow velocity
                                Output_ELData(Output_ELRow,10)=CT                       ! Element tangential force coefficient (per span) based on local chord and flow velocity
                                Output_ELData(Output_ELRow,11)=Fx                       ! Element global x force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,12)=Fy                       ! Element global y force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,13)=Fz                       ! Element global z force coefficient based on freestream flow and turbine area
                                Output_ELData(Output_ELRow,14)=te                       ! Element torque coefficient contribution based on freestream flow, turbine area, and Rmax       
                                
                                ! BV Logic
                                Output_ELData(Output_ELRow,15)=BVLogicOutputs(1)          ! BV Dynamic stall flag for lift coefficient
                                Output_ELData(Output_ELRow,16)=BVLogicOutputs(2)          ! BV Dynamic stall flag for drag coefficient
                                
                                ! LB Logic 
                                do Loop=1,9
                                        Output_ELData(Output_ELRow,16+Loop)=LBLogicOutputs(nej,Loop)   
                                end do
                                Output_ELData(Output_ELRow,26)=LBCheck           
                                
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
                                Reg_DFL=BVLogicOutputs(1)
                                Reg_LBC=LBCheck
                                Reg_ElemAOA=alpha*180.0/3.14159 
                                Reg_ElemCirc=GB(nej1)
                                Reg_dElemCirc=dgb   
                                Call WriteRegTOutput(1)                  
                        end if                                                
                                                                  
		end do

                ! Calculate torque and power loss due to strut drag
                If (Istrut == 1) Then 
                  Do j=1,nbe
                     nej=nei+j
				
                     ! Strut element center location
                     xe=0.5*(xSe(nGeom,nej)+xSe(nGeom,nej-1))
                     ye=0.5*(ySe(nGeom,nej)+ySe(nGeom,nej-1))
                     ze=0.5*(zSe(nGeom,nej)+zSe(nGeom,nej-1))

                     rade = 0.5*(real(j)/real(nbe) + real(j-1)/real(nbe)) ! element radius 
                     dr = 1.0 / real(nbe) ! element width

                     ! Freestream velocity at strut location
                     Call CalcFreestream(ySe(nGeom,nej),uFSs,vFSs,wFSs,ygcerr)

                     ! Blade velocity due to rotation                                                      
                     CALL CalcBladeVel(wRotX,wRotY,wRotZ,xe,ye,ze,uBlade,vBlade,wBlade)

                     ! Induced velocity at strut element
                     Call CalcIndVel(NT,ntTerm,NBE,NB,NE,xe,ye,ze,us,vs,ws)

                     ! Calculate relative velocity magnitude at strut element
                     uTot = us+uFSs-uBlade
                     vTot = vs+vFSs-vBlade
                     wTot = ws+wFSs-wBlade
                     ur = sqrt(uTot*uTot + vTot*vTot + wTot*wTot)

                     ! Calculate strut element profile drag
                     Restrut = ReM*ur*eChord(nbe/2)  ! Strut chord Reynolds number
                     Cflam = 2.66 / SQRT(Restrut)  ! Laminar friction drag coefficient
                     Cdlam = 2.0 * Cflam * (1 + sthick) + sthick*sthick ! Laminar drag coefficient
                     Cfturb = 0.044 / Restrut**(1.0/6.0) ! Turbulent friction drag coefficient
                     Cdturb = 2.0 * Cfturb * (1.0 + 2.0*sthick + 60*sthick*sthick*sthick*sthick) ! Turbulent drag coefficient
                     Fblend = 0.5 * (1.0 + TANH((LOG10(Restrut)-LOG10(Recrit))/0.2)) ! Blending function for transition between laminar and turbulent drag 
                     Cd0 = (1.0-Fblend) * Cdlam + Fblend * Cdturb ! Profile drag coefficient

                     Delem = Cd0 * dr * eChord(nbe/2) / at * ur*ur
                     Dtorq = Delem * rade
                     BladeTorque(i)=BladeTorque(i)-Dtorq
                     ctr = ctr - Dtorq
                     cp  = cp - Dtorq*ut
                  End Do

                  ! Blade/strut junction interference drag
                  t_ave = 0.5 * (sthick + tc(iSect(nbe/2)))
                  Cdj = t_ave*t_ave * (17.0 * t_ave*t_ave - 0.05)
                  !Cdj = 0.0112  ! t/c_avg = 0.165
                  !Cdj = 0.0535  ! t/c_avg = 0.24
                  Cdj = Cdj + Cdpar  ! Additional user-specified parasitic drag
                  Djunc = Cdj * eChord(nbe/2) * eChord(nbe/2) / at * ur*ur
                  Dtorq = Djunc ! (* 1.0) junction drag acts at r=1
                  BladeTorque(i)=BladeTorque(i)-Dtorq
                  ctr = ctr - Dtorq
                  cp  = cp - Dtorq*ut
               End If
	end do
                                                                 
        ! Apply any user specified machine level excrescence torque. CTExcrM = TorqueExcr / (1/2*rho*Utip^2*Rmax^3)
        CTExcr = CTExcrM*ut**2/at 
        ctr = ctr - CTExcr 
        cp  = cp - CTExcr*ut                                                     
                                                                 
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
