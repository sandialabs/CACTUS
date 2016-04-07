subroutine BladeLoads(NLTol,iConv)

    use parameters
    use pidef
    use blade
    use wake
    use configr
    use regtest
    use output
    use element
    use airfoil
    use dystl

    Implicit None

    integer iConv
    integer i, j, nei, nej, nej1, IsBE, Loop, LBCheck, nElem
    real alpha, alpha5, alpha75, adotnorm, Re, umach, ur, CL, CD, CM25, CLCirc, CN, CT, te, NLTol, dgb, Fx, Fy, Fz
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
            CALL bsload(nej,IsBE,alpha,alpha5,alpha75,adotnorm,Re,umach,ur,CL,CD,CM25,CLCirc,CN,CT,Fx,Fy,Fz,te)

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

            if (BladeElemOutFlag == 1) then
                BladeElemOutRow=(i-1)*nbe+j
                BladeElemOutData(BladeElemOutRow,1)=TimeN                     ! Normalized simulation time (t*Uinf/Rmax)
                BladeElemOutData(BladeElemOutRow,2)=Theta                     ! Phase angle
                BladeElemOutData(BladeElemOutRow,3)=i                         ! Blade Number
                BladeElemOutData(BladeElemOutRow,4)=j                         ! Element number
                BladeElemOutData(BladeElemOutRow,5)=irev                      ! rotation number
                BladeElemOutData(BladeElemOutRow,6)=xBC(nej)                  ! Element quarter chord position x
                BladeElemOutData(BladeElemOutRow,7)=yBC(nej)                  ! Element quarter chord position y
                BladeElemOutData(BladeElemOutRow,8)=zBC(nej)                  ! Element quarter chord position z
                BladeElemOutData(BladeElemOutRow,9)=alpha*condeg              ! Element angle of attack @ 25% chord
                BladeElemOutData(BladeElemOutRow,10)=alpha5*condeg            ! Element angle of attack @ 50% chord
                BladeElemOutData(BladeElemOutRow,11)=alpha75*condeg           ! Element angle of attack @ 75% chord
                BladeElemOutData(BladeElemOutRow,12)=adotnorm                 ! Normalized AOA rate
                BladeElemOutData(BladeElemOutRow,13)=Re                       ! Element Reynolds number based on local chord and flow velocity
                BladeElemOutData(BladeElemOutRow,14)=umach                    ! Element Mach number based on local flow velocity
                BladeElemOutData(BladeElemOutRow,15)=ur                       ! Element velocity ratio with freestream
                BladeElemOutData(BladeElemOutRow,16)=(UB(nej)+UB(nej-1))/2.0  ! Element induced x velocity (average of the endpoint induced velocities)
                BladeElemOutData(BladeElemOutRow,17)=(VB(nej)+VB(nej-1))/2.0  ! Element induced y velocity (average of the endpoint induced velocities)
                BladeElemOutData(BladeElemOutRow,18)=(WB(nej)+WB(nej-1))/2.0  ! Element induced z velocity (average of the endpoint induced velocities)
                BladeElemOutData(BladeElemOutRow,19)=(GB(nej1))               ! Element bound vorticity (note the last element is just a placeholder)
                BladeElemOutData(BladeElemOutRow,20)=CL                       ! Element lift coeff (defined with alpha5 flow direction)
                BladeElemOutData(BladeElemOutRow,21)=CD                       ! Element drag coeff (defined with alpha5 flow direction)
                BladeElemOutData(BladeElemOutRow,22)=CM25                     ! Element moment coeff (about quarter-chord point)
                BladeElemOutData(BladeElemOutRow,23)=CLCirc                   ! Element circulatory lift coeff (defined with alpha5 flow direction and specifies element bound circ strength)
                BladeElemOutData(BladeElemOutRow,24)=CN                       ! Element normal force coefficient (per span) based on local chord and flow velocity
                BladeElemOutData(BladeElemOutRow,25)=CT                       ! Element tangential force coefficient (per span) based on local chord and flow velocity
                BladeElemOutData(BladeElemOutRow,26)=Fx                       ! Element global x force coefficient based on freestream flow and turbine area
                BladeElemOutData(BladeElemOutRow,27)=Fy                       ! Element global y force coefficient based on freestream flow and turbine area
                BladeElemOutData(BladeElemOutRow,28)=Fz                       ! Element global z force coefficient based on freestream flow and turbine area
                BladeElemOutData(BladeElemOutRow,29)=te                       ! Element torque coefficient contribution based on freestream flow, turbine area, and Rmax
            end if


            ! Dynamic stall diagnostic output
            if (DynStallOutFlag == 1 .AND. DynStallOutType == 1) then
                DynStallOutBVRow=(i-1)*nbe+j
                DynStallOutBVData(DynStallOutBVRow,1)=TimeN         ! Normalized simulation time (t*Uinf/Rmax)
                DynStallOutBVData(DynStallOutBVRow,2)=Theta         ! Phase angle
                DynStallOutBVData(DynStallOutBVRow,3)=i
                DynStallOutBVData(DynStallOutBVRow,4)=j
                DynStallOutBVData(DynStallOutBVRow,5)=irev
                DynStallOutBVData(DynStallOutBVRow,6)=BV_alpha
                DynStallOutBVData(DynStallOutBVRow,7)=BV_adotnorm
                DynStallOutBVData(DynStallOutBVRow,8)=BV_alrefL                   ! ref lift AOA
                DynStallOutBVData(DynStallOutBVRow,9)=BV_alrefD                   ! ref drag AOA
                DynStallOutBVData(DynStallOutBVRow,10)=BV_DynamicFlagL(nej)       ! Lift flag
                DynStallOutBVData(DynStallOutBVRow,11)=BV_DynamicFlagD(nej)       ! Drag flag
            else if (DynStallOutFlag == 1 .AND. DynStallOutType == 2) then
                DynStallOutLBRow=(i-1)*nbe+j
                DynStallOutLBData(DynStallOutLBRow,1)=TimeN         ! Normalized simulation time (t*Uinf/Rmax)
                DynStallOutLBData(DynStallOutLBRow,2)=Theta         ! Phase angle
                DynStallOutLBData(DynStallOutLBRow,3)=i
                DynStallOutLBData(DynStallOutLBRow,4)=j
                DynStallOutLBData(DynStallOutLBRow,5)=irev
                ! LB Logic
                do Loop=1,9
                    DynStallOutLBData(DynStallOutLBRow,5+Loop)=LB_LogicOutputs(nej,Loop)
                end do
                BladeElemOutData(BladeElemOutRow,15)=LBCheck
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
                Reg_DFL=BV_DynamicFlagL(nej)
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
    CTR_B=CTR_B - CTExcr
    ! CTR_B=CTR_B + - CTExcr
    CP_B=CP_B - CTExcr*ut

    return
end subroutine BladeLoads
