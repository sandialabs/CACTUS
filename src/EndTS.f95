subroutine EndTS()

    use util
    use configr
    use output
    use element
    use strut
    use regtest

    Implicit None

    real :: CP, CTR, CFx, CFy, CFz
    integer :: i, offset

    ! Collect timestep results and compile output

    ! Machine output
    CP=CP_B+CP_S
    CTR=CTR_B+CTR_S
    CFx=CFx_B+CFx_S
    CFy=CFy_B+CFy_S
    CFz=CFz_B+CFz_S

    ! Sums for rev averages
    CPSum=CPSum+CP
    CTRSum=CTRSum+CTR
    CFxSum=CFxSum+CFx
    CFySum=CFySum+CFy
    CFzSum=CFzSum+CFz

    TSOutData(1,1)=TimeN         ! Normalized simulation time (t*Uinf/Rmax)
    TSOutData(1,2)=Theta         ! Turbine phase angle (rad)
    TSOutData(1,3)=irev
    TSOutData(1,4)=CTR               ! Torque coeff
    TSOutData(1,5)=CP                ! Power coeff
    TSOutData(1,6)=CFx                ! Fx coeff
    TSOutData(1,7)=CFy                ! Fy coeff
    TSOutData(1,8)=CFz                ! Fz coeff
    ! Blade output
    do i=1,nb
        offset=8+(i-1)*4
        TSOutData(1,offset+1)=Blades(i)%CFx     ! Blade Fx coeff
        TSOutData(1,offset+2)=Blades(i)%CFy     ! Blade Fy coeff
        TSOutData(1,offset+3)=Blades(i)%CFz     ! Blade Fz coeff
        TSOutData(1,offset+4)=Blades(i)%CTR     ! Blade torque coeff
    end do
    do i=1,NStrut
        offset=8+nb*4+(i-1)*4
        TSOutData(1,offset+1)=Struts(i)%CFx     ! Strut Fx coeff
        TSOutData(1,offset+2)=Struts(i)%CFy     ! Strut Fy coeff
        TSOutData(1,offset+3)=Struts(i)%CFz     ! Strut Fz coeff
        TSOutData(1,offset+4)=Struts(i)%CTR     ! Strut torque coeff
    end do

    ! Write to timestep data csv file
    Call csvwrite(10,TSOutHead,TSOutData,0,1)

    ! Write to element loads data csv file
    if (BladeElemOutFlag == 1) then
        Call csvwrite(11,BladeElemOutHead,BladeElemOutData,0,-1)
    end if

    ! Write to dynamic stall diagnostic data csv file
    if (DynStallOutFlag == 1) then
        if (DynStallOutType == 1) then
            Call csvwrite(16,DynStallOutBVHead,DynStallOutBVData,0,-1)
        else if (DynStallOutType == 2) then
            Call csvwrite(16,DynStallOutLBHead,DynStallOutLBData,0,-1)
        end if
    end if

    ! Reg test
    if (RegTFlag == 1) then
        Reg_CPOut=CP
    end if

    return
end subroutine EndTS
