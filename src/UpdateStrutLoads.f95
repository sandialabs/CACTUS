SUBROUTINE UpdateStrutLoads()  

    use strut
    use configr

    Implicit None

    integer :: i, j, NElem
    integer :: ygcerr
    real :: xs, ys, zs, xj, yj, zj
    real :: uFSs, vFSs, wFSs, uBlade, vBlade, wBlade, us, vs, ws
    real :: uTot, vTot, wTot, ur, ReStrut
    real :: Delem
    real :: Fx, Fy, Fz, TRx, TRy, TRz, te
    real :: t_ave, carea, cdj, djunc

    ! Updates loads on struts

    ! Zero out current strut loads sum
    CP_S=0.0
    CTR_S=0.0
    CFx_S=0.0
    CFy_S=0.0
    CFz_S=0.0                          

    do i=1,NStrut

        ! Zero out current strut loads
        Struts(i)%CP=0.0
        Struts(i)%CTR=0.0
        Struts(i)%CFx=0.0
        Struts(i)%CFy=0.0
        Struts(i)%CFz=0.0

        NElem=Struts(i)%NElem
        do j=1,NElem 

            xs=Struts(i)%PEx(j)
            ys=Struts(i)%PEy(j)
            zs=Struts(i)%PEz(j)

            ! Freestream velocity at strut location
            Call CalcFreestream(xs,ys,zs,uFSs,vFSs,wFSs,ygcerr)

            ! Blade velocity due to rotation                                                      
            CALL CalcBladeVel(wRotX,wRotY,wRotZ,xs,ys,zs,uBlade,vBlade,wBlade)

            ! Induced velocity at strut element
            Call CalcIndVel(NT,ntTerm,NBE,NB,NE,xs,ys,zs,us,vs,ws)

            ! Calculate relative velocity magnitude at strut element
            uTot = us+uFSs-uBlade
            vTot = vs+vFSs-vBlade
            wTot = ws+wFSs-wBlade
            ur = sqrt(uTot*uTot + vTot*vTot + wTot*wTot)
            ReStrut = ReM*ur*Struts(i)%ECtoR(j)  ! Strut chord Reynolds number

            ! Fill current flow quantities at strut element
            Struts(i)%u(j)=uTot
            Struts(i)%v(j)=vTot
            Struts(i)%w(j)=wTot
            Struts(i)%ur(j)=ur
            Struts(i)%ReStrut(j)=ReStrut

            ! Calculate strut element coeffs
            Call StrutElemCoeffs(i,j)

            ! Drag coeff vector from this strut element, re-referenced to full turbine scale 
            ! (D/(1/2*rho*Uinf^2*At)
            Delem = Struts(i)%Cd0(j) * Struts(i)%EAreaR(j) / at * ur**2
            Fx=Delem*uTot/ur
            Fy=Delem*vTot/ur
            Fz=Delem*wTot/ur
            ! Corresponding torque coeff. (T/(1/2*rho*Uinf^2*At*R)) 
            CALL cross(xs-RotPX,ys-RotPY,zs-RotPZ,Fx,Fy,Fz,TRx,TRy,TRz)
            te=(TRx*RotX+TRy*RotY+TRz*RotZ)

            ! Add to strut output
            Struts(i)%CTR=Struts(i)%CTR + te
            Struts(i)%CP=Struts(i)%CP + te*ut
            Struts(i)%CFx=Struts(i)%CFx + Fx
            Struts(i)%CFy=Struts(i)%CFy + Fy
            Struts(i)%CFz=Struts(i)%CFz + Fz
        end do

        ! Blade/strut junction interference drag
        ! First strut element
        if (Struts(i)%BIndS > 0) then
            carea=Struts(i)%ECtoR(1)**2
            xj=Struts(i)%PEx(1)
            yj=Struts(i)%PEy(1)
            zj=Struts(i)%PEz(1)
            t_ave = 0.5 * (Struts(i)%TtoC + Struts(i)%tcS)
            uTot=Struts(i)%u(1)
            vTot=Struts(i)%v(1)
            wTot=Struts(i)%w(1)
            ur=Struts(i)%ur(1)
            Cdj = t_ave*t_ave * (17.0 * t_ave*t_ave - 0.05)
            !Cdj = 0.0112  ! t/c_avg = 0.165
            !Cdj = 0.0535  ! t/c_avg = 0.24
            Cdj = Cdj + Cdpar  ! Additional user-specified parasitic drag
            ! Drag coeff vector, re-referenced to full turbine scale 
            ! (D/(1/2*rho*Uinf^2*At)
            Djunc = Cdj * carea / at * ur**2
            Fx=Djunc*uTot/ur
            Fy=Djunc*vTot/ur
            Fz=Djunc*wTot/ur
            ! Corresponding torque coeff. (T/(1/2*rho*Uinf^2*At*R)) 
            CALL cross(xj-RotPX,yj-RotPY,zj-RotPZ,Fx,Fy,Fz,TRx,TRy,TRz)
            te=(TRx*RotX+TRy*RotY+TRz*RotZ)
            ! Add to strut output
            Struts(i)%CTR=Struts(i)%CTR + te
            Struts(i)%CP=Struts(i)%CP + te*ut
            Struts(i)%CFx=Struts(i)%CFx + Fx
            Struts(i)%CFy=Struts(i)%CFy + Fy
            Struts(i)%CFz=Struts(i)%CFz + Fz
        end if
        ! Last strut element
        if (Struts(i)%BIndE > 0) then
            carea=Struts(i)%ECtoR(NElem)**2
            xj=Struts(i)%PEx(NElem)
            yj=Struts(i)%PEy(NElem)
            zj=Struts(i)%PEz(NElem)
            t_ave = 0.5 * (Struts(i)%TtoC + Struts(i)%tcE)
            uTot=Struts(i)%u(NElem)
            vTot=Struts(i)%v(NElem)
            wTot=Struts(i)%w(NElem)
            ur=Struts(i)%ur(NElem)
            Cdj = t_ave*t_ave * (17.0 * t_ave*t_ave - 0.05)
            !Cdj = 0.0112  ! t/c_avg = 0.165
            !Cdj = 0.0535  ! t/c_avg = 0.24
            Cdj = Cdj + Cdpar  ! Additional user-specified parasitic drag
            ! Drag coeff vector, re-referenced to full turbine scale 
            ! (D/(1/2*rho*Uinf^2*At)
            Djunc = Cdj * carea / at * ur**2
            Fx=Djunc*uTot/ur
            Fy=Djunc*vTot/ur
            Fz=Djunc*wTot/ur
            ! Corresponding torque coeff. (T/(1/2*rho*Uinf^2*At*R)) 
            CALL cross(xj-RotPX,yj-RotPY,zj-RotPZ,Fx,Fy,Fz,TRx,TRy,TRz)
            te=(TRx*RotX+TRy*RotY+TRz*RotZ)
            ! Add to strut output
            Struts(i)%CTR=Struts(i)%CTR + te
            Struts(i)%CP=Struts(i)%CP + te*ut
            Struts(i)%CFx=Struts(i)%CFx + Fx
            Struts(i)%CFy=Struts(i)%CFy + Fy
            Struts(i)%CFz=Struts(i)%CFz + Fz
        end if

        ! Add to total struts output
        CTR_S=CTR_S + Struts(i)%CTR
        CP_S=CP_S + Struts(i)%CP
        CFx_S=CFx_S + Struts(i)%CFx
        CFy_S=CFy_S + Struts(i)%CFy
        CFz_S=CFz_S + Struts(i)%CFz

    end do

    Return
End SUBROUTINE UpdateStrutLoads
