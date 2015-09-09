SUBROUTINE CalcFreestream(xElem,yElem,zElem,u,v,w,ygcErr)

    use shear
    use configr
    use iecgust
    use tower

    real u, v, w
    real xElem,yElem,zElem
    integer ygcErr

    real IECGustVel, Vxtower

    ! Freestream velocity (with ground shear model)
    ! At y/R = 0, u/Uinf = 0. At y/R = yref, u/Uinf = 1
    ! slex = 0 : Constant freestream
    ! slex = 1/2 : Laminar shear layer (approx)
    ! slex = 1/7 : Turbulent shear layer (approx)

    if (Igust .EQ. 1) then
        u = IECGustVel((nt-1)*dt,xElem)
        v = 0.0
        w = 0.0
        !Write(*,'(F15.8)') u
        Return
    end if

    if ((yElem+ygc) <= 0.0) then
        ! reflect across zero...
        ygcErr=1
        u=(-(yElem+ygc)/yref)**slex
        u=max(u,.01)  ! limit to some non zero value...
        v=0.0
        w=0.0
    else
        u=((yElem+ygc)/yref)**slex
        u=max(u,.01)  ! limit to some non zero value...
        v=0.0
        w=0.0
    end if

    if (Itower .EQ. 1) then
       Vxtower = wake_defect_velocity(xElem,yElem,zElem)
       !Write(20,'(4F20.12)') xElem,yElem,zElem,Vxtower
       u = u - Vxtower
       u = max(u,.01)
       !If (Vxtower .GT. .001) Write(20,'(F20.12)') u+Vxtower, u
    end if

    Return
End SUBROUTINE CalcFreestream


REAL FUNCTION IECGustVel(time,x)

    Use iecgust
    Use pidef

    real time, x

    real tr

    tr = time - x - gustX0

    if ((tr .LE. gustT) .AND. (tr .GT. 0)) then
        IECGustVel = 1.0 - 0.37 * gustA * sin(3*pi*tr/gustT) &
            * (1.0 - cos(2*pi*tr/gustT))
        !Write(*,'(2F20.12)') tr,IECGustVel
    else
        IECGustVel = 1.0
    end if

    Return

End FUNCTION IECGustVel
