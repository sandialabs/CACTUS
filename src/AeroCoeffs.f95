SUBROUTINE AeroCoeffs(nElem,alpha75,alpha5,Re,wPNorm,adotnorm,umach,SectInd,IsBE,CL,CD,CN,CT,CLCirc,CM25) 

    use airfoil
    use dystl
    use pidef

    implicit none

    integer :: SectInd, IsBE, nElem
    real :: alpha75, alpha5, adotnorm, wPNorm, Re, umach, CL, CD, CN, CT, CLCirc
    real :: CLstat75, CLstat5, CDstat75, CLdyn5, CDdyn5, dCLAD, dCTAM, dCNAM, CL5, CD5, C, C1, CM25stat, CM25
    real :: alphaL, alphaD, aref, Fac


    ! Calc static characteristics
    CALL intp(Re,alpha75*condeg,CLstat75,CDstat75,CM25stat,SectInd) 
    CALL intp(Re,alpha5*condeg,CLstat5,C,C1,SectInd)

    ! Apply pitch rate effects by analogy to pitching flat plate potential flow theory (SAND report)
    CL5=CLstat75
    CD5=CDstat75
    CM25=CM25stat+cos(alpha5)*(CLstat75-CLstat5)/4.0
    alphaL=alpha75
    CLCirc=CLstat75

    ! If no dynamic stall, use static values, else calc dynamic stall
    if (DSFlag/=0) then

        if (DSFlag==1) then
            ! Modified Boeing-Vertol approach  
            Call BV_DynStall(nElem,CL5,CD5,alphaL,adotnorm,umach,Re,SectInd,CLdyn5,CDdyn5)   
        else
            ! Leishman-Beddoes model
            Call LB_DynStall(nElem,CL5,CD5,alphaL,alpha5,umach,Re,SectInd,CLdyn5,CDdyn5)
        end if

        CL5=CLdyn5
        CD5=CDdyn5  
        CLCirc=CLdyn5  

    end if

    ! Tangential and normal coeffs
    CN=CL5*cos(alpha5)+CD5*sin(alpha5)                                   
    CT=-CL5*sin(alpha5)+CD5*cos(alpha5) 

    ! Calc tangential added mass increment by analogy to pitching flat plate potential flow theory (SAND report) 
    dCTAM=2.0/cos(alpha5)*wPNorm*CM25stat-CLstat5/2.0*wPNorm
    ! Add in alphadot added mass effects (Theodorsen flat plate approx., Katz ch. 13)
    dCLAD=pi*adotnorm
    dCTAM=dCTAM-dCLAD*sin(alpha5)
    dCNAM=dCLAD*cos(alpha5)       

    ! Add in added mass effects at low AOA (models not accurate at high AOA)
    Fac=1.0
    aref=abs(alpha5)
    if ((aref > pi/4.0) .AND. (aref < 3.0*pi/4.0)) then
        Fac = abs(1-4.0/pi*(aref-pi/4.0))
    end if
    CT=CT+Fac*dCTAM
    CN=CN+Fac*dCNAM

    ! Calc total lift and drag coefficient based on flow direction at half-chord for reference
    CL=CN*cos(alpha5)-CT*sin(alpha5)
    CD=CN*sin(alpha5)+CT*cos(alpha5)

    Return
End SUBROUTINE AeroCoeffs
