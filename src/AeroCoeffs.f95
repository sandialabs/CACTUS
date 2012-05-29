SUBROUTINE AeroCoeffs(nElem,alpha75,alpha5,alpha,Re75,Re5,Re,wPNorm,adotnorm,umach,SectInd,IsBE,CL,CD,CN,CT,CLCirc,CM25) 

        use airfoil
        use dystl
        use pidef

        implicit none
        
        integer :: SectInd, IsBE, nElem
        real :: alpha75, alpha5, alpha, adotnorm, wPNorm, Re75, Re5, Re, umach, CL, CD, CN, CT, CLCirc
        real :: CLstat75, CLstat5, CDstat5, CLdyn5, CDdyn5, dCLAD, dCTAM, dCNAM, CL5, CD5, C, C1, CM25stat, CM25
        real :: alphaL, alphaD, ReL, ReD, aref, Fac
        
        
        ! Calc static characteristics
        CALL intp(Re75,alpha75*condeg,CLstat75,C,C1,SectInd) 
        CALL intp(Re5,alpha5*condeg,CLstat5,CDstat5,CM25stat,SectInd)
        
        ! Apply pitch rate effects by analogy to pitching flat plate potential flow theory 
        ! (SAND81-7017) 
        CL5=CLstat75
        CD5=CDstat5
        CM25=CM25stat+cos(alpha5)*(CLstat75-CLstat5)/4.0
        alphaL=alpha75
        alphaD=alpha5
        ReL=Re75
        ReD=Re5
        CLCirc=CLstat75
        
        ! If blade end segment or no dynamic stall, use static values, else calc dynamic stall
        if (IsBE == 0 .AND. DSFlag/=0) then
            
            if (DSFlag==1) then
                ! Modified Boeing-Vertol approach
                Call BV_DynStall(CL5,CD5,alphaL,alphaD,adotnorm,umach,ReL,ReD,Re,SectInd,CLdyn5,CDdyn5)      
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
        
        ! Calc tangential added mass increment by analogy to pitching flat plate potential flow theory 
        ! (SAND81-7017) 
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
        
        ! Calc quarter-chord lift and drag coefficient for reference
        CL=CN*cos(alpha)-CT*sin(alpha)
        CD=CN*sin(alpha)+CT*cos(alpha)
                      		
Return
End
