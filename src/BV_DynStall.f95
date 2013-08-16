SUBROUTINE BV_DynStall(nElem,CLstat,CDstat,alpha,adotnorm,umach,Re,SectInd,CL,CD)

    use airfoil
    use dystl
    use pidef 

    implicit none

    real :: CLstat, CDstat, CL, CD, alpha, adotnorm, umach, Re
    integer :: SectInd, nElem

    integer :: isgn
    real :: dalphaRefMax, TransA, dalphaLRef, dalphaDRef, dalphaL, dalphaD, Fac
    real :: diff, smachl, hmachl, gammaxl, dgammal, smachm, hmachm, gammaxm, dgammam, gammal, gammam
    real :: alssn, alssp, alrefL, alLagD, alrefD, delN, delP, C, C1, AOA0

    ! Calculate the static stall envelope limits and other params 
    CALL CalcBVStallAOALim(Re,SectInd,alssp,alssn)                                                                                     
    AOA0=alzer(SectInd)
    diff=0.06-tc(SectInd)                                                   
    smachl=0.4+5.0*diff                                            
    hmachl=0.9+2.5*diff                                            
    gammaxl=1.4-6.0*diff                                           
    dgammal=gammaxl/(hmachl-smachl)                       
    smachm=0.2                                                     
    hmachm=0.7+2.5*diff                                            
    gammaxm=1.0-2.5*diff                                           
    dgammam=gammaxm/(hmachm-smachm)          

    ! Limit reference dalpha to a maximum to keep sign of CL the same for
    ! alpha and lagged alpha (considered a reasonable lag...). Note:
    ! magnitude increasing and decreasing effect ratios are maintained.
    Fac=.9 ! Margin to ensure that dalphaRef is never large enough to make alrefL == AOA0 (blows up linear expansion model)
    dalphaRefMax=Fac*min(abs(alssp-AOA0),abs(alssn-AOA0))/max(k1pos,k1neg)
    TransA=.5*dalphaRefMax ! transition region for fairing lagged AOA in pure lag model

    isgn=sign(1.0,adotnorm)                                                        

    ! Modified Boeing-Vertol approach

    ! Lift                                                                                                            
    gammal=gammaxl-(umach-smachl)*dgammal               
    dalphaLRef=gammal*sqrt(abs(adotnorm))                              
    dalphaLRef=min(dalphaLRef,dalphaRefMax)

    if ((adotnorm*(alpha-AOA0)) < 0.0) then     
        ! Magnitude of CL decreasing                                     
        dalphaL=k1neg*dalphaLRef                                                     
        alrefL=alpha-dalphaL*isgn                                           

        ! Only switch DS off using lagged alpha
        if (BV_DynamicFlagL(nElem) == 1 .AND. (alrefL > alssn .AND. alrefL < alssp)) then
            BV_DynamicFlagL(nElem)=0     
        end if

    else                
        ! Magnitude of CL increasing                                     
        dalphaL=dalphaLRef*k1pos                                               
        alrefL=alpha-dalphaL*isgn                                           

        ! switch DS on or off using alpha
        if (alpha <= alssn .OR. alpha >= alssp) then
            BV_DynamicFlagL(nElem)=1
        else     
            BV_DynamicFlagL(nElem)=0
        end if
    end if

    ! Drag                                                                                                       
    gammam=gammaxm-(umach-smachm)*dgammam              
    if (umach < smachm) then
        gammam=gammaxm  
    end if
    dalphaDRef=gammam*sqrt(abs(adotnorm))                              
    dalphaDRef=min(dalphaDRef,dalphaRefMax)

    if ((adotnorm*(alpha-AOA0)) < 0.0) then     
        ! Magnitude of CL decreasing                                                                                   
        dalphaD=k1neg*dalphaDRef                                                 
        alLagD=alpha-dalphaD*isgn

        ! Only switch DS off using lagged alpha
        if (BV_DynamicFlagD(nElem) == 1) then
            delN=alssn-alLagD
            delP=alLagD-alssp
        else
            delN=0.0
            delP=0.0
        end if
    else                
        ! Magnitude of CL increasing                                                                                    
        dalphaD=dalphaDRef*k1pos                                           
        alLagD=alpha-dalphaD*isgn 

        ! switch DS on or off using alpha
        delN=alssn-alpha
        delP=alpha-alssp
    end if

    if (delN > TransA .OR. delP > TransA) then
        alrefD=alLagD
        BV_DynamicFlagD(nElem)=1
    elseif (delN > 0 .AND. delN < TransA) then
        ! Transition region (fairing effect...)
        alrefD=alpha+(alLagD-alpha)*delN/TransA
        BV_DynamicFlagD(nElem)=1 
    elseif (delP > 0 .AND. delP < TransA) then
        ! Transition region (fairing effect...)
        alrefD=alpha+(alLagD-alpha)*delP/TransA
        BV_DynamicFlagD(nElem)=1 
    else
        BV_DynamicFlagD(nElem)=0
    end if

    ! Static or dynamic model       
    if (BV_DynamicFlagL(nElem) == 1) then
        ! Dynamic stall characteristics                                 
        ! Linear expansion model for linear region coeffs  
        CALL intp(Re,alrefL*condeg,CL,C,C1,SectInd)                                                                                      
        CL=CL/(alrefL-AOA0)*(alpha-AOA0) 
    else        
        ! Static characteristics  
        CL=CLstat                                      
    end if

    if (BV_DynamicFlagD(nElem) == 1) then
        ! Dynamic characteristics                                
        ! Pure lag model for drag        
        CALL intp(Re,alrefD*condeg,C,CD,C1,SectInd)                                                                                    
    else    
        ! Static characteristics    
        CD=CDstat                                                                              
    end if


    ! Diagnostic output
    BV_alpha=alpha*condeg
    BV_adotnorm=adotnorm
    BV_alrefL=alrefL*condeg           
    BV_alrefD=alrefD*condeg              

    Return
End SUBROUTINE BV_DynStall
