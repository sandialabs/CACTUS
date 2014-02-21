SUBROUTINE bsload(nElem,IsBE,alpha,alpha5,alpha75,adotnorm,Re,umach,ur,CL,CD,CM25,CLCirc,CN,CT,Fx,Fy,Fz,te) 

    use element
    use blade
    use pidef       
    use configr
    use airfoil
    use dystl
    use util

    implicit none

    integer nElem, IsBE
    real alpha, alpha5, alpha75, adotnorm, Re, umach, ur, CL, CD, CM25, CLCirc, CN, CT, Fx, Fy, Fz, te

    integer SectInd, nElem1
    real ElemAreaR, ElemChordR, xe, ye, ze, nxe, nye, nze, txe, tye, tze, sxe, sye, sze
    real dal, wP, wPNorm
    real uAve, vAve, wAve, uFSAve, vFSAve, wFSAve, uBlade, vBlade, wBlade, urdn, urdc
    real xe5, ye5, ze5, xe75, ye75, ze75, uBlade5, vBlade5, wBlade5, uBlade75, vBlade75, wBlade75
    real urdn5, ur5, urdn75, ur75   
    real FN, FT, MS, TRx, TRy, TRz, CircDir  


    ! Calculates aero loads on a blade element. Static and dynamic airfoil characteristics calculated here...

    nElem1=nElem-1  ! Blade element is referenced by its upper end location index nElem, with nElem-1 being the lower end location index.                                             


    ! Retrieve the blade segment geometric information                  

    ! Element span                                                      
    ElemAreaR=eArea(nElem)  
    ElemChordR=eChord(nElem)                                                    

    ! Quarter chord location
    xe=xBC(nElem)
    ye=yBC(nElem)
    ze=zBC(nElem)

    ! Element normal, tangential, and spanwise vectors
    nxe=nxBC(nElem)                                                  
    nye=nyBC(nElem)                                                  
    nze=nzBC(nElem)
    txe=txBC(nElem)                                                  
    tye=tyBC(nElem)                                                  
    tze=tzBC(nElem)
    sxe=sxBC(nElem)                                                  
    sye=syBC(nElem)                                                  
    sze=szBC(nElem)

    ! Direction of circulation in wake grid at positive lift
    CircDir=CircSign(nElem)

    ! Airfoil section                                                 
    SectInd=isect(nElem)                                                     


    ! Calculate the local blade segment angle of attack                 

    ! Wall and wake induced velocity                                                   
    uAve=(uB(nElem)+uB(nElem1))/2.0
    vAve=(vB(nElem)+vB(nElem1))/2.0
    wAve=(wB(nElem)+wB(nElem1))/2.0
    ! Freestream velocity
    uFSAve=(uFSB(nElem)+uFSB(nElem1))/2.0
    vFSAve=(vFSB(nElem)+vFSB(nElem1))/2.0
    wFSAve=(wFSB(nElem)+wFSB(nElem1))/2.0
    ! Blade velocity due to rotation                                                      
    CALL CalcBladeVel(wRotX,wRotY,wRotZ,xe,ye,ze,uBlade,vBlade,wBlade)

    ! Calc element normal and tangential velocity components. Calc element pitch rate.
    urdn = (nxe*(uAve+uFSAve-uBlade)+nye*(vAve+vFSAve-vBlade)+nze*(wAve+wFSAve-wBlade))     ! Normal
    urdc = (txe*(uAve+uFSAve-uBlade)+tye*(vAve+vFSAve-vBlade)+tze*(wAve+wFSAve-wBlade))     ! Tangential
    wP = sxe*wRotX+sye*wRotY+sze*wRotZ

    ur=sqrt(urdn**2+urdc**2)                                          
    alpha=atan2(urdn,urdc) 
    wPNorm=wP*ElemChordR/(2.0*max(ur,0.001))        ! wP*c/(2*U)

    Re=ReM*ElemChordR*ur                                                         
    umach=ur*Minf                                                 

    !---------
    ! These .5c and .75c locations are used to calc pitch rate effects
    if (PRFlag/=0) then
        xe5=xe+0.25*ElemChordR*txe
        ye5=ye+0.25*ElemChordR*tye
        ze5=ze+0.25*ElemChordR*tze
        xe75=xe+0.5*ElemChordR*txe
        ye75=ye+0.5*ElemChordR*tye
        ze75=ze+0.5*ElemChordR*tze
        CALL CalcBladeVel(wRotX,wRotY,wRotZ,xe5,ye5,ze5,uBlade5,vBlade5,wBlade5)
        CALL CalcBladeVel(wRotX,wRotY,wRotZ,xe75,ye75,ze75,uBlade75,vBlade75,wBlade75)
        urdn5 = (nxe*(uAve+uFSAve-uBlade5)+nye*(vAve+vFSAve-vBlade5)+nze*(wAve+wFSAve-wBlade5)) 
        ur5=sqrt(urdn5**2+urdc**2)                                                                                             
        alpha5=atan2(urdn5,urdc) 
        urdn75 = (nxe*(uAve+uFSAve-uBlade75)+nye*(vAve+vFSAve-vBlade75)+nze*(wAve+wFSAve-wBlade75))     
        ur75=sqrt(urdn75**2+urdc**2)                                                                                             
        alpha75=atan2(urdn75,urdc)
    else
        alpha5=alpha
        alpha75=alpha
        wPNorm=0.0
    end if

    dal=alpha75-AOA_Last(nelem1)                                                                                       
    adotnorm=dal/DT*ElemChordR/(2.0*max(ur,0.001))  ! adot*c/(2*U)
    !--------

    ! Evaluate aero coefficients and dynamic stall effects as appropriate
    Call AeroCoeffs(nElem,alpha75,alpha5,Re,wPNorm,adotnorm,umach,SectInd,IsBE,CL,CD,CN,CT,CLCirc,CM25)

    ! Bound vortex strength from CL via Kutta-Joukowski analogy. 
    ! Save corresponding AOA as well                                                                                          
    GB_Raw(nElem1)=CircDir*(CLCirc*ElemChordR*ur/2.0)
    !         AOA(nElem1)=alpha
    AOA(nElem1)=alpha75
    ! normalized time step used to update states in the LB model
    ds(nElem)=2.0*ur*DT/ElemChordR

    ! Force and moment coeff. from this blade element, re-referenced to full turbine scale 
    ! (F/(1/2*rho*Uinf^2*At) and M/(1/2*rho*Uinf^2*At*R)                                         
    FN=CN*(ElemAreaR/at)*ur**2                                                       
    FT=CT*(ElemAreaR/at)*ur**2   
    MS=CM25*ElemChordR*(ElemAreaR/at)*ur**2      
    ! Corresponding torque coeff. (T/(1/2*rho*Uinf^2*At*R))                                              
    Fx=FN*nxe+FT*txe
    Fy=FN*nye+FT*tye
    Fz=FN*nze+FT*tze
    CALL cross(xe-RotPX,ye-RotPY,ze-RotPZ,Fx,Fy,Fz,TRx,TRy,TRz)
    te=(TRx*RotX+TRy*RotY+TRz*RotZ)+MS*(sxe*RotX+sye*RotY+sze*RotZ)                           

    Return                                                                                                                                                           
End SUBROUTINE bsload
