SUBROUTINE bsload(nElem,nGeom,IsBE,alpha,Re,umach,ur,CN,CT,Fx,Fy,Fz,te) 
       
        use element
        use blade
        use pidef       
        use configr
        use airfoil
        use dystl
        use cltab
        use freestream
        use test  
                                                      
        integer SectInd, nElem, nGeom, IsBE
        real ElemSpanR, ElemChordR, xe, ye, ze, nxe, nye, nze, txe, tye, tze, alpha, alpdot, adotnorm, te
        real uAve, vAve, wAve, uBlade, vBlade, wBlade
        
        
        ! Calculates aero loads on a blade element. Static and dynamic airfoil characteristics calculated here...
                                                
        nElem1=nElem-1  ! Blade element is referenced by its upper end location index nElem, with nElem-1 being the lower end location index.                                             
         
        
        ! Retrieve the blade segment geometric information                  
                                                                        
        ! Element span                                                      
        ElemSpanR=eSpan(nElem)  
        ElemChordR=eChord(nElem)                                                    
        
        ! Quarter chord location
        xe=0.5*(xBE(nGeom,nElem)+xBE(nGeom,nElem1))
        ye=0.5*(yBE(nGeom,nElem)+yBE(nGeom,nElem1))
        ze=0.5*(zBE(nGeom,nElem)+zBE(nGeom,nElem1))
        
        ! Element normal and tangential vectors
        nxe=nx(nGeom,nElem)                                                  
        nye=ny(nGeom,nElem)                                                  
        nze=nz(nGeom,nElem)
        txe=tx(nGeom,nElem)                                                  
        tye=ty(nGeom,nElem)                                                  
        tze=tz(nGeom,nElem)
        
        ! Airfoil section                                                 
        SectInd=isect(nElem)                                                     
                        
                                                        
        ! Calculate the local blade segment angle of attack                 

        ! Wall and wake induced velocity                                                   
        uAve=(u(nt,nElem)+u(nt,nElem1))/2.0                                   
        vAve=(v(nt,nElem)+v(nt,nElem1))/2.0                                   
        wAve=(w(nt,nElem)+w(nt,nElem1))/2.0    
        ! Freestream velocity
        uFSAve=(uFS(nt,nElem)+uFS(nt,nElem1))/2.0                                   
        vFSAve=(vFS(nt,nElem)+vFS(nt,nElem1))/2.0                                   
        wFSAve=(wFS(nt,nElem)+wFS(nt,nElem1))/2.0 
        ! Blade velocity due to rotation                                                      
        CALL CalcBladeVel(wRotX,wRotY,wRotZ,xe,ye,ze,uBlade,vBlade,wBlade)
        
        ! Calc element normal and tangential velocity components
        urdn = (nxe*(uAve+uFSAve-uBlade)+nye*(vAve+vFSAve-vBlade)+nze*(wAve+wFSAve-wBlade))     ! Normal
        urdc = (txe*(uAve+uFSAve-uBlade)+tye*(vAve+vFSAve-vBlade)+tze*(wAve+wFSAve-wBlade))     ! Tangential

        ur=sqrt(urdn**2+urdc**2)                                          
        alpha=atan2(urdn,urdc) 
        dal=alpha-AOA_Last(nelem1)                                                                                       
        adotnorm=dal/DT*ElemChordR/(2.0*ur)  ! adot*c/(2*U)
                                     
        Re=ReM*ElemChordR*ur                                                         
        umach=ur*Minf                                                 
        
        !---------
        ! JCM: These .5c and .75c locations are used to approx. pitch rate and alphadot effects in the 
        ! static coeff. calculation. Could be exchanged for tables...
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
        Re5=ReM*ElemChordR*ur5                                                       
        alpha5=atan2(urdn5,urdc) 
        urdn75 = (nxe*(uAve+uFSAve-uBlade75)+nye*(vAve+vFSAve-vBlade75)+nze*(wAve+wFSAve-wBlade75))     
        ur75=sqrt(urdn75**2+urdc**2)                                        
        Re75=ReM*ElemChordR*ur75                                                       
        alpha75=atan2(urdn75,urdc)
        !--------
        
        ! Evaluate aero coefficients and dynamic stall effects as appropriate
        Call AeroCoeffs(nElem,alpha75,alpha5,Re75,Re5,Re,adotnorm,umach,SectInd,IsBE,CL,CD,CN,CT)
                
        ! Bound vortex strength from CL via Kutta-Joukowski analogy. 
        ! Save corresponding AOA as well                                                                                          
        GB(nElem1)=CL*ElemChordR*ur/2.0  
        AOA(nElem1)=alpha
        ! normalized time step used to update states in the LB model
        ds(nElem)=2.0*ur*DT/ElemChordR
        
        ! Force coeff. from this blade element, re-referenced to full turbine scale (F/(1/2*rho*Uinf^2*At))                                         
        FN=CN*(ElemChordR*ElemSpanR/at)*ur**2                                                       
        FT=CT*(ElemChordR*ElemSpanR/at)*ur**2         
        ! Corresponding torque coeff. (T/(1/2*rho*Uinf^2*At*R))                                              
        Fx=FN*nxe+FT*txe
        Fy=FN*nye+FT*tye
        Fz=FN*nze+FT*tze
        CALL cross(xe,ye,ze,Fx,Fy,Fz,TRx,TRy,TRz)
        te=TRx*RotX+TRy*RotY+TRz*RotZ                             
      
Return                                                                                                                                                           
End                                                               
