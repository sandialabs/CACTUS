SUBROUTINE AeroCoeffs(alphaL,alphaD,ReL,ReD,Re,adotnorm,umach,SectInd,IsBE,CL,CD,CN,CT,DynamicFlagL,DynamicFlagD) 

	use dystl
	use pidef 
	
        implicit none
        
	integer :: SectInd, IsBE, DynamicFlagL, DynamicFlagD
	real :: alphaL, alphaD, adotnorm, ReL, ReD, Re, umach, CL, CD, CN, CT
        
        integer :: k, BoundInd, isgn
        real :: dalphaRefMax, TransA, dalphaLRef, dalphaDRef, gammal, gammam, dalphaL, dalphaD
        real :: alssn, alssp, alrefL, alLagD, alrefD, delN, delP, CL_L, CL_D, CD_L, CD_D
	
        ! Calculate the static stall enevelope limits                                                                                       
        k=1                                                               
        BoundInd=0
        do while ((k < nstl(SectInd)) .AND. BoundInd==0)
                if (restl(k+1,SectInd) >= Re) then
                        BoundInd=1
                else
                        k=k+1
                end if
        end do                                                               
        alssp=alstlp(k,SectInd)+dapdre(k,SectInd)*(Re-restl(k,SectInd))            
        alssn=alstln(k,SectInd)+dandre(k,SectInd)*(Re-restl(k,SectInd))
          
        ! Limit reference dalpha to a maximum to keep sign of CL the same for
        ! alpha and lagged alpha (considered a reasonable lag...). Note:
        ! magnitude increasing and decreasing effect ratios are maintained.
        dalphaRefMax=min(abs(alssp-alzer(SectInd)),abs(alssn-alzer(SectInd)))/max(k1pos,k1neg)
        TransA=.5*dalphaRefMax ! transition region for fairing lagged AOA in pure lag model
       
        ! Dynamic flags
        DynamicFlagL=0
        DynamicFlagD=0
        
        ! If blade end segment, force into static mode, else test dynamic stall     
        isgn=sign(1.0,adotnorm)                                                        
        if (IsBE == 0) then
        
                ! Modified Boeing-Vertol approach
                
                ! Lift                                                                                                            
                gammal=gammaxl(SectInd)-(umach-smachl(SectInd))*dgammal(SectInd)                 
                dalphaLRef=gammal*sqrt(abs(adotnorm))                              
                dalphaLRef=min(dalphaLRef,dalphaRefMax)
                
                if ((adotnorm*(alphaL-alzer(SectInd))) < 0.0) then     
                        ! Magnitude of CL decreasing                                     
                        dalphaL=k1neg*dalphaLRef                                                     
                        alrefL=alphaL-dalphaL*isgn                                           
                
                        ! switch on lagged alpha
                        if (alrefL <= alssn .OR. alrefL >= alssp) then
                                DynamicFlagL=1     
                        end if
                
                else                
                        ! Magnitude of CL increasing                                     
                        dalphaL=dalphaLRef*k1pos                                               
                        alrefL=alphaL-dalphaL*isgn                                           
                
                        ! switch on alpha
                        if (alphaL <= alssn .OR. alphaL >= alssp) then
                                DynamicFlagL=1    
                        end if
                end if 
                
                ! Drag                                                                                                       
                gammam=gammaxm(SectInd)-(umach-smachm(SectInd))*dgammam(SectInd)               
                if (umach < smachm(SectInd)) then
                        gammam=gammaxm(SectInd)  
                end if
                dalphaDRef=gammam*sqrt(abs(adotnorm))                              
                dalphaDRef=min(dalphaDRef,dalphaRefMax)
                
                if ((adotnorm*(alphaD-alzer(SectInd))) < 0.0) then     
                        ! Magnitude of CL decreasing                                                                                   
                        dalphaD=k1neg*dalphaDRef                                                 
                        alLagD=alphaD-dalphaD*isgn
                        
                        ! switch on lagged alpha
                        delN=alssn-alLagD
                        delP=alLagD-alssp 
                else                
                        ! Magnitude of CL increasing                                                                                    
                        dalphaD=dalphaDRef*k1pos                                           
                        alLagD=alphaD-dalphaD*isgn 
                
                        ! switch on alpha
                        delN=alssn-alphaD
                        delP=alphaD-alssp
                end if
                
                if (delN > TransA .OR. delP > TransA) then
                        alrefD=alLagD
                        DynamicFlagD=1
                elseif (delN >= 0 .AND. delN < TransA) then
                        ! Transition region (fairing effect...)
                        alrefD=alphaD+(alLagD-alphaD)*delN/TransA
                        DynamicFlagD=1 
                elseif (delP >= 0 .AND. delP < TransA) then
                        ! Transition region (fairing effect...)
                        alrefD=alphaD+(alLagD-alphaD)*delP/TransA
                        DynamicFlagD=1 
                end if
        
        end if                                                     
                                                                        
        ! Static or dynamic model       
        if (DynamicFlagL == 1) then
                ! Dynamic stall characteristics                                 
                ! Linear expansion model for linear region coeffs  
                CALL intp(ReL,alrefL*condeg,CL_L,CD_L,SectInd)                                                                                      
                CL_L=CL_L/(alrefL-alzer(SectInd))*(alphaL-alzer(SectInd)) 
        else        
                ! Static characteristics  
                CALL intp(ReL,alphaL*condeg,CL_L,CD_L,SectInd)                                        
        end if
        CL=CL_L
        
        if (DynamicFlagD == 1) then
                ! Dynamic characteristics                                
                ! Pure lag model for drag        
                CALL intp(ReD,alrefD*condeg,CL_D,CD_D,SectInd)                                                                                    
        else    
                ! Static characteristics
                CALL intp(ReD,alphaD*condeg,CL_D,CD_D,SectInd)                                                                              
        end if
        CD=CD_D
        
        ! Tangential and normal coeffs
        CN=CL_L*cos(alphaL)+CD_L*sin(alphaL)                                   
        CT=-CL_D*sin(alphaD)+CD_D*cos(alphaD) 
		
Return
End
