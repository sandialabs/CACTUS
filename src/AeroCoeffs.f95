SUBROUTINE AeroCoeffs(nElem,alpha75,alpha5,Re75,Re5,Re,adotnorm,umach,SectInd,IsBE,CL,CD,CN,CT) 

        use dystl
        use pidef

        implicit none
        
	integer :: SectInd, IsBE, nElem
	real :: alpha75, alpha5, adotnorm, Re75, Re5, Re, umach, CL, CD, CN, CT
        
        real :: CLstat, CLstat5, CDstat, C, CTAM
        
        ! Static characteristics  
        CALL intp(Re75,alpha75*condeg,CLstat,C,SectInd) 
        CALL intp(Re5,alpha5*condeg,CLstat5,C,SectInd)
        CALL intp(Re5,alpha5*condeg,C,CDstat,SectInd) 
        
        ! If blade end segment or no dynamic stall, use static values, else calc dynamic stall
        if (IsBE == 0 .AND. DSFlag/=0) then
                
                if (DSFlag==1) then
                        ! Modified Boeing-Vertol approach
                        Call BV_DynStall(CLstat,CDstat,alpha75,alpha5,adotnorm,umach,Re75,Re5,Re,SectInd,CL,CD)      
                else
                        ! Leishman-Beddoes model
                        Call LB_DynStall(nElem,CLstat,CDstat,alpha75,alpha5,umach,Re,SectInd,CL,CD)
                end if
        
        else
                CL=CLstat
                CD=CDstat
        end if
        
        ! Tangential and normal coeffs
        CN=CL*cos(alpha5)+CD*sin(alpha5)                                   
        CT=-CL*sin(alpha5)+CD*cos(alpha5) 
        
        ! Add in added mass effects (flat plate approx)
        CTAM=-CLstat5/2.0*adotnorm
        CT=CT+CTAM
              		
Return
End
