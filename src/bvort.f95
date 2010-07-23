SUBROUTINE bvort(nGeom,NLTol,iConv)  

	use gam
	use test
	use configr
        use regtest       
	
        integer IsBE, DynamicFlag
        real alpha, Re, umach, ur, CN, CT, te, NLTol, dgb
	
	! Calculates bound and new shed vorticity
	iConv=0											                                           
	do i=1,nb                                                                                                     
		nei=1+(i-1)*(nbe+1)                                               
		do j=1,nbe                                                     
			nej=nei+j                                                         
			nej1=nej-1 
			
			IsBE=0
			if (j==1 .OR. j==nbe) then
				IsBE=1
			end if
			                                                                                                                             
			! Calculate the loads on the blade segment                                                                                              
			CALL bsload(nej,nGeom,IsBE,DynamicFlag,alpha,Re,umach,ur,CN,CT,te) 
			                                                                                  
			! Calculate the bound vortex strength change                                                                                                   
			dgb=abs((GB(nej1)-GS(nt,nej1))/GB(nej1))                          
			
			! If change outside tolerance for any element, set flag
			if (dgb .gt. NLTol) iConv=1 
			 
			! Regression test
                        if (RegTFlag == 1) then
                                Reg_ElemNum=nej1 
                                Reg_DF=DynamicFlag
                                Reg_ElemAOA=alpha*180.0/3.14159 
                                Reg_ElemCirc=GB(nej1)
                                Reg_dElemCirc=dgb   
                                Call WriteRegTOutput(1)                  
			end if 
			
			! Set the bound circulation as the current entry in the spanwise
			! vorticity array for velocity calculation (and eventual wake convection)                                  
			GS(nt,nej1)=GB(nej1)                                              
		end do 
	end do                                                         
	
Return
End                                                               
