SUBROUTINE SGeomSetup() 

	use parameters
	
        use configr       
	use element
        use strut 
        use airfoil            
	use pidef      
	
        Implicit None
        
        integer i, j, NElem, BInd, EInd, iSc      

	! Sets up strut geometry arrays.        
           
	! Set initial geometry (zero theta)                                                            
	do i=1,NStrut  
            NElem=Struts(i)%NElem 
            BInd=Struts(i)%BInd
            EInd=Struts(i)%EInd  
                                                                                                                                                                     
            ! Get blade thickness to chord
            iSc=Blades(BInd)%iSect(EInd)
            Struts(i)%tc=tc(iSc) 
                                                                                                                                                                              
            do j=1,NElem 
                Struts(i)%CRm(j)=0.5*(Struts(i)%CRe(j)+Struts(i)%CRe(j+1))                 
                
                ! Center locations 
                Struts(i)%SCx(j)=0.5*(Struts(i)%SEx(j)+Struts(i)%SEx(j+1))                                      
                Struts(i)%SCy(j)=0.5*(Struts(i)%SEy(j)+Struts(i)%SEy(j+1))                                                                                              
                Struts(i)%SCz(j)=0.5*(Struts(i)%SEz(j)+Struts(i)%SEz(j+1))
            end do   
        end do
        
return                                                            
end 
