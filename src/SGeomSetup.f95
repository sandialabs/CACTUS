SUBROUTINE SGeomSetup() 

	use parameters
	
        use configr       
	use element
        use strut 
        use airfoil            
	use pidef      
	
        Implicit None
        
        integer :: i, j, NElem, BIndS, EIndS, BIndE, EIndE, iSc     
        real :: sx, sy, sz, VMag, LR

	! Sets up strut geometry arrays.        
           
	! Set initial geometry (zero theta)                                                            
	do i=1,NStrut  
            NElem=Struts(i)%NElem 
            BIndS=Struts(i)%BIndS
            EIndS=Struts(i)%EIndS
            BIndE=Struts(i)%BIndE
            EIndE=Struts(i)%EIndE 
                                                                                                                                                                     
            ! Get blade thickness to chord
            if (BIndS > 0) then
                iSc=Blades(BIndS)%iSect(EIndS)
                Struts(i)%tcS=tc(iSc) 
            else 
                Struts(i)%tcS=0.0
            end if
            if (BIndE > 0) then
                iSc=Blades(BIndE)%iSect(EIndE)
                Struts(i)%tcE=tc(iSc) 
            else 
                Struts(i)%tcE=0.0
            end if

            ! Init length sum
            LR=0.0
                                                                                                                                                                              
            do j=1,NElem 
                Struts(i)%CRm(j)=0.5*(Struts(i)%CRe(j)+Struts(i)%CRe(j+1))                 
                
                ! Center locations 
                Struts(i)%SCx(j)=0.5*(Struts(i)%SEx(j)+Struts(i)%SEx(j+1))                                      
                Struts(i)%SCy(j)=0.5*(Struts(i)%SEy(j)+Struts(i)%SEy(j+1))                                                                                              
                Struts(i)%SCz(j)=0.5*(Struts(i)%SEz(j)+Struts(i)%SEz(j+1))
                
                ! Spanwise vector 
                sx=Struts(i)%SEx(j+1)-Struts(i)%SEx(j)  
                sy=Struts(i)%SEy(j+1)-Struts(i)%SEy(j)  
                sz=Struts(i)%SEz(j+1)-Struts(i)%SEz(j)
                VMag=sqrt(sx**2+sy**2+sz**2) 
                sx=sx/VMag 
                sy=sy/VMag    
                sz=sz/VMag 
                Struts(i)%sx(j)=sx                                  
                Struts(i)%sy(j)=sy                                                                                         
                Struts(i)%sz(j)=sz
                
                ! Sum length
                LR=LR+VMag

            end do 
            Struts(i)%LR=LR 
              
        end do
        
return                                                            
end 
