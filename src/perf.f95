SUBROUTINE perf(nGeom,cpl) 

	use dystl
	use wakeloc
	use element
	use vel
	use gam
	use shear
	use pidef
	use test
	use rvrr
	use configr
					
	integer IsBE, DynamicFlag
	real alpha, Re, umach, ur, CN, CT, te                         
	real BladeTorque(3) ! Total torque on each blade calculated for output...                                              
									
	! Calculate the local lift on each blade segment and therefore get the bound vortex strengths 
	! and the turbine performance for this timestep					   
																	
	tr=0.0                                                            
	cpl=0.0                                                                                                                        
	write(6,601) nt,irev,Theta(nGeom)                                         
	write(9,901) nt,irev,Theta(nGeom)                                                                           
	do i=1,nb                                                      
		BladeTorque(i)=0.0                                                       
										
		! Retrieve the blade segment geometric information                 							                                              
		nei=1+(i-1)*(nbe+1)  
		
		! Calculate the loads on the blade segment                                              
		do j=1,nbe                                                     
											
			nej=nei+j                                                         
			nej1=nej-1 
			
			IsBE=0
			if (j == 1 .OR. j == nbe) then
				IsBE=1
			end if
			
			CALL bsload(nej,nGeom,IsBE,DynamicFlag,alpha,Re,umach,ur,CN,CT,te)                                                                     
											
			! Calculate the bound vortex strength and the turbine performance                                                      
			alphaDeg=alpha*condeg                                                                                                           
			write(9,904) DynamicFlag,alphaDeg,Re,umach,CN,CT,ur,te						
			alfold(nej1)=alpha                                                
			GS(nt,nej1)=GB(nej1)                                              
			BladeTorque(i)=BladeTorque(i)+te                                                
			tr=tr+te                                                          
			cpl=cpl+te*ut ! Increment to the power coeff. integration for this timestep                                                    
		end do                                                                                
	end do                                                          
	write(9,903) BladeTorque(1),BladeTorque(2),BladeTorque(3),tr,cpl                              
	write(6,609) tr,cpl                                              

Return                                                            
601 FORMAT('0',20X,'NT=',I5,'   IREV=',I5,10X,'ANGLE=',F9.2,' DEG.')                   
609 FORMAT(' ',10X,'ROTOR TORQUE COEFFICIENT=',E15.6,/,' ',10X,'ROTOR POWER COEFFICIENT=',E15.6)                                 
901 FORMAT(2I5,F10.2)                                                 
903 FORMAT(10E12.5)
904 FORMAT(I3,7E12.5)                                                    
End                                                               
