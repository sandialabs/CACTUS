SUBROUTINE DynStall(alpha,adotnorm,sgn,Re,umach,SectInd,IsBE,IsZX,DynamicFlag,CL,CD,CN,CT) 

	use dystl
	use pidef 
	use test
	
	integer k, BoundInd, SectInd, IsBE, IsZX, DynamicFlag
	real alpha, adotnorm, Re, umach, CL, CD, CN, CT, c1, sgn
	
	! Evaluates dynamic stall model, if active.
	
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
									
	! If crossing zero alpha or blade end segment, force into static mode, else test dynamic stall                                                            
	if (IsZX == 0 .AND. IsBE == 0) then 
							
		! See if alpha is in dynamic stall region                           
										                                                                                   
		gammal=gammaxl(SectInd)-(umach-smachl(SectInd))*dgammal(SectInd)                 
		gammam=gammaxm(SectInd)-(umach-smachm(SectInd))*dgammam(SectInd)                 
		if (umach < smachm(SectInd)) gammam=gammaxm(SectInd)                        
		dalpha=gammal*sqrt(adotnorm)                              
		dalpham=dalpha*gammam/gammal                                      
		
		if ((dal*alpha) < 0.0) then	
									
			! Magnitude of alpha decreasing                                     
										
			dalpha=k1neg*dalpha                                               
			dalpham=k1neg*dalpham                                             
			alssp1=alssp+dalpha                                               
			alssp2=alssp-dalpha                                               
			alssn1=alssn+dalpha                                               
			alssn2=alssn-dalpha       
			alref=alpha-dalpha*sgn                                           
			alrefm=alpha-dalpham*sgn
							
			if (alpha <= alssn1 .OR. alpha >= alssp2) then
				DynamicFlag=1	
			end if
			
		else		
						
			! Magnitude of alpha increasing                                     
										
			dalpha=dalpha*k1pos                                               
			dalpham=dalpham*k1pos                                             
			alssp1=alssp+dalpha                                               
			alssn2=alssn-dalpha
			alref=alpha-dalpha*sgn                                           
			alrefm=alpha-dalpham*sgn 
									
			if (alpha <= alssn .OR. alpha >= alssp) then
				
				DynamicFlag=1     
				                                         						
				! If alref and alpha are different signs, set alref equal to the static stall limit                                			
				if ((alref*alpha) <= 0.0) then                            
					                                                      
					if (alpha < 0.0) then
						alref=alssn  
					else
						alref=alssp 
					end if
					
					! Check for bad alrefm					
					if ((alrefm*alpha) <= 0.0) then                                         
						alrefm=alpha                                                      
						ierr0=1 
					end if  
					
				end if  
				                                                        
			end if                                                        
		end if ! Magnitude alpha decreasing                                                         									
	end if   ! Potentially dynamic stall                                                      
								
	
	if (DynamicFlag == 1) then 
										
		! Get dynamic stall characteristics                                 
										
		alpr=alref*condeg                                                 
		CALL intp(Re,alpr,CL,CD,SectInd)                                       
		CL=CL/(alref-alzer(SectInd))*(alpha-alzer(SectInd))                         
		alprm=alrefm*condeg                                               
		CD1=CD                                                            
		CALL intp(Re,alprm,CL1,CD,SectInd)                                     
		CN=CL*cos(alpha)+CD*sin(alpha)                                   
		CT=-CL*sin(alpha)+CD*cos(alpha)                                    
	else
		CL=0.0
		CD=0.0
		CN=0.0
		CT=0.0		
	end if
		
Return
End
