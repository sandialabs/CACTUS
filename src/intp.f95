SUBROUTINE intp(RE,ALPHA,CY,CZ,KK)   

	use parameters
	
	use cltab
	use test
				
	real :: CYA(2),CZA(2)                                           
	integer :: U1, X1, iUB, iLB
	logical :: NotDone                                               
	                                                                      
	!    INTERPOLATE ON RE NO. AND ANGLE OF ATTACK TO GET AIRFOIL CHARACTERISTICS                                            
	                                                                      
	CYA(:)=0.0                                                        
	CZA(:)=0.0  
	                                                      
      	if (RE >= TRE(1,KK)) then                                                                                                   
      		
      		! Find Re upper and lower bounds.                                     
		
		NotDone=.true.    
		iUB=2                                                                 
		do while (NotDone)   
			
			if (RE <= TRE(iUB,KK)) then
				! Done
				NotDone=.false.
				if (RE == TRE(iUB,KK)) then
					iLB=iUB
				else
					iLB=iUB-1                                                           
					XRE=(RE-TRE(iLB,KK))/(TRE(iUB,KK)-TRE(iLB,KK))
				end if
			else
				if (iUB == nRET(KK)) then	
					! warning: no upper bound in table, take last point and set warning...
					NotDone=.false.                                                       
					iLB=iUB                                                           
					XRE=0.0                                                           
					IUXTP=1
				else	
					! No upper bound, increment and continue				
					iUB=iUB+1
				end if
			end if

		end do                                                    
	
	else        
		! warning: no lower bound in table, take first point and set warning                                               
		iLB=1                                                             
		iUB=1                                                             
		XRE=0.0                                                                                                 
		ILXTP=1
	end if                                                       
                                                                       
                                                                       
	! INTERPOLATE ON THE ANGLE OF ATTACK                               
                                                                       
      I=1                                                               
      do J=iLB,iUB                                                  
		
		NTB=NTBL(J,KK) ! # of alpha values in table for this section                                                
		
		! Find upper and lower bound indicies on alpha                                                                     
                                                     
		! DO INTERVAL HALVING LOOK UP                                      
																
		U1=NTB                                                                                              
		L1=1                                                              
		X1=NTB/2 
		NotDone=.true. 
		do while (NotDone)                                                        
			if (ALPHA < TA(X1,J,KK)) then
				U1=X1                                                                                                            
			else    
				L1=X1 
			end if    
										
			if ((U1-L1) == 1) then
				NotDone=.false.
			else 
				X1=L1+(U1-L1)/2
			end if                                                                                              
		end do                                         
                                                                      
		! DO STRAIGHT LINE INTERPOLATION ON ALPHA                          
                                                                    
		XA=(ALPHA-TA(L1,J,KK))/(TA(U1,J,KK)-TA(L1,J,KK))                  
		CYA(I)=TCL(L1,J,KK)+XA*(TCL(U1,J,KK)-TCL(L1,J,KK))                
		CZA(I)=TCD(L1,J,KK)+XA*(TCD(U1,J,KK)-TCD(L1,J,KK))     
		           
		I=I+1           
	end do                                                          
                                                                       
	! DO STRAIGHT LINE INTERPOLATION ON RE NO.                         
                                                                       
	CY=CYA(1)+XRE*(CYA(2)-CYA(1))                                     
	CZ=CZA(1)+XRE*(CZA(2)-CZA(1))  
                                   
RETURN                                                            
6002  FORMAT(1H ,5X,10HREY. NO. =,E12.4,27H ABOVE MAX. TABLE VALUE OF ,E12.4,18H. MAX. VALUE USED.)   
6005  FORMAT(1H ,5X,10HREY. NO. =,E12.4,27H BELOW MIN. TABLE VALUE OF ,E12.4,21H. MINIMUM VALUE USED.)                                
END                                                               
