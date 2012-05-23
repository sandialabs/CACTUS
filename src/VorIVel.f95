subroutine VorIVel(VFlag,G,X1,Y1,Z1,X2,Y2,Z2,XP,YP,ZP,UP,VP,WP)      

	use vortex
			
        integer :: VFlag       ! 0 for bound vortex element, 1 for trailing wake element, 2 for spanwise wake element (sets core radius)
                        
	AX = X2-X1 
	AY = Y2-Y1 
	AZ = Z2-Z1
	A2 = AX**2 + AY**2 + AZ**2
							
	BX = X2-XP                                                                                             
	BY = Y2-YP                                                 
	BZ = Z2-ZP                                        
	B = SQRT(BX**2 + BY**2 + BZ**2)                                   
	
	CX = X1-XP                                                     
	CY = Y1-YP                                                     
	CZ = Z1-ZP                                                                                                                        
	C = SQRT(CX**2 + CY**2 + CZ**2) 
	
	if (B == 0.0) then 
		B = 1.0
	end if 
	if (C == 0.0) then
		C = 1.0
	end if
	
	ADBDB = (AX*BX + AY*BY + AZ*BZ)/B
	ADCDC = (AX*CX + AY*CY + AZ*CZ)/C                                 
									
	CCAX = CY*AZ - AY*CZ                                              
	CCAY = AX*CZ - CX*AZ                                              
	CCAZ = CX*AY - AX*CY                                              
														
	CCAV = CCAX**2 + CCAY**2 + CCAZ**2                                
        
        ! Select the vortex core radius for bound vortex, trailing
        ! or spanwise wake elements
        if (VFlag == 1) then
            VRAD2=VRAD2_T
        else if (VFlag == 2) then
            VRAD2=VRAD2_S
        else
            VRAD2=VRAD2_B
        end if

	if (CCAV >= 1.0E-07) then
		if (CCAV < A2*VRAD2) then
			CCAV = A2*VRAD2
		end if
		VF   = (ADBDB-ADCDC)*G/(12.56637*CCAV)
	else
		VF   = 0.0
	end if
				
	UP = UP+CCAX*VF                                                 
	VP = VP+CCAY*VF                                                 
	WP = WP+CCAZ*VF
                                                                                                               
RETURN                                                            
END                                                               
