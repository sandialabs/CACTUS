SUBROUTINE FIVEL(X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,G,UU,VV,WW,NTS,NTT,N10)      

      use vortex
                                                                             

! For some reason, the original code had dimensions of 1 for the arrays
! below. 300 has been arbitrarily defined and seems to work.
      DIMENSION X1(300),X2(300),Y1(300),Y2(300),Z1(300),Z2(300),G(300)
!                                                                       
!    CALCULATE THE VELOCITY INDUCED AT A POINT DUE TO A VORTEX AT       
!    SOME OTHER POINT                                                   
!                                                                       
!                  DO LOOP, L FROM 1 TO NTT                             
!                  CALL TO FIVEL---                                     
!                       X1(1,K) TO USE X1(L,K)                          
                                       

      DO 20 L=NTS,NTT                                                   
!                              X2T = X2(L+N10)                          
      AX = X2(L+N10)-X1(L)                                              
      BX = X2(L+N10)-X3                                                 
!                              Y2T = Y2(L+N10)                          
      AY = Y2(L+N10)-Y1(L)                                              
      BY = Y2(L+N10)-Y3                                                 
!                              Z2T = Z2(L+N10)                          
      AZ = Z2(L+N10)-Z1(L)                                              
      BZ = Z2(L+N10)-Z3                                                 
      B = SQRT(BX**2 + BY**2 + BZ**2)                                   

! Cray vector instruction.
!     B = CVMGT(B,1.0,B.NE.0.0)

! Non-vectorized version.
      IF ( B .EQ. 0.0 ) B = 1.0

      ADBDB = (AX*BX + AY*BY + AZ*BZ)/B                                 
!                                                                       
      CX = X1(L)-X3                                                     
      CY = Y1(L)-Y3                                                     
      CZ = Z1(L)-Z3                                                     
!                                                                       
      C = SQRT(CX**2 + CY**2 + CZ**2)                                   

! Cray vector instruction.
!     C = CVMGT(C,1.0,C.NE.0.0)

! Non-vectorized version.
      IF ( C .EQ. 0.0 ) C = 1.0

      ADCDC = (AX*CX + AY*CY + AZ*CZ)/C                                 
!                                                                       
      CCAX = CY*AZ - AY*CZ                                              
      CCAY = AX*CZ - CX*AZ                                              
      CCAZ = CX*AY - AX*CY                                              
      A2 =AX**2 + AY**2 + AZ**2                                         
!                                                                       
      CCAV = CCAX**2 + CCAY**2 + CCAZ**2                                

! Cray vector instruction.
!     CCAV = CVMGT(A2*VRAD2, CCAV, CCAV.GE.1.0E-07 .AND.
!                                  CCAV.LT.A2*VRAD2)

! Non-vectorized version.
      IF ((CCAV .GE. 1.0E-07).AND.(CCAV .LT. A2*VRAD2)) CCAV = A2*VRAD2

! Cray vector instruction.
!     CCAV = CVMGT(CCAV, 1.0, CCAV.GE.1.0E-07)
!     VF   = CVMGT((ADBDB-ADCDC)*G(L)/(12.56637*CCAV), 0.0,
!                                                       CCAV.GE.1.0E-07)

! Non-vectorized version.
      IF ( CCAV .GE. 1.0E-07 ) THEN
	VF   = (ADBDB-ADCDC)*G(L)/(12.56637*CCAV)
      ELSE
	CCAV = 1.0
	VF   = 0.0
      ENDIF

!                  MOVE  20 CONTINUE TO HERE                            

      UU = UU + CCAX*VF                                                 
      VV = VV + CCAY*VF                                                 
      WW = WW + CCAZ*VF                                                 
   20 CONTINUE                                                          
RETURN                                                            
END                                                               
