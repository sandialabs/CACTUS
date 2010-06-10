SUBROUTINE INTERP3(XP,YP,ZP,UIN,VIN,WIN) 

      ! JCM: Note that this is only used in the fixed wake velocity calculation in swivel...

      use xwake
      use uwake


      DIMENSION UX(2,2),VX(2,2),WX(2,2),UQ(2,2,2),VQ(2,2,2),WQ(2,2,2)   
      DIMENSION UY(2),VY(2),WY(2)                                       
!                                                                       
!      DO 3-DIMENSIONAL LINEAR INTERPOLATION TO DETERMINE THE VELOCITY  
!    AT POINT XP,YP,ZP                                                  
!                                                                       
      IW=1                                                              
      JW=1                                                              
      KW=1                                                              
!                                                                       
!      BRACKET THE POINT WITHIN THE KNOWN ARRAYS                        
!                                                                       
5     CONTINUE                                                          
      IF (XP .LE. XFW(IW)) GO TO 10                                     
      IW=IW+1                                                           
      GO TO 5                                                           
10    CONTINUE                                                          
      IF (YP .LE. YFW(JW)) GO TO 15                                     
      JW=JW+1                                                           
      GO TO 10                                                          
15    CONTINUE                                                          
      IF (ZP .LE. ZFW(KW)) GO TO 20                                     
      KW=KW+1                                                           
      GO TO 15                                                          
20    CONTINUE                                                          
      IW=IW-1                                                           
      JW=JW-1                                                           
      KW=KW-1                                                           
      HX=XFW(IW+1)-XFW(IW)                                              
      HY=YFW(JW+1)-YFW(JW)                                              
      HZ=ZFW(KW+1)-ZFW(KW)                                              
!                                                                       
!      LOAD THE INTERPOLATION ARRAYS WITH THE VALUES OF THE 8 KNOWN     
!      POINTS IN THE CUBE SURROUNDING THE DESIRED POINT                 
!                                                                       
      DO 30 IQ=1,2                                                      
      DO 30 JQ=1,2                                                      
      DO 30 KQ=1,2                                                      
      IS=IW+IQ-1                                                        
      JS=JW+JQ-1                                                        
      KS=KW+KQ-1                                                        
      UQ(IQ,JQ,KQ)=UFW(IS,JS,KS)                                        
      VQ(IQ,JQ,KQ)=VFW(IS,JS,KS)                                        
      WQ(IQ,JQ,KQ)=WFW(IS,JS,KS)                                        
30    CONTINUE                                                          
      DHX=(XP-XFW(IW))/HX                                               
      DHY=(YP-YFW(JW))/HY                                               
      DHZ=(ZP-ZFW(KW))/HZ                                               
      DO 50 J=1,2                                                       
!                                                                       
!      INTERPOLATE TO CORRECT X LOCATION                                
!                                                                       
      DO 40 I=1,2                                                       
      UX(I,J)=UQ(1,I,J)*(1.0-DHX)+UQ(2,I,J)*DHX                         
      VX(I,J)=VQ(1,I,J)*(1.0-DHX)+VQ(2,I,J)*DHX                         
      WX(I,J)=WQ(1,I,J)*(1.0-DHX)+WQ(2,I,J)*DHX                         
40    CONTINUE                                                          
!                                                                       
!      INTERPOLATE TO CORRECT Y LOCATION                                
!                                                                       
      UY(J)=UX(1,J)*(1.0-DHY)+UX(2,J)*DHY                               
      VY(J)=VX(1,J)*(1.0-DHY)+VX(2,J)*DHY                               
      WY(J)=WX(1,J)*(1.0-DHY)+WX(2,J)*DHY                               
50    CONTINUE                                                          
!                                                                       
!      FINISH UP BY INTERPOLATING TO CORRECT Z LOCATION                 
!                                                                       
      UIN=UY(1)*(1.0-DHZ)+UY(2)*DHZ                                     
      VIN=VY(1)*(1.0-DHZ)+VY(2)*DHZ                                     
      WIN=WY(1)*(1.0-DHZ)+WY(2)*DHZ                                     
RETURN                                                            
END                                                               
