SUBROUTINE SWIVEL(NT,ntTerm,NBE,NB,NE,IUT,NSW,NPW)    

      use parameters
      
      use xwake
      use uwake
      use wakeloc
      use vel
      use veo
      use gam
      use test
      use freestream
   
      	integer ygcErr
   
   
!    JCM: Update fixed wake grid velocities at each point in the grid (if on an update iteration)
!    and interpolate these velocities onto the free wake nodes. By default, the velocities are only 
!    produced for the fixed wake grid for enough x stations to capture all of the free wake points.
!    If free wake points are observed to be outside the max x station currently calculated, the 
!    ifw value is incremented (adding x stations to the calculate list) until all of the free wake 
!    are included. If any point lies outside the fixed wake grid, it is calculated in the direct manner.

!                                                                       
!    UPDATE THE OLD WAKE VELOCITIES AND COMPUTE NEW ONES, IF NECESSARY  
!                                                                       
      NPW=0    
      NT1=NT-1                                                          
                                                                      
!                                                                       
!      UPDATE THE OLD NODE POINT VELOCITIES                             
!                                                                       
                                                               
      DO 10 I=1,NE                                                      
      DO 10 J=ntTerm,NT1                                                   
      UO(J,I)=U(J,I)+UFS(J,I)                                                  
      VO(J,I)=V(J,I)+VFS(J,I)                                                    
      WO(J,I)=W(J,I)+WFS(J,I)                                                    
      IF (X(J,I) .GE. XFW(1)) NPW=NPW+1                                 
10    CONTINUE  
                                                        
        ! Calculate freestream velocity at wake locations
	ygcErr=0
	do I=1,NE                                                      
		do J=ntTerm,NT1
			CALL CalcFreestream(Y(J,I),UFS(J,I),VFS(J,I),WFS(J,I),ygcErr) 	                                
		end do
	end do  
	                                                        
	if (ygcErr .ne. 0) then
		ierr1=1
		WRITE (6,607)  
	end if                                                   
                                                        
!                                                                       
!      IF VELOCITIES ARE NOT TO BE RECALCULATED, RETURN                 
!                                                                       
      IF (NT .NE. NSW) RETURN                                          
!                                                                       
!      UPDATE THE GRID POINT VELOCITIES                                 
!                                                                       
      XMAX=0.0                                                          
      YMAX=0.0                                                          
      YMIN=0.0                                                          
      ZMAX=0.0                                                          
      ZMIN=0.0                                                          
      JNPT=0                                                            
      KNPT=0                                                            
      DO 15 IW=1,IFW                                                    
      DO 15 JW=1,JFW                                                    
      DO 15 KW=1,KFW                                                    
      CALL PIVEL(NT,ntTerm,NBE,NB,NE,XFW(IW),YFW(JW),ZFW(KW),UFW(IW,JW,KW),VFW(IW,JW,KW),WFW(IW,JW,KW),1)                           
15    CONTINUE                                                          
!                                                                       
!      UPDATE THE NODE POINT VELOCITIES                                 
!                                                                       
      DO 75 I=1,NE                                                      
      DO 75 J=ntTerm,NT1                                                   
      XV=X(J,I)                                                         
      YV=Y(J,I)                                                         
      ZV=Z(J,I)                                                         
      ZAB=ABS(Z(J,I))                                                   
!                                                                       
!      IF NODE POINT IS OUTSIDE GRID, COMPUTE DIRECTLY                  
!                                                                       
      IF (XV .LT. XFW(1)) GO TO 70                                      
      XMAX=AMAX1(XMAX,XV)                                               
      IF (XV .GE. XFW(IFW)) GO TO 35                                    
      IF (YV .LT. YFW(1)) GO TO 60                                      
      IF (YV .GT. YFW(JFW)) GO TO 60                                    
      IF (ZAB .GT. ZFW(KFW)) GO TO 65                                   
      GO TO 50                                                          
35    CONTINUE                                                          
!                                                                       
!      IF BEYOND MAX X GRID POINT, INCREASE GRID IN X-DIRECTION         
!                                                                       
      IF (IFW .EQ. MaxFixWakeX) GO TO 70                                     
      IFW=IFW+1                                                         

!                                                                       
!      CALCULATE VELOCITIES FOR NEW GRID POINTS                         
!                                                                       
      DO 40 JW=1,JFW                                                    
      DO 40 KW=1,KFW                                                    
      CALL PIVEL(NT,ntTerm,NBE,NB,NE,XFW(IFW),YFW(JW),ZFW(KW),UFW(IFW,JW,KW),VFW(IFW,JW,KW),WFW(IFW,JW,KW),1)                       
40    CONTINUE                                                          
      IF (XV .GE. XFW(IFW)) GO TO 35                                    
      IF (YV .LT. YFW(1)) GO TO 70                                      
      IF (YV .GT. YFW(JFW)) GO TO 70                                    
      IF (ZAB .GT. ZFW(KFW)) GO TO 70                                   
50    CONTINUE                                                          
!                                                                       
!      GET INDUCED VELOCITY FROM STRAIGHT-LINE INTERPOLATION            
!      (INSIDE GRID)                                                    
!                                                                       
      CALL INTERP3(X(J,I),Y(J,I),Z(J,I),U(J,I),V(J,I),W(J,I))           
      GO TO 75                                                          
60    CONTINUE                                                          
      IF (XV .LE. XFW(1)) GO TO 70                                      
      JNPT=JNPT+1                                                       
      YMAX=AMAX1(YMAX,YV)                                               
      YMIN=AMIN1(YMIN,YV)                                               
      GO TO 70                                                          
65    CONTINUE                                                          
      IF (XV .LE. XFW(1)) GO TO 70                                      
      KNPT=KNPT+1                                                       
      ZMAX=AMAX1(ZMAX,ZV)                                               
      ZMIN=AMIN1(ZMIN,ZV)                                               
70    CONTINUE                                                          
!                                                                       
!      COMPUTE VORTEX INDUCED VELOCITY DIRECTLY                         
!                                                                       
      CALL PIVEL(NT,ntTerm,NBE,NB,NE,X(J,I),Y(J,I),Z(J,I),U(J,I),V(J,I),W(J,I),1)                                                
75    CONTINUE                                                          
      NSW=NT+IUT                                                       
RETURN                                                            
605   FORMAT(1H ,14HBUMP IFW. IFW=,I3,6H, JFW=,I4,6H, KFW=,I4,9H, IFWMAX=,I4,11H, XFW(IFW)=,F6.2,5H, XV=,F8.3,4H, I=,I3,4H, J=,I3) 
607   FORMAT(1H ,5X,92H*** (Y+YGC) LESS THAN ZERO FOR AT LEAST ONE NODE. URS SET EQUAL TO ZERO FOR THOSE NODES. ***)        
609   FORMAT(1H ,5HJNPT=,I4,7H, YMIN=,F6.3,7H, YMAX=,F6.3,7H, KNPT=,I3, 7H, ZMIN=,F6.3,7H, ZMAX=,F6.3,7H, XMAX=,F7.3)                     
END                                                               
