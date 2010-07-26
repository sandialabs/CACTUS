SUBROUTINE LISTIT (IOPT,IDEV)     

      use dystl
      use varscale
      use cltab
      use shear
      use airfoil
      use configr
      use ioption
      use pidef
      use vortex
      use time
      use test                                            
                                    
      DIMENSION ALPOS(2),ALNEG(2),NSECTON(2)                            
      CHARACTER NSECTON*6
      DATA NSECTON/'SINGLE','DOUBLE'/                                   
!                                                                       
!    TELL THE USER WHAT THIS TURBINE CONFIGURATION IS AND WHAT THE      
!    OPERATING CONDITIONS ARE                                           
!                                                                       
      YREF1=YREF*RMAX                                                   
      YGC1=YGC*RMAX                                                                                                        
      WRITE (IDEV,601)                                      
      WRITE (IDEV,602) DMY,HMS,JBTITLE                                 
      WRITE (IDEV,603) AREAT   
      
      IF (IFWG .NE. 0) WRITE (IDEV,604)                                 
      WRITE (IDEV,605) VRAD                                             
      IF (IXTERM .NE. 0) WRITE (IDEV,606) XSTOP                         
      IF (ifc .EQ. 0) WRITE (IDEV,608) IUT,NTI  
      IF (ifc .NE. 0) WRITE (IDEV,614) IUT,NTI,IUTF,NTIF                                                      
      IF (NSECT .GT. 1) WRITE (IDEV,616) (AFTITLE(ISECT(J)),J=2,6)
      IF(NSECT.EQ.2) WRITE(IDEV,618) AFTITLE(1),AFTITLE(2),DFTITLE(1),DFTITLE(2)                 
      IF(NSECT.EQ.1) WRITE(IDEV,619) AFTITLE(1),DFTITLE(1)
      ISTP1=NSTL(1)                                                     
      ISTP2=NSTL(2)                                                     
      I=1                                                               
10    CONTINUE                                                          
      DO 15 KK=1,NSECT                                                  
      IF (I .GT. NSTL(KK)) GO TO 15                                     
      ALPOS(KK)=ALSTLP(I,KK)*CONDEG                                     
      ALNEG(KK)=ALSTLN(I,KK)*CONDEG                                     
15    CONTINUE                                                          
      IF(I.LE.ISTP1) WRITE(IDEV,623) RESTL(I,1),ALPOS(1),ALNEG(1)       
      IF(I.LE.ISTP1.AND.I.LE.ISTP2) WRITE(IDEV,624) RESTL(I,2),ALPOS(2),ALNEG(2)                                                        
      IF(I.LE.ISTP2.AND.I.GT.ISTP1) WRITE(IDEV,625) RESTL(I,2),ALPOS(2),ALNEG(2)                                                        
      IF(I.GE.ISTP1.AND.I.GE.ISTP2) GO TO 20                            
      I=I+1                                                             
      GO TO 10                                                          
20    CONTINUE                                                          
      ALZER1=ALZER(1)*CONDEG                                            
      IF (NSECT .EQ. 2) GO TO 25                                        
      WRITE (IDEV,628) TC(1),ALZER1                                     
      IF (CAMB(1) .NE. 0.0) WRITE (IDEV,630) CAMBER(1)                  
      GO TO 40                                                          
25    CONTINUE                                                          
      ALZER2=ALZER(2)*CONDEG                                            
      WRITE (IDEV,629) TC(1),TC(2),ALZER1,ALZER2                        
      IF (CAMB(1) .EQ. 0.0) GO TO 35                                    
      WRITE(IDEV,630) CAMBER(1)                                         
      IF (CAMB(2) .NE. 0.0) WRITE(IDEV,631) CAMBER(2)                   
      GO TO 40                                                          
35    CONTINUE                                                          
      IF (CAMB(2) .NE. 0.0) WRITE(IDEV,632) CAMBER(2)                   
40    CONTINUE                                                          
      WRITE (IDEV,634) RHO,VIS,SLEX,TEMPR,YGC1,YREF1                     
      IF (IOPT .NE. 0) GO TO 50                                         
                                                        
      WRITE (6,601)                                                     
      RETURN                                                            
50    CONTINUE                                                          
      IF (IREV .LT. 1) GO TO 80                                         
      WRITE (IDEV,640)                                                  
      DO 70 I=1,IREV                                                    
      WRITE (IDEV,642) I,CP(I),KP(I),DTIME(I),ETIME(I),CPF(I),KPF(I)    
70    CONTINUE                                                          
      WRITE (IDEV,645) POWER,uMPH                                       
80    CONTINUE                                                          
      IF (IERR0 .NE. 0) WRITE (IDEV,635)                                
      IF (IERR1 .NE. 0) WRITE (IDEV,637)                                
      WRITE (IDEV,650) TimeF                                          
RETURN                                                            
601   FORMAT('1')                                                       
602   FORMAT(' ',39X,'TURBINE CODE',/,' ',46X,'EXECUTED ON ',A9,', AT ',A8,/,' ',24X,'JOB: ',A80,/,' ')       
603   FORMAT(' ',2X,'TURBINE PROJECTED AREA IS ',F10.2,' SQ FT',/,'0',2X,'OPTIONS IN EFFECT',/,' ')  
 	                       
604   FORMAT(' ',5X,'FIXED WAKE GRID.')                                
605   FORMAT(' ',5X,'FINITE CORE VORTEX WITH',F6.4,' CORE RADIUS-TO-MAX ROTOR RADIUS RATIO.')                                           
606   FORMAT(' ',5X,'IGNORE ALL NODE POINTS MORE THAN ',F5.2,' RADII DOWNSTREAM OF TOWER CENTERLINE.')                        
608   FORMAT(' ',5X,'RECALCULATE WAKE VELOCITIES EVERY ',I2,' TIME STEPS (',I2,' THETA INCREMENTS PER REVOLUTION).')        

614   FORMAT(' ',5X,'THIS CALCULATION INVOLVES TWO THETA INCREMENT SIZES.',/,' ',6X,'RECALCULATE WAKE VELOCITIES EVERY ',I2,' TIME STEPS, INITIALLY (',I2, &
	' THETA INCREMENTS PER REVOLUTION).',/,' ',6X,'THEN RECALCULATE WAKE VELOCITIES EVERY ',I2,' TIME STEPS (',I2,' THETA INCREMENTS PER REVOLUTION).')                            
616   FORMAT('0',2X,'BLADE SEGMENT SECTION PROFILES ARE',/,' ',5X,'SEG 1-',4A4,'  SEG 2-',4A4,'  SEG 3-',4A4,'  SEG 4-',4A4,'  SEG 5-',4A4,/,' ')                                               
618   FORMAT(' ',2X,'AIRFOIL STALL CHARACTERISTICS ARE (ANGLES IN DEGREES)',//,' ',26X,4A4,50X,4A4,/,' ',2X,16A4,2X,16A4,/,'0','REYNOLDS NO.',  &
	5X,'POS. STALL ANGLE',5X,'NEG. STALL ANGLE',12X,'REYNOLDS NO.',5X,'POS. STALL ANGLE',5X,'NEG. STALL ANGLE',/,' ')                
619   FORMAT(' ',2X,'AIRFOIL STALL CHARACTERISTICS ARE (ANGLES IN DEGREES)',//,' ',26X,4A4,/,' ',2X,16A4,/,'0','REYNOLDS NO.',5X,'POS. STALL ANGLE',5X,'NEG. STALL ANGLE',/,' ')                   
623   FORMAT(' ',3X,F10.0,11X,F5.1,15X,F5.1)                            
624   FORMAT('+',69X,F10.0,11X,F5.1,15X,F5.1)                           
625   FORMAT(' ',69X,F10.0,11X,F5.1,15X,F5.1)                           
628   FORMAT('0',1X,'THICKNESS-TO-CHORD RATIO IS ',F5.3/2X,'ZERO LIFT ANGLE IS ',F5.2,' DEGREES.')                                      
629   FORMAT('0',1X,'THICKNESS-TO-CHORD RATIO IS ',F5.3,33X,'THICKNESS-TO-CHORD RATIO IS ',F5.3,/,2X,'ZERO LIFT ANGLE IS ',F5.2,' DEGREES.',33X,'ZERO LIFT ANGLE IS ',F5.2,' DEGREES.')                    
630   FORMAT('0',1X,'THIS AIRFOIL IS CAMBERED ',A3,' WITH RESPECT TO THE TOWER.')                                                       
631   FORMAT('+',67X,'THIS AIRFOIL IS CAMBERED ',A3,' WITH RESPECT TO THE TOWER.')                                                      
632   FORMAT('0',67X,'THIS AIRFOIL IS CAMBERED ',A3,' WITH RESPECT TO THE TOWER.')                                                      
634   FORMAT('0',15X,'ATMOSPHERIC MODEL DATA ARE',/,20X,'AIR DENSITY=',F8.6,' SLUG/FT3',/,20X,'ABS. VISCOSITY=',E12.5,' LBF-SEC/FT2',/,20X,'SHEAR EXPONENT=',F4.2,/,20X, &
	'AMBIENT TEMPERATURE=',F5.0,' DEGREES F.',/,20X,'GROUND CLEARANCE HT.=',F5.1,' FT.',/,20X,'REFERENCE HT.=',F5.1,' FT.')                                    
635   FORMAT('0',5X,'***** ALREFM IS BAD FOR AT LEAST ONE BLADE SEGMENT. *****')                    
637   FORMAT('0',5X,'***** (Y+YGC) IS LESS THAN ZERO FOR AT LEAST ONE CALCULATION. CHECK PRINTOUT FOR DETAILS. *****')                   
640   FORMAT('0',15X,'REV',7X,'CP',13X,'KP',11X,'DELTA CPU TIME',4X,'TOTAL CPU TIME',5X,'CPF',12X,'KPF',/,' ')                      
642   FORMAT(' ',13X,I5,2E15.5,8X,F8.2,10X,F8.2,2E15.5)                 
645   FORMAT('0',15X,'TURBINE POWER IS ',F5.0,' KW AT WIND VELOCITY OF ',F6.2,' MPH.')                                                  
650   FORMAT('0',/,'0',10X,'THIS EXECUTION REQUIRED ',F8.2,' SECONDS OF CPU TIME.')                                                       
1101  FORMAT('TITLE=',20A4,', TSR=',F4.1,', RPM=',F4.0)                 
END                                                               
