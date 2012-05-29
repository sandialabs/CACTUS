MODULE airfoil
	
	! Airfoil section data
	
	character*80, allocatable :: aftitle(:)  	! Title for each airfoil section
	integer, allocatable :: camb(:)			! Camber flag for each section (
	real, allocatable :: tc(:)			! Thickness to chord ratio for each section 
        real, allocatable :: alzer(:)                   ! Zero lift AOA for each section    
             
        ! Airfoil section coefficient data
        real, allocatable :: TA(:,:,:)          ! Table AOA values
        real, allocatable :: TCL(:,:,:)         ! Table CL values
        real, allocatable :: TCD(:,:,:)         ! Table CD values
        real, allocatable :: TCM(:,:,:)         ! Table Cm values
        real, allocatable :: TRE(:,:)           ! Table Re values
        integer, allocatable :: nTBL(:,:)       ! Number of AOA values for each Re number, in each section data table
        integer, allocatable :: nRET(:)         ! Number of Re number values in each section data table       
            
        ! Interpolation warning flags    
        integer :: ilxtp
        integer :: iuxtp
             
        ! Airfoil params for BV dyn stall     
        real, allocatable :: alstlp(:,:)                ! Stall AOA (positive) at all Re numbers for each section
        real, allocatable :: alstln(:,:)                ! Stall AOA (negative) at all Re numbers for each section    
             
        ! Airfoil params for LB dyn stall
        real, allocatable :: CLaData(:,:)
        real, allocatable :: CLCritPData(:,:)
	real, allocatable :: CLCritNData(:,:)

	CONTAINS

	SUBROUTINE airfoil_cns(MaxAOAVals,MaxReVals,MaxAirfoilSect)

		! Constructor for the arrays in this module

		integer :: MaxAOAVals,MaxReVals,MaxAirfoilSect
		
		allocate(aftitle(MaxAirfoilSect))
		allocate(camb(MaxAirfoilSect))
		allocate(tc(MaxAirfoilSect))
                allocate(alzer(MaxAirfoilSect))    
                allocate(TA(MaxAOAVals,MaxReVals,MaxAirfoilSect))
                allocate(TCL(MaxAOAVals,MaxReVals,MaxAirfoilSect))
                allocate(TCD(MaxAOAVals,MaxReVals,MaxAirfoilSect))
                allocate(TCM(MaxAOAVals,MaxReVals,MaxAirfoilSect))              
                allocate(TRE(MaxReVals,MaxAirfoilSect))         
                allocate(nTBL(MaxReVals,MaxAirfoilSect))
                allocate(nRET(MaxAirfoilSect)) 
                allocate(alstlp(MaxReVals,MaxAirfoilSect))
                allocate(alstln(MaxReVals,MaxAirfoilSect))         
                allocate(CLaData(MaxReVals,MaxAirfoilSect))
                allocate(CLCritPData(MaxReVals,MaxAirfoilSect))
                allocate(CLCritNData(MaxReVals,MaxAirfoilSect))              	
		
	End SUBROUTINE


        SUBROUTINE intp(RE,ALPHA,CL,CD,CM25,KK)   
                                
            real :: RE, ALPHA, CL, CD, CM25
            integer :: KK                    
            
            real :: CLA(2),CDA(2),CM25A(2)                                      
            integer :: U1, X1, iUB, iLB
            logical :: NotDone                                               
                                                                                
            !    INTERPOLATE ON RE NO. AND ANGLE OF ATTACK TO GET AIRFOIL CHARACTERISTICS                                            
                                                                                
            CLA(:)=0.0                                                        
            CDA(:)=0.0  
            CM25A(:)=0.0
                                                                
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
                CLA(I)=TCL(L1,J,KK)+XA*(TCL(U1,J,KK)-TCL(L1,J,KK))                
                CDA(I)=TCD(L1,J,KK)+XA*(TCD(U1,J,KK)-TCD(L1,J,KK))    
                CM25A(I)=TCM(L1,J,KK)+XA*(TCM(U1,J,KK)-TCM(L1,J,KK)) 
                        
                I=I+1           
            end do                                                          
                                                                            
            ! DO STRAIGHT LINE INTERPOLATION ON RE NO.                         
                                                                        
            CL=CLA(1)+XRE*(CLA(2)-CLA(1))                                     
            CD=CDA(1)+XRE*(CDA(2)-CDA(1))  
            CM25=CM25A(1)+XRE*(CM25A(2)-CM25A(1))
                                                                   
        END SUBROUTINE                                                               

        
        Subroutine CalcBVStallAOALim(Re,SectInd,alssp,alssn)
        
            ! Get stall AOA limits for BV model from airfoil data
            integer :: SectInd
            real :: Re, alssp, alssn
            
            integer :: iUB, iLB
            logical :: NotDone 
            
            ! Find Re upper and lower bounds.                                     
                
            if (RE >= TRE(1,SectInd)) then                                                                                                                                        
                
                NotDone=.true.    
                iUB=2                                                                 
                do while (NotDone)   
                        
                    if (RE <= TRE(iUB,SectInd)) then
                        ! Done
                        NotDone=.false.
                        if (RE == TRE(iUB,SectInd)) then
                            iLB=iUB
                        else
                            iLB=iUB-1                                                           
                            XRE=(RE-TRE(iLB,SectInd))/(TRE(iUB,SectInd)-TRE(iLB,SectInd))
                        end if
                    else
                        if (iUB == nRET(SectInd)) then       
                            ! No upper bound in table, take last point...
                            NotDone=.false.                                                       
                            iLB=iUB                                                           
                            XRE=0.0                                                           
                        else    
                            ! No upper bound, increment and continue                                
                            iUB=iUB+1
                        end if
                    end if

                end do                                                    
            
            else        
                ! No lower bound in table, take first point.                                            
                iLB=1                                                             
                iUB=1                                                             
                XRE=0.0                                                                                                 
            end if
            
            ! Interp
            alssp=alstlp(iLB,SectInd)+xRE*(alstlp(iUB,SectInd)-alstlp(iLB,SectInd))            
            alssn=alstln(iLB,SectInd)+xRE*(alstln(iUB,SectInd)-alstln(iLB,SectInd))  
            
        End Subroutine
        
        
        Subroutine CalcLBStallAOALim(Re,SectInd,CLa,CLCritP,CLCritN)
        
            ! Get stall data for LB model from airfoil data
            integer :: SectInd
            real :: Re, CLa, CLCritP, CLCritN
            
            integer :: iUB, iLB
            logical :: NotDone 
            
            ! Find Re upper and lower bounds.                                     
                
            if (RE >= TRE(1,SectInd)) then                                                                                                                                        
                
                NotDone=.true.    
                iUB=2                                                                 
                do while (NotDone)   
                        
                    if (RE <= TRE(iUB,SectInd)) then
                        ! Done
                        NotDone=.false.
                        if (RE == TRE(iUB,SectInd)) then
                            iLB=iUB
                        else
                            iLB=iUB-1                                                           
                            XRE=(RE-TRE(iLB,SectInd))/(TRE(iUB,SectInd)-TRE(iLB,SectInd))
                        end if
                    else
                        if (iUB == nRET(SectInd)) then       
                            ! No upper bound in table, take last point...
                            NotDone=.false.                                                       
                            iLB=iUB                                                           
                            XRE=0.0                                                           
                        else    
                            ! No upper bound, increment and continue                                
                            iUB=iUB+1
                        end if
                    end if

                end do                                                    
            
            else        
                ! No lower bound in table, take first point.                                            
                iLB=1                                                             
                iUB=1                                                             
                XRE=0.0                                                                                                 
            end if
            
            ! Interp
            CLa=CLaData(iLB,SectInd)+xRE*(CLaData(iUB,SectInd)-CLaData(iLB,SectInd))            
            CLCritP=CLCritPData(iLB,SectInd)+xRE*(CLCritPData(iUB,SectInd)-CLCritPData(iLB,SectInd))  
            CLCritN=CLCritNData(iLB,SectInd)+xRE*(CLCritNData(iUB,SectInd)-CLCritNData(iLB,SectInd)) 
            
        End Subroutine
End