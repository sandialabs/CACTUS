subroutine UpdateWall()

        use configr
      	use blade
      	use wallsoln
        use regtest 
      
        integer :: ygcErr, i, IBCInd
        real :: Point(3), dVel(3), NVelSum, TVelSum, dUdX, dUdXSum                                                              
        
        ! If this is a wall update timestep...                                                               
        if (nt == nsWall) then                                                               

            if (GPFlag == 1) then
                ! Ground plane                                                               
                
                ! Calculate the velocities at wall panels from wake (including bound vorticity), and freestream. 
                ! Update wall RHS and calc new panel source strengths
                do i=1,NumWP
        
                    ! Calculate freestream velocity at panel locations
                    CALL CalcFreestream(WCPoints(i,2),dVel(1),dVel(2),dVel(3),ygcErr) 	                                
                    NVelSum=sum(WZVec(i,1:3)*dVel)															
                    
                    ! Calc wake induced velocity at wall panel locations                                                                            
                    CALL BladeIndVel(NT,ntTerm,NBE,NB,NE,WCPoints(i,1),WCPoints(i,2),WCPoints(i,3),dVel(1),dVel(2),dVel(3),dUdX,0,0)                  	
                    NVelSum=NVelSum+sum(WZVec(i,1:3)*dVel)
                    
                    ! Calc FS induced velocity
                    if (FSFlag == 1) then
                        Point=[WCPoints(i,1),WCPoints(i,2),WCPoints(i,3)]
                        Call FSIndVel(Point,0,dVel,dUdX)
                        NVelSum=NVelSum+sum(WZVec(i,1:3)*dVel)
                    end if
                    
                    ! Set RHS 
                    WRHS(i,1)=-NVelSum
                
                end do
                
                ! Calc new wall panel source strengths
                WSource=matmul(WSMatI,WRHS)
                
            end if

            
            if (FSFlag == 1) then
                ! Free Surface                                                               
                
                ! Calculate the velocities at wall panels from wake (including bound vorticity), and freestream. 
                ! Update wall RHS and calc new panel source strengths
                if (UseFSWall) then
                    ! Wall BC
                    do i=1,NumFSCP
            
                        ! Calculate freestream velocity at panel locations
                        CALL CalcFreestream(FSCPPoints(i,2),dVel(1),dVel(2),dVel(3),ygcErr)                                     
                        NVelSum=sum(FSCZVec(i,1:3)*dVel)    
                        TVelSum=sum(FSCXVec(i,1:3)*dVel)            
                        
                        ! Calc wake induced velocity at wall panel locations                                                                            
                        CALL BladeIndVel(NT,ntTerm,NBE,NB,NE,FSCPPoints(i,1),FSCPPoints(i,2),FSCPPoints(i,3),dVel(1),dVel(2),dVel(3),dUdX,0,0)                          
                        NVelSum=NVelSum+sum(FSCZVec(i,1:3)*dVel)
                        TVelSum=TVelSum+sum(FSCXVec(i,1:3)*dVel)
    
                        ! Calc GP induced velocity
                        if (GPFlag == 1) then
                            Point=[FSCPPoints(i,1),FSCPPoints(i,2),FSCPPoints(i,3)]
                            Call GPIndVel(Point,0,dVel,dUdX)
                            NVelSum=NVelSum+sum(FSCZVec(i,1:3)*dVel)
                            TVelSum=TVelSum+sum(FSCXVec(i,1:3)*dVel)
                        end if
                        
                        ! Set RHS in list for running averages
                        FSRHS(i,FSRHSInd)=-NVelSum
                        ! Set new average RHS value. Note: Must build average over at least one revolution...
                        FSRHSAve(i,1)=sum(FSRHS(i,1:NFSRHSAve))/real(NFSRHSAve)
                        
                        ! Set induced tangent velocity at colocation points in list for running averages (used for output)
                        FSVT(i,FSRHSInd)=TVelSum
                        ! Set new average VT value. Note: Must build average over at least one revolution...
                        FSVTAve(i,1)=sum(FSVT(i,1:NFSRHSAve))/real(NFSRHSAve)
                                            
                    end do
                                    
                else
                    ! Free surface BC
                    do i=1,NumFSCP
            
                        ! Calculate freestream velocity at panel locations
                        CALL CalcFreestream(FSCPPoints(i,2),dVel(1),dVel(2),dVel(3),ygcErr) 	                                
                        NVelSum=sum(FSCZVec(i,1:3)*dVel)	
                        TVelSum=sum(FSCXVec(i,1:3)*dVel)	
                        dUdXSum=0.0
                        
                        ! Calc wake induced velocity at wall panel locations                                                                            
                        CALL BladeIndVel(NT,ntTerm,NBE,NB,NE,FSCPPoints(i,1),FSCPPoints(i,2),FSCPPoints(i,3),dVel(1),dVel(2),dVel(3),dUdX,0,1)                   	
                        NVelSum=NVelSum+sum(FSCZVec(i,1:3)*dVel)
                        TVelSum=TVelSum+sum(FSCXVec(i,1:3)*dVel)
                        dUdXSum=dUdXSum+dUdX
    
                        ! Calc GP induced velocity
                        if (GPFlag == 1) then
                            Point=[FSCPPoints(i,1),FSCPPoints(i,2),FSCPPoints(i,3)]
                            Call GPIndVel(Point,1,dVel,dUdX)
                            NVelSum=NVelSum+sum(FSCZVec(i,1:3)*dVel)
                            TVelSum=TVelSum+sum(FSCXVec(i,1:3)*dVel)
                            dUdXSum=dUdXSum+dUdX
                        end if
                        
                        ! Set RHS in list for running averages
                        FSRHS(i,FSRHSInd)=(FnR**2)*dUdXSum-NVelSum
                        ! Set new average RHS value. Note: Must build average over at least one revolution...
                        FSRHSAve(i,1)=sum(FSRHS(i,1:NFSRHSAve))/real(NFSRHSAve)
                        
                        ! Set BC RHS if this is a BC colocation point
                        if (FSBCRow(i)>0) then
                            FSRHS(FSBCRow(i),FSRHSInd)=-dUdXSum
                            FSRHSAve(FSBCRow(i),1)=sum(FSRHS(FSBCRow(i),1:NFSRHSAve))/real(NFSRHSAve)
                        end if
                        
                        ! Set induced tangent velocity at colocation points in list for running averages (used for output)
                        FSVT(i,FSRHSInd)=TVelSum
                        ! Set new average VT value. Note: Must build average over at least one revolution...
                        FSVTAve(i,1)=sum(FSVT(i,1:NFSRHSAve))/real(NFSRHSAve)
                    
                    end do
                    
                end if
                
                ! Calc new free surface panel source strengths with running average (over approx 360 deg turbine rotation)
                ! RHS values. 
                FSSource=matmul(FSSMatI,FSRHSAve) 
                
                ! Reset running average index (once entire buffer is filled, overwrites old data at beginning of average buffer)
                if (FSRHSInd<NFSRHSAve) then
                    FSRHSInd=FSRHSInd+1
                else
                    FSRHSInd=1
                end if         

            end if
            
            ! Set next wall update timestep
            nsWall=nt+iWall
        
        end if
        
        ! Regression Test
        if (RegTFlag == 1) then
            Reg_MaxWS=maxval(abs(WSource))
        end if 

return
end

