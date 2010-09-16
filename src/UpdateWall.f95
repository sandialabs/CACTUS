subroutine UpdateWall()

        use configr
      	use blade
      	use wallsoln
        use regtest 
      
        integer :: ygcErr, i
        real :: dVel(3), NVelSum                                                               
                                                                       
	! Calculate the velocities at wall panels from wake (including bound vorticity), and freestream. 
	! Update wall RHS and calc new panel source strengths                                                                                                                             
	   
	                                                                                                
	do i=1,NumWP

		! Calculate freestream velocity at panel locations
		CALL CalcFreestream(WCPoints(i,2),dVel(1),dVel(2),dVel(3),ygcErr) 	                                
		NVelSum=sum(WZVec(i,1:3)*dVel)															
                				
		! If this is proper time step, update bound and wake vorticity influence on wall                                                           
		if (NT .eq. NSW) then                                         	
			! Calc wake induced velocity at wall panel locations                                                                            
			CALL PIVEL(NT,ntTerm,NBE,NB,NE,WCPoints(i,1),WCPoints(i,2),WCPoints(i,3),dVel(1),dVel(2),dVel(3),1)
			WWakeNVel(i)=sum(WZVec(i,1:3)*dVel)                    	
                                    
		end if
		NVelSum=NVelSum+WWakeNVel(i)
		
		! Set RHS 
		WRHS(i,1)=-NVelSum
	
  	end do
  	
  	! Calc new wall panel source strengths
  	WSource=matmul(WSMatI,WRHS)
        
        ! Regression Test
        if (RegTFlag == 1) then
                Reg_MaxWS=maxval(abs(WSource))
        end if 

return
end

