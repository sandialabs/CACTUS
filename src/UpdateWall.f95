subroutine UpdateWall(nt,ntTerm,nbe,nb,ne,iut,nsw)

      	use wakeloc
      	use wallsoln
      
        integer :: ygcErr, i
        real :: Point(3,1), dVel(3,1), NVelSum                                                               
                                                                       
	! Calculate the velocities at wall panels from wake (including bound vorticity), and freestream. 
	! Update wall RHS and calc new panel source strengths                                                                                                                             
	   
	                                                                                                
	do i=1,NumWP

		! Calculate freestream velocity at panel locations
		CALL CalcFreestream(WCPoints(2,i),dVel(1,1),dVel(2,1),dVel(3,1),ygcErr) 	                                
		NVelSum=sum(WZVec(1:3,i)*dVel(1:3,1))															
																	
		! If this is proper time step, update bound and wake vorticity influence on wall                                                           
		if (NT .eq. NSW) then                                         	
			! Calc wake induced velocity at wall panel locations                                                                            
			CALL PIVEL(NT,ntTerm,NBE,NB,NE,WCPoints(1,i),WCPoints(2,i),WCPoints(3,i),dVel(1,1),dVel(2,1),dVel(3,1),1)
			WWakeNVel(i)=sum(WZVec(1:3,i)*dVel(1:3,1))	
                                    
		end if
		NVelSum=NVelSum+WWakeNVel(i)
		
		! Set RHS 
		WRHS(i,1)=-NVelSum
	
  	end do
  	
  	! Calc new wall panel source strengths
  	WSource=matmul(WSMatI,WRHS)

return
end

