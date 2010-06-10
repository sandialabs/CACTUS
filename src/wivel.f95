SUBROUTINE wivel(NT,ntTerm,NBE,NB,NE,IUT,NSW,NPW) 

      	use xwake  
      	use wakeloc
      	use vel
      	use veo
      	use test 
      	use freestream
      	use wallsoln
      
        integer :: ygcErr 
        real :: Point(3,1), IndVel(3,1)                                                                
                                                                       
	! Calculate the induced velocity at each lattice point in the wake from wake (including bound vorticity), wall, and freestream        
                                                                      
      	NPW=0                                                             
      	if (NT .ge. 1) then                                           
      	
		NT1=NT-1                                                          
                                                  
		! Update the old wake velocity values                                                                                
      		do I=1,NE                                                      
      			do J=ntTerm,NT1                                                   
      				UO(J,I)=U(J,I)+UFS(J,I)                                                  
      				VO(J,I)=V(J,I)+VFS(J,I)                                                    
      				WO(J,I)=W(J,I)+WFS(J,I)                                                    
      				if (X(J,I) .ge. XFW(1)) then
					NPW=NPW+1 
				end if
			end do                                
		end do                                                         
                                                                                                                                                                                          
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
                                                                                                                                                                           
		! If this is proper time step, calculate wall and wake influence on wake velocities                                                             
      		if (NT .eq. NSW) then                                         

			do I=1,NE                                                      
      				do J=ntTerm,NT1

					! Calculate wall and wake induced velocities at wake locations
					Point=reshape([X(J,I),Y(J,I),Z(J,I)],[3,1])
					Call CalcIndVel(NT,ntTerm,NBE,NB,NE,Point,IndVel)
					U(J,I)=IndVel(1,1)
					V(J,I)=IndVel(2,1)
					W(J,I)=IndVel(3,1)
					
      				end do
			end do                                                                           

		end if

	end if 
                                                     
RETURN  
607   FORMAT(1H ,5X,92H*** (Y+YGC) LESS THAN ZERO FOR AT LEAST ONE NODE. URS SET EQUAL TO ZERO FOR THOSE NODES. ***)                                                         
END                                                               
