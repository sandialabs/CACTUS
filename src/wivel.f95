SUBROUTINE wivel(NT,ntTerm,NBE,NB,NE,IUT,NSW,NPW) 

      	use xwake  
      	use wakeloc
      	use vel
      	use veo
      	use test 
      	use freestream
      	use wallsoln
        use regtest 
      
        integer :: ygcErr 
        real :: Point(3), IndVel(3)                                                                
                                                                       
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
					Point=[X(J,I),Y(J,I),Z(J,I)]
					Call CalcIndVel(NT,ntTerm,NBE,NB,NE,Point,IndVel)
					U(J,I)=IndVel(1)
					V(J,I)=IndVel(2)
					W(J,I)=IndVel(3)
					
      				end do
			end do                                                                           

		end if

	end if
         
        ! Regression test 
        if (RegTFlag == 1) then
                Reg_MaxWVM=0.0
                ! Max wake velocity mag in first wake shed
                do I=1,NE
                        Reg_MaxWVM=max(Reg_MaxWVM,sqrt(U(1,I)**2+V(1,I)**2+W(1,I)**2))       
                end do
        end if
                                                     
RETURN  
607   FORMAT(1H ,5X,92H*** (Y+YGC) LESS THAN ZERO FOR AT LEAST ONE NODE. URS SET EQUAL TO ZERO FOR THOSE NODES. ***)                                                         
END                                                               
