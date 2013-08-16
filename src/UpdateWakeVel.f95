SUBROUTINE UpdateWakeVel() 

    use configr
    use blade
    use regtest 

    integer :: ygcErr                                                                


	! Calculate the induced velocity at each lattice point in the wake from wake (including bound vorticity), wall, and freestream        

    if (NT .ge. 1) then                                           

		NT1=NT-1                                                          

        ! Update the old wake velocity values
        do I=1,NE                                                      
            do J=ntTerm,NT1                                                   
                UO(J,I)=U(J,I)+UFS(J,I)                                                  
                VO(J,I)=V(J,I)+VFS(J,I)                                                    
                WO(J,I)=W(J,I)+WFS(J,I)                                                    
            end do
		end do

        ! Calculate freestream velocity at wake locations
        ygcErr=0
		do I=1,NE                                                      
            do J=ntTerm,NT1
                CALL CalcFreestream(X(J,I),Y(J,I),UFS(J,I),VFS(J,I),WFS(J,I),ygcErr) 	                                
            end do
		end do

        ! If this is proper time step, calculate system influence on wake velocities
        if (NT .eq. NSW) then                                         

            do I=1,NE                                                      
                do J=ntTerm,NT1

                    ! Calculate wall and wake induced velocities at wake locations
                    Call CalcIndVel(NT,ntTerm,NBE,NB,NE,X(J,I),Y(J,I),Z(J,I),U(J,I),V(J,I),W(J,I))                                                                   

                end do
            end do

            ! Set the next update timestep                                        
            nsw=nt+iut 

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
END SUBROUTINE UpdateWakeVel
