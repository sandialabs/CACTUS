SUBROUTINE conlp()     

    use configr
	use blade

	logical NotDone

    ! Calculate the convection of shed vortex lattice points

	NT1=NT-1                                                          
	do I=1,NE                                                      
		if (NT > 1) then                                           
			do J=ntTerm,NT1                                                                                         

                if (iut == 1) then
					! Temporally "upwinded" calculation for the velocity to use in
					! convecting the wake (use velocity extrapolated to t=t+.5*dt)
					X(J,I)=X(J,I)+(3.0*(U(J,I)+UFS(J,I))-UO(J,I))*DT/2.0
					Y(J,I)=Y(J,I)+(3.0*(V(J,I)+VFS(J,I))-VO(J,I))*DT/2.0
					Z(J,I)=Z(J,I)+(3.0*(W(J,I)+WFS(J,I))-WO(J,I))*DT/2.0
				else
					! If wake velocities are not updated every timestep, simply
					! assume constant velocity over a single timestep.
					X(J,I)=X(J,I)+(U(J,I)+UFS(J,I))*DT
					Y(J,I)=Y(J,I)+(V(J,I)+VFS(J,I))*DT
					Z(J,I)=Z(J,I)+(W(J,I)+WFS(J,I))*DT
				end if

			end do
		end if

        ! Use straight integration for newest wake points for which there is no old data
		X(NT,I)=X(NT,I)+(U(NT,I)+UFS(NT,I))*DT                                  
		Y(NT,I)=Y(NT,I)+(V(NT,I)+VFS(NT,I))*DT                                        
		Z(NT,I)=Z(NT,I)+(W(NT,I)+WFS(NT,I))*DT                                        
	end do

    ! If ixterm is 1, update ntTerm as the farthest index in the wake that has at least one point below XSTOP
	if (IXTERM .eq. 1) then
		NotDone = .true.
		do while (NotDone .AND. ntTerm < NT1)
            ! Check
			do I=1,NE                                                      
				if (X(ntTerm,I) < XSTOP) then
					NotDone=.false.  
				end if
			end do

			if (NotDone) then
				ntTerm=ntTerm+1 
			end if
		end do
	end if

    Return                   
End SUBROUTINE conlp
