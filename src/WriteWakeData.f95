SUBROUTINE WriteWakeData()   

	! Write wake data outputs

	use wakedata
	use blade
	use wallsoln 
	use configr

    implicit none

	integer :: tCount, tCountMax, wcount, xcount, ycount, zcount
    real :: dxgrid, dygrid, dzgrid       
    real :: vx, vy, vz 

    ! Wake deficit setup
    if (NT==1) then

        ntcount=0

        if (WakeOutFlag==2) then

            ! Setup horizontal grid
            dxgrid=(xhgridU-xhgridL)/(nxhgrid-1)
            dzgrid=(zhgridU-zhgridL)/(nzhgrid-1)
            do xcount=1,nxhgrid                                                    
                do zcount=1,nzhgrid
                    XHGrid(xcount,zcount)=xhgridL+(xcount-1)*dxgrid
                    ZHGrid(xcount,zcount)=zhgridL+(zcount-1)*dzgrid
                    VXIndH(xcount,zcount)=0.0
                    VYIndH(xcount,zcount)=0.0
                    VZIndH(xcount,zcount)=0.0
                end do
            end do
        else if (WakeOutFlag==3) then

            ! Setup vertical grid
            dxgrid=(xvgridU-xvgridL)/(nxvgrid-1)
            dygrid=(yvgridU-yvgridL)/(nyvgrid-1)
            do xcount=1,nxvgrid                                                    
                do ycount=1,nyvgrid
                    XVGrid(xcount,ycount)=xvgridL+(xcount-1)*dxgrid
                    YVGrid(xcount,ycount)=yvgridL+(ycount-1)*dygrid
                    VXIndV(xcount,ycount)=0.0
                    VYIndV(xcount,ycount)=0.0
                    VZIndV(xcount,ycount)=0.0
                end do
            end do
        else if (WakeOutFlag==4) then

            ! Setup cross-section grid
            dygrid=(ycgridU-ycgridL)/(nycgrid-1)
            dzgrid=(zcgridU-zcgridL)/(nzcgrid-1)
            do ycount=1,nycgrid
                do zcount=1,nzcgrid
                    YCGrid(ycount,zcount)=ycgridL+(ycount-1)*dygrid
                    ZCGrid(ycount,zcount)=zcgridL+(zcount-1)*dzgrid
                    VXIndC(ycount,zcount)=0.0
                    VYIndC(ycount,zcount)=0.0
                    VZIndC(ycount,zcount)=0.0
                end do
            end do
        end if
    end if

	! Write wake positions and velocity for each wake line on last 5 revs
      if (irev > nr-nRevWake) then
        tCountMax=NT
        do wcount=1,NWakeInd 
            do tCount=1,tCountMax
	          write(12,'(I8,",",$)') NT
                write(12,'(I8,",",$)') WakeLineInd(wcount)
                write(12,'(E13.7,",",$)') X(tCount,WakeLineInd(wcount)) 
                write(12,'(E13.7,",",$)') Y(tCount,WakeLineInd(wcount))
                write(12,'(E13.7,",",$)') Z(tCount,WakeLineInd(wcount))
                write(12,'(E13.7,",",$)') U(tCount,WakeLineInd(wcount)) 
                write(12,'(E13.7,",",$)') V(tCount,WakeLineInd(wcount))
                ! Dont suppress carriage return on last column
                write(12,'(E13.7)') W(tCount,WakeLineInd(wcount))
            end do
        end do
  end if
	  

	if (irev > nr-nRevWake) then
        if (WakeOutFlag == 2) then

            ! Output blade, wake, and wall induced streamwise velocity deficit on a plane.

            ! Averaged over last revolution
            do zcount=1,nzhgrid
                do xcount=1,nxhgrid                                                    

                    ! Calculate wall and wake induced velocities at grid
                    Call CalcIndVel(NT,ntTerm,NBE,NB,NE,XHGrid(xcount,zcount),yhgrid,ZHGrid(xcount,zcount),vx,vy,vz)
                    VXIndH(xcount,zcount)=VXIndH(xcount,zcount)+vx/nti
                    VYIndH(xcount,zcount)=VYIndH(xcount,zcount)+vy/nti
                    VZIndH(xcount,zcount)=VZIndH(xcount,zcount)+vz/nti

                end do
            end do
            ntcount=ntcount+1

            ! Write on last iter
            if (ntcount == nti) then
            do zcount=1,nzhgrid  
                do xcount=1,nxhgrid        

                        write(13,'(E13.7,",",$)') XHGrid(xcount,zcount) 
                        write(13,'(E13.7,",",$)') ZHGrid(xcount,zcount) 
                        write(13,'(E13.7,",",$)') VXIndH(xcount,zcount)
                        write(13,'(E13.7,",",$)') VYIndH(xcount,zcount)  
                        ! Dont suppress carriage return on last column
                        write(13,'(E13.7)') VZIndH(xcount,zcount)
                    end do
                end do
            end if

        else if (WakeOutFlag == 3) then

            ! Output blade, wake, and wall induced streamwise velocity deficit on a plane.

            ! Averaged over last revolution
            do ycount=1,nyvgrid
                do xcount=1,nxvgrid                                                    
                    ! Calculate wall and wake induced velocities at grid
                    Call CalcIndVel(NT,ntTerm,NBE,NB,NE,XVGrid(xcount,ycount),YVGrid(xcount,ycount),zvgrid,vx,vy,vz)
                    VXIndV(xcount,ycount)= vx
                    VYIndV(xcount,ycount)= vy
                    VZIndV(xcount,ycount)= vz
			  write(13,'(I8,",",$)') NT
                    write(13,'(E13.7,",",$)') XVGrid(xcount,ycount) 
                    write(13,'(E13.7,",",$)') YVGrid(xcount,ycount) 
                    write(13,'(E13.7,",",$)') VXIndV(xcount,ycount)
                    write(13,'(E13.7,",",$)') VYIndV(xcount,ycount)  
                    ! Dont suppress carriage return on last column
                    write(13,'(E13.7)') VZIndV(xcount,ycount)
                end do
            end do
            ntcount=ntcount+1

            ! ! Write on last iter
!             if (ntcount == nti) then
!                 do xcount=1,nxvgrid        
!                     do ycount=1,nyvgrid  
! 				write(13,'(I8,",",$)') NT
!                         write(13,'(E13.7,",",$)') XVGrid(xcount,ycount) 
!                         write(13,'(E13.7,",",$)') YVGrid(xcount,ycount) 
!                         write(13,'(E13.7,",",$)') VXIndV(xcount,ycount)
!                         write(13,'(E13.7,",",$)') VYIndV(xcount,ycount)  
!                         ! Dont suppress carriage return on last column
!                         write(13,'(E13.7)') VZIndV(xcount,ycount)
!                     end do
!                 end do
!             end if

        else if (WakeOutFlag == 4) then

            ! Output blade, wake, and wall induced streamwise velocity deficit on a plane.

            ! Averaged over last revolution
            do zcount=1,nzcgrid
                do ycount=1,nycgrid

                    ! Calculate wall and wake induced velocities at grid
                    Call CalcIndVel(NT,ntTerm,NBE,NB,NE,xcgrid,YCGrid(ycount,zcount),ZCGrid(ycount,zcount),vx,vy,vz)
                    VXIndC(ycount,zcount)=VXIndC(ycount,zcount)+vx/nti
                    VYIndC(ycount,zcount)=VYIndC(ycount,zcount)+vy/nti
                    VZIndC(ycount,zcount)=VZIndC(ycount,zcount)+vz/nti

                end do
            end do
            ntcount=ntcount+1

            ! Write on last iter
            if (ntcount == nti) then
            do zcount=1,nzcgrid
                do ycount=1,nycgrid
                        write(13,'(E13.7,",",$)') YCGrid(ycount,zcount)
                        write(13,'(E13.7,",",$)') ZCGrid(ycount,zcount)
                        write(13,'(E13.7,",",$)') VXIndC(ycount,zcount)
                        write(13,'(E13.7,",",$)') VYIndC(ycount,zcount)
                        ! Dont suppress carriage return on last column
                        write(13,'(E13.7)') VZIndC(ycount,zcount)
                    end do
                end do
            end if

        end if

    end if

    Return
End SUBROUTINE WriteWakeData
