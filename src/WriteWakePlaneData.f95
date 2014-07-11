SUBROUTINE WriteWakePlaneData()

    ! Write wake plane data
    
    use wakedata
    use blade
    use wallsoln 
    use configr
    
    implicit none
    
    integer :: xcount, ycount, zcount
    real :: dxgrid, dygrid, dzgrid       
    real :: vx, vy, vz 
 
    ! Set up grids for plane wake output (only on first iteration)
    if (NT==1) then
        if (WakePlaneOutFlag==1) then
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
        else if (WakePlaneOutFlag==2) then
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
        else if (WakePlaneOutFlag==3) then
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

    ! Compute and write wake data on specified cartesian grid
    if (WakePlaneOutFlag == 1) then
        ! Output blade, wake, and wall induced streamwise velocity deficit on a horizontal plane.
        do zcount=1,nzhgrid
            do xcount=1,nxhgrid                                                    
                ! Calculate wall and wake induced velocities at grid
                Call CalcIndVel(NT,ntTerm,NBE,NB,NE,XHGrid(xcount,zcount),yhgrid,ZHGrid(xcount,zcount),vx,vy,vz)
                VXIndH(xcount,zcount) = vx
                VYIndH(xcount,zcount) = vy
                VZIndH(xcount,zcount) = vz

                ! Write to file
                write(13,'(I8,",",$)') NT
                write(13,'(E13.7,",",$)') XHGrid(xcount,zcount) 
                write(13,'(E13.7,",",$)') ZHGrid(xcount,zcount) 
                write(13,'(E13.7,",",$)') VXIndH(xcount,zcount)
                write(13,'(E13.7,",",$)') VYIndH(xcount,zcount)  
                write(13,'(E13.7)') VZIndH(xcount,zcount) ! Dont suppress carriage return on last column
            end do
        end do

    else if (WakePlaneOutFlag == 2) then
        ! Output blade, wake, and wall induced streamwise velocity deficit on a vertical plane.
        do ycount=1,nyvgrid
            do xcount=1,nxvgrid                                                    
                ! Calculate wall and wake induced velocities at grid
                Call CalcIndVel(NT,ntTerm,NBE,NB,NE,XVGrid(xcount,ycount),YVGrid(xcount,ycount),zvgrid,vx,vy,vz)
                VXIndV(xcount,ycount) = vx
                VYIndV(xcount,ycount) = vy
                VZIndV(xcount,ycount) = vz

                ! Write to file
                write(13,'(I8,",",$)') NT
                write(13,'(E13.7,",",$)') XVGrid(xcount,ycount) 
                write(13,'(E13.7,",",$)') YVGrid(xcount,ycount) 
                write(13,'(E13.7,",",$)') VXIndV(xcount,ycount)
                write(13,'(E13.7,",",$)') VYIndV(xcount,ycount)  
                write(13,'(E13.7)') VZIndV(xcount,ycount) ! Dont suppress carriage return on last column
            end do
        end do

    else if (WakePlaneOutFlag == 3) then
        ! Output blade, wake, and wall induced streamwise velocity deficit on a cross-section plane.
        do zcount=1,nzcgrid
            do ycount=1,nycgrid
                ! Calculate wall and wake induced velocities at grid
                Call CalcIndVel(NT,ntTerm,NBE,NB,NE,xcgrid,YCGrid(ycount,zcount),ZCGrid(ycount,zcount),vx,vy,vz)
                VXIndC(ycount,zcount) = vx
                VYIndC(ycount,zcount) = vy
                VZIndC(ycount,zcount) = vz

                ! Write to file
                write(13,'(I8,",",$)') NT
                write(13,'(E13.7,",",$)') YCGrid(ycount,zcount)
                write(13,'(E13.7,",",$)') ZCGrid(ycount,zcount)
                write(13,'(E13.7,",",$)') VXIndC(ycount,zcount)
                write(13,'(E13.7,",",$)') VYIndC(ycount,zcount)
                write(13,'(E13.7)') VZIndC(ycount,zcount) ! Dont suppress carriage return on last column

            end do
        end do  
    end if

    Return
End SUBROUTINE WriteWakePlaneData
