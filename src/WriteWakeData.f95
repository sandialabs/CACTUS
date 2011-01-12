SUBROUTINE WriteWakeData()   

	! Write wake data outputs

	use wakedata
	use blade
	use wallsoln 
	use configr
        
        implicit none
        
	integer :: tCount, tCountMax, wcount

        ! Write header
        if (NT==1) then
                write(12,*) trim(WakeOutHead)
                
                ! JCM test: wake deficit setup
                ntcount=0
                
                dxgrid=(xgridU-xgridL)/(nxgrid-1)
                dzgrid=(zgridU-zgridL)/(nzgrid-1)
                do xcount=1,nxgrid                                                    
                        do zcount=1,nzgrid
                                XGrid(xcount,zcount)=xgridL+(xcount-1)*dxgrid
                                ZGrid(xcount,zcount)=zgridL+(zcount-1)*dzgrid
                                SVDef(xcount,zcount)=0.0
                        end do
                end do
        end if

	! Write wake positions and velocity for each wake line on last rev
        if (irev == nr) then      
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
                 
                
                if (WakeOutFlag > 1) then
                
                        ! Output blade, wake, and wall induced streamwise velocity deficit on a plane.
                        
                        ! Averaged over last revolution
                        do xcount=1,nxgrid                                                    
                                do zcount=1,nzgrid
        
                                        ! Calculate wall and wake induced velocities at grid
                                        PointGrid=[XGrid(xcount,zcount),ygrid,ZGrid(xcount,zcount)]
                                        Call CalcIndVel(NT,ntTerm,NBE,NB,NE,PointGrid,IndVelGrid)
                                        SVDef(xcount,zcount)=SVDef(xcount,zcount)+IndVelGrid(1)/nti
                                        
                                end do
                        end do
                        ntcount=ntcount+1
                        
                        ! Write on last iter
                        if (ntcount == nti) then
                        do xcount=1,nxgrid        
                                        do zcount=1,nzgrid  
                                                if (zcount < nzgrid) then
                                                        write(13,'(E13.7,",",$)') SVDef(xcount,zcount)
                                                else 
                                                        ! Dont suppress carriage return on last column
                                                        write(13,'(E13.7)') SVDef(xcount,zcount)
                                                end if
                                        end do
                                end do 
                        end if
                
                end if
                 
        end if     
				
Return
End	    
