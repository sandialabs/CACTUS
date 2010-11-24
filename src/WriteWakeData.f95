SUBROUTINE WriteWakeData()   

	! Write wake data outputs

	use wakedata
	use blade
	use wallsoln 
	use configr
        
	integer :: tCount, tCountMax, wcount

        ! Write header
        if (NT==1) then
                write(12,*) trim(WakeOutHead)
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
        end if     
				
Return
End	    
