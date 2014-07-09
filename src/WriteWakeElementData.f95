SUBROUTINE WriteWakeElementData() 

    ! Write wake element positions and velocity

    use wakedata
    use blade
    use wallsoln 
    use configr
    
    implicit none

    integer :: tCount, tCountMax, wcount

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
            write(12,'(E13.7)') W(tCount,WakeLineInd(wcount)) ! Dont suppress carriage return on last column
        end do
    end do


Return
End SUBROUTINE WriteWakeElementData
