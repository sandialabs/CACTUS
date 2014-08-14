SUBROUTINE WriteWakeElementData() 

    ! Write wake element positions and velocity

    use wakedata
    use blade
    use wake
    use wallsoln 
    use configr
    
    implicit none

    integer :: tCount, tCountMax, wcount

    tCountMax=NT
    do wcount=1,NWakeInd 
        do tCount=1,tCountMax
            write(12,'(E13.7,",",$)') TimeN                             ! Normalized simulation time (t*Uinf/Rmax)                                   ! Timestep number
            write(12,'(I0,",",$)') WakeLineInd(wcount)                  ! Blade element number that wake element originated from
            write(12,'(E13.7,",",$)') X(tCount,WakeLineInd(wcount))     ! Wake element X position
            write(12,'(E13.7,",",$)') Y(tCount,WakeLineInd(wcount))     ! Wake element Y position
            write(12,'(E13.7,",",$)') Z(tCount,WakeLineInd(wcount))     ! Wake element Z position
            write(12,'(E13.7,",",$)') U(tCount,WakeLineInd(wcount))     ! Wake element X velocity
            write(12,'(E13.7,",",$)') V(tCount,WakeLineInd(wcount))     ! Wake element Y velocity
            ! Dont suppress carriage return on last column
            write(12,'(E13.7)') W(tCount,WakeLineInd(wcount))           ! Wake element Z velocity
        end do
    end do


Return
End SUBROUTINE WriteWakeElementData
