SUBROUTINE WriteWakeElementData() 

    ! Write wake element positions and velocity

    use wakedata
    use blade
    use wake
    use wallsoln 
    use configr
    use fnames
    
    implicit none

    integer :: tCount, tCountMax, wcount, node_id
    character(len=10) :: nt_str

    ! Optional wake element data output
    write(nt_str,'(I0)') nt
    WakeOutputFN=trim(FNBase)//'_WakeData_'//trim(nt_str)//'.csv'
    OPEN(12, FILE=WakeOutputFN)
    write(12,'(A)') trim(WakeOutHead)

    tCountMax=nt
    do wcount=1,NWakeInd+nb
        do tCount=1,tCountMax
            ! compute the unique ID number of the current wake filament
            ! (higher numbers correspond to newer elements)
            node_id = (NWakeInd+nb)*(tcount-1) + wcount

            write(12,'(E13.7,",",$)') TimeN                ! Normalized simulation time (t*Uinf/Rmax)
            write(12,'(I0,",",$)') node_id                 ! Unique node ID
            write(12,'(I0,",",$)') wcount                  ! Node number that wake element originated from
            write(12,'(E13.7,",",$)') X(tCount,wcount)     ! Wake element X position
            write(12,'(E13.7,",",$)') Y(tCount,wcount)     ! Wake element Y position
            write(12,'(E13.7,",",$)') Z(tCount,wcount)     ! Wake element Z position
            write(12,'(E13.7,",",$)') U(tCount,wcount)     ! Wake element X velocity
            write(12,'(E13.7,",",$)') V(tCount,wcount)     ! Wake element Y velocity
            ! Dont suppress carriage return on last column
            write(12,'(E13.7)') W(tCount,wcount)           ! Wake element Z velocity
        end do
    end do

    ! close the output file
    CLOSE(12)

    Return
End SUBROUTINE WriteWakeElementData
