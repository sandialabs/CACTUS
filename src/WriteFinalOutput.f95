SUBROUTINE WriteFinalOutput()

    use util
    use airfoil, only : ilxtp, iuxtp

    implicit none

    ! Check error flags and write notifications to stdout
    if (ilxtp .gt. 0) write (6,615)
    if (iuxtp .gt. 0) write (6,618)

    Return
618 FORMAT ('AT LEAST ONE BLADE REYNOLDS NUMBER WAS ABOVE TABLE LIMIT.  UPPER LIMIT USED.')
615 FORMAT ('AT LEAST ONE BLADE REYNOLDS NUMBER WAS BELOW TABLE LIMIT.  LOWER LIMIT USED.')
End SUBROUTINE WriteFinalOutput
