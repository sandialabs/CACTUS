SUBROUTINE WriteWallData()   

	! Write wall data outputs

	use wallsoln 
	use configr

    implicit none

	integer :: i, OutIter  
    real :: TVel, dH 
    real, allocatable :: TVelIFS(:,:)   

    ! Setup
    if (NT==1) then

        ! Iteration on which to write (FS only, GP always last iter)
        OutIterFS=nti

        ! Initialize output counter
        OutCount=0

    end if

	! Write wake positions and velocity for each wake line on last rev
    if (irev == nr) then      
        OutCount=OutCount+1 

        if (GPFlag == 1) then
            ! Output ground plane source density averaged over the last revolution

            ! Average over last revolution
            do i=1,NumWP                                                  
                WSourceOut(i)=WSourceOut(i)+WSource(i,1)/real(nti)
            end do

            ! Write on last iter
            if (OutCount == nti) then
                do i=1,NumWP     
                    write(14,'(E13.7,",",$)') WCPoints(i,1) 
                    write(14,'(E13.7,",",$)') WCPoints(i,2) 
                    write(14,'(E13.7,",",$)') WCPoints(i,3) 
                    ! Dont suppress carriage return on last column
                    write(14,'(E13.7)') WSourceOut(i)
                end do
            end if

        end if

        if (FSFlag == 1) then
            ! Output wall tangent velocity and free surface height, dH/R=1/2*FnR^2*(1-(Ux@FS/UInf)^2)

            ! Write on appropriate iter
            if (OutCount == OutIterFS) then

                ! Calc tangential velocity induced by free surface 
                allocate(TVelIFS(NumFSCP,1))
                TVelIFS=matmul(FSInCoeffT,FSSource)

                ! Calc tangential velocity at each control point (already an average over last revolution)
                do i=1,NumFSCP

                    ! Get average tangential velocity on free surface
                    TVel=FSVTAve(i,1)+TVelIFS(i,1)

                    ! Calc dH (over radius)
                    dH=0.5*FnR**2*(1-TVel**2)

                    ! Write
                    write(15,'(E13.7,",",$)') FSCPPoints(i,1) 
                    write(15,'(E13.7,",",$)') FSCPPoints(i,2) 
                    write(15,'(E13.7,",",$)') FSCPPoints(i,3) 
                    write(15,'(E13.7,",",$)') TVel
                    ! Dont suppress carriage return on last column
                    write(15,'(E13.7)') dH

                end do

            end if

        end if

    end if

    Return
End SUBROUTINE WriteWallData
