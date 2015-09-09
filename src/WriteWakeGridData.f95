SUBROUTINE WriteWakeGridData()

    ! Write wake grid data

    use wakedata
    use blade
    use wake
    use wallsoln
    use configr
    use fnames

    implicit none

    real :: xnode, ynode, znode
    integer :: ygcErr
    character(len=10) :: nt_str

    ! Open file for writing - a new file at each timestep
    write(nt_str,'(I5.5)') nt
    WakeDefOutputFN=trim(FNBase)//'_WakeDefData_'//trim(nt_str)//'.csv'
    OPEN(13, FILE=WakeDefOutputFN)
    write(13,'(A)') trim(GridVelOutHead)

    !! Compute wake data on specified cartesian grid
    ! Compute blade, wake, and wall induced streamwise velocity deficit
    do zcount=1,nzgrid
!$omp parallel do private(xcount,xnode,ynode,znode)
        do ycount=1,nygrid
            do xcount=1,nxgrid
                ! Get the grid node location
                xnode = XGrid(xcount,ycount,zcount)
                ynode = YGrid(xcount,ycount,zcount)
                znode = ZGrid(xcount,ycount,zcount)

                ! Calculate wall and wake induced velocities at grid locations
                Call CalcIndVel(NT,ntTerm,NBE,NB,NE, &
                    xnode,ynode,znode, &
                    VXInd(xcount,ycount,zcount),VYInd(xcount,ycount,zcount),VZInd(xcount,ycount,zcount))

                ! Calculate free stream velocities at grid locations
                Call CalcFreestream(xnode,ynode,znode, &
                    UfsGrid(xcount,ycount,zcount),VfsGrid(xcount,ycount,zcount),WfsGrid(xcount,ycount,zcount), &
                    ygcErr)
            end do
        end do
!$omp end parallel do
    end do


    ! Output blade, wake, and wall induced streamwise velocity deficit
    do zcount=1,nzgrid
        do ycount=1,nygrid
            do xcount=1,nxgrid
                ! Write to file
                write(13,'(E13.7,",",$)') TimeN                            ! Normalized simulation time (t*Uinf/Rmax)
                write(13,'(E13.7,",",$)') XGrid(xcount,ycount,zcount)
                write(13,'(E13.7,",",$)') YGrid(xcount,ycount,zcount)      ! grid node locations
                write(13,'(E13.7,",",$)') ZGrid(xcount,ycount,zcount)
                write(13,'(E13.7,",",$)') VXInd(xcount,ycount,zcount)
                write(13,'(E13.7,",",$)') VYInd(xcount,ycount,zcount)      ! induced velocities
                write(13,'(E13.7,",",$)') VZInd(xcount,ycount,zcount)
                write(13,'(E13.7,",",$)') UfsGrid(xcount,ycount,zcount)
                write(13,'(E13.7,",",$)') VfsGrid(xcount,ycount,zcount)    ! freestream velocities
                write(13,'(E13.7)'      ) WfsGrid(xcount,ycount,zcount)    ! Dont suppress carriage return on last column
            end do
        end do
    end do

    ! close the output file
    CLOSE(13)

    Return
End SUBROUTINE WriteWakeGridData
