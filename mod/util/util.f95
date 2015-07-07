MODULE util

    ! Utilities used in other modules and elsewhere in the code

CONTAINS

    SUBROUTINE cross(ax,ay,az,bx,by,bz,cx,cy,cz) 

        real ax,ay,az,bx,by,bz,cx,cy,cz 

        cx = ay*bz - az*by
        cy = az*bx - ax*bz
        cz = ax*by - ay*bx

    End SUBROUTINE cross


    subroutine QuatRot(vx,vy,vz,Theta,nRx,nRy,nRz,Ox,Oy,Oz,vRx,vRy,vRz)

        ! % Perform rotation of vector v around normal vector nR using the
        ! % quaternion machinery.
        ! % v: input vector
        ! % Theta: rotation angle (rad)
        ! % nR: normal vector around which to rotate
        ! % Origin: origin point of rotation
        ! %
        ! % vR: Rotated vector

        real :: vx,vy,vz,Theta,nRx,nRy,nRz,Ox,Oy,Oz,vRx,vRy,vRz     

        real :: p(4,1), pR(4,1), q(4), qbar(4), nRMag, vOx, vOy, vOz
        real :: QL(4,4), QbarR(4,4)

        ! Force normalize nR
        nRMag=sqrt(nRx**2+nRy**2+nRz**2)
        nRx=nRx/nRMag
        nRy=nRy/nRMag
        nRz=nRz/nRMag

        ! Quaternion form of v
        vOx=vx-Ox
        vOy=vy-Oy
        vOz=vz-Oz
        p=reshape((/0.0,vOx,vOy,vOz/),(/4,1/))

        ! Rotation quaternion and conjugate
        q=(/cos(Theta/2),nRx*sin(Theta/2),nRy*sin(Theta/2),nRz*sin(Theta/2)/)
        qbar=(/q(1),-q(2),-q(3),-q(4)/)

        QL=transpose(reshape((/q(1), -q(2), -q(3), -q(4), &
            q(2),  q(1), -q(4),  q(3), &
            q(3),  q(4),  q(1), -q(2), &
            q(4), -q(3),  q(2),  q(1)/),(/4,4/)))

        QbarR=transpose(reshape((/qbar(1), -qbar(2), -qbar(3), -qbar(4), &
            qbar(2),  qbar(1),  qbar(4), -qbar(3), &
            qbar(3), -qbar(4),  qbar(1),  qbar(2), &
            qbar(4),  qbar(3), -qbar(2),  qbar(1)/),(/4,4/)))

        ! Rotate p
        pR=matmul(matmul(QbarR,QL),p)
        vRx=pR(2,1)+Ox
        vRy=pR(3,1)+Oy
        vRz=pR(4,1)+Oz

    end subroutine QuatRot


    SUBROUTINE csvwrite(FID,Header,Data,WriteHead,NRows)

        integer :: FID, WriteHead, NRows
        character(10000) :: Header
        real :: Data(:,:)

        integer :: nRow, nCol, i, j

        ! Writes comma separated data to file specified by FID (assumed already opened with the open command).
        ! If WriteHead is 1, will write the input header line, else will skip header line.
        ! NRows is the number of rows of Data to write. If NRows<0, all rows of Data will be written...

        ! Write header
        if (WriteHead>0) then
            write (FID,'(A)') trim(Header)
        end if

        ! Write data
        if (NRows>=0) then
            nRow=NRows
        else
            nRow=size(Data,1)
        end if
        nCol=size(Data,2)
        do i=1,nRow
            do j=1,nCol
                if (j<nCol) then
                    write(FID,10) Data(i,j)
                else
                    write(FID,'(E13.7)') Data(i,j)
                end if
            end do
        end do

        Return
10      format(E13.7,',',$) 
    End SUBROUTINE Csvwrite

    SUBROUTINE file_to_stdout(filename)
    ! file_to_stdout() : Read a file in and write it's output to stdout.
    !   (Useful for appending input files to the output)

    integer, parameter :: max_linewidth = 1000    
    character(80) :: filename
    integer :: reason

    character(max_linewidth) :: readline

    ! open file for reading
    open(30, file=filename)

    ! read in line by line
    do
        read(30, '(A)', iostat=reason) readline

        ! if end of file is reached, break out
        if (reason < 0) then
            exit
        else
            write(*,*) trim(readline)
        end if
    end do

    ! close file
    close(30)
    
    End SUBROUTINE file_to_stdout

    SUBROUTINE read_p3d_grid(xyz_filename,x,y,z)
    ! read_p3d_surf() : read a formatted (ASCII) Plot3D multi-block grid file.
    !   Coordinates must be specified one number per line.

        ! maximum grid dimensions
        integer, parameter :: imax  = 200
        integer, parameter :: jmax  = 200
        integer, parameter :: kmax  = 200   ! for surfaces, can probably set this to 1, since surfaces have k-index = 1
        integer, parameter :: nbmax = 100

        ! set the IOunit
        integer, parameter :: IOunit = 20

        integer :: i
        integer :: j
        integer :: m
        integer :: n
        integer :: nblocks
        integer :: ni(nbmax)
        integer :: nj(nbmax)
        integer :: nk(nbmax)

        character(len=*), intent(in) :: xyz_filename 

        real, intent(out) :: x(imax,jmax,kmax,nbmax)
        real, intent(out) :: y(imax,jmax,kmax,nbmax)
        real, intent(out) :: z(imax,jmax,kmax,nbmax)

        ! open file
        open (unit=IOunit, form='formatted', file=xyz_filename)

        ! read number of blocks
        read(IOunit,*) nblocks
        
        ! read dimensions of each blocks
        do m=1,nblocks
            read(IOunit,*) ni(m),nj(m),nk(m)
        end do
        
        ! read in data for each block
        do  m=1,nblocks
            do i=1,ni(m)
                do j=1,nj(m)
                    do k=1,nk(m)
                        read(IOunit,*) x(i,j,k,m)
                    end do
                end do
            end do

            do i=1,ni(m)
                do j=1,nj(m)
                    do k=1,nk(m)
                        read(IOunit,*) y(i,j,k,m)
                    end do
                end do
            end do


            do i=1,ni(m)
                do j=1,nj(m)
                    do k=1,nk(m)
                        read(IOunit,*) z(i,j,k,m)
                    end do
                end do
            end do
        enddo

        return

        ! close file
        close(IOunit)

    end SUBROUTINE read_p3d_surf

End MODULE Util
