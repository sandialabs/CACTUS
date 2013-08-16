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
            write (FID,*) trim(Header)
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


End MODULE Util
