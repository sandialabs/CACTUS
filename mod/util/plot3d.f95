module plot3d
    
    ! plot3d : subroutines for reading plot3d mesh data

    implicit none

    integer, parameter :: IOunit = 20

contains

    subroutine read_p3d_multiblock(xyz_filename,nblocks,ni,nj,nk,x,y,z)
        ! read_p3d_multiblock() : Read a formatted (ASCII) Plot3D multi-block mesh file.
        !   Coordinates must be specified one number per line.
        !
        !   Format Example:
        !   [nblocks]
        !   [nx] [ny] [nz]
        !   [x_1]
        !   [x_2]
        !   ...
        !   [x_nx]
        !   [y_1]
        !   [y_2]
        !   ...
        !   [y_ny]
        !   [z_1]
        !   [z_2]
        !   ...
        !   [z_nz]
        !   EOF

        character(len=*), intent(in) :: xyz_filename 

        integer :: i,j,k
        integer :: m,n
        integer :: nimax, njmax, nkmax

        integer, intent (out)             :: nblocks                             ! number of blocks (out)
        integer, allocatable, intent(out) :: ni(:), nj(:), nk(:)                 ! mesh dimensions (out)
        real, allocatable, intent(out)    :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)  ! mesh coordinates (out)

        ! open file
        open(unit=IOunit, form='formatted', file=xyz_filename)

        ! read number of blocks
        read(IOunit,*) nblocks
        
        ! allocate for the dimensions of the mesh blocks
        allocate(ni(nblocks))
        allocate(nj(nblocks))
        allocate(nk(nblocks))

        ! read dimensions of each block
        do m=1,nblocks
            read(IOunit,*) ni(m),nj(m),nk(m)
        end do

        ! compute maximum mesh size and allocate
        nimax = maxval(ni)
        njmax = maxval(nj)
        nkmax = maxval(nk)

        allocate(x(nimax,njmax,nkmax,nblocks))
        allocate(y(nimax,njmax,nkmax,nblocks))
        allocate(z(nimax,njmax,nkmax,nblocks))
        
        ! read in data for each block
        do  m=1,nblocks
            do k=1,nk(m)
                do j=1,nj(m)
                    do i=1,ni(m)
                        read(IOunit,*) x(i,j,k,m)
                    end do
                end do
            end do

            do k=1,nk(m)
                do j=1,nj(m)
                    do i=1,ni(m)
                        read(IOunit,*) y(i,j,k,m)
                    end do
                end do
            end do

            do k=1,nk(m)
                do j=1,nj(m)
                    do i=1,ni(m)
                        read(IOunit,*) z(i,j,k,m)
                    end do
                end do
            end do
        enddo

        return

        ! close file
        close(IOunit)

    end subroutine read_p3d_multiblock

end module plot3d