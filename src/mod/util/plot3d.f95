module plot3d

    ! plot3d : subroutines for reading plot3d mesh data

    implicit none

    private
    public read_p3d_multiblock

    integer, parameter :: IOunit = 20

contains

    subroutine read_p3d_multiblock(xyz_filename,nblocks,ni,nj,nk,x,y,z)
        ! read_p3d_multiblock() : Read a formatted (ASCII) Plot3D multi-block mesh file.
        !   Coordinates must be specified one number per line.
        !
        !   Inputs
        !   ======
        !   xyz_filename : Plot3d mesh filename (*.xyz)
        !
        !   Outputs
        !   =======
        !   nblocks      : number of blocks
        !   ni           : mesh dimensions in first direction
        !   nj           : mesh dimensions in second direciton
        !   nk           : mesh dimensions in third direction

        character(len=*), intent(in) :: xyz_filename                             ! Plot3d mesh filename (*.xyz)

        integer, intent (out)             :: nblocks                             ! number of blocks
        integer, allocatable, intent(out) :: ni(:), nj(:), nk(:)                 ! mesh dimensions
        real, allocatable, intent(out)    :: x(:,:,:,:), y(:,:,:,:), z(:,:,:,:)  ! mesh coordinates

        integer :: i,j,k
        integer :: m,n
        integer :: nimax, njmax, nkmax

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
        do  m = 1, nblocks
            read(IOunit,*) &
            ((( x(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
            ((( y(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m)), &
            ((( z(i,j,k,m), i=1,ni(m)), j=1,nj(m)), k=1,nk(m))
         enddo

        ! close file
        close(IOunit)

    end subroutine read_p3d_multiblock


end module plot3d