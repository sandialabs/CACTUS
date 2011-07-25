SUBROUTINE csvwrite(FID,Header,Data,NRows)

        integer :: FID, NRows
        character(1000) :: Header
        real :: Data(:,:)
        
        integer :: nRow, nCol, i, j
        
        ! Writes a csv file to file specified by FID (assumed already opened with the open command).
        ! NRows is the number of rows of Data to write. If NRows<1, all rows of Data will be written...
        
        ! Write header
        write (FID,*) trim(Header)
        
        ! Write data
        if (NRows>0) then
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
10 format(E13.7,',',$) 
End
