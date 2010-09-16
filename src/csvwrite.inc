interface
        subroutine csvwrite(FID,Header,Data,NRows)
                integer :: FID, NRows
                character(1000) :: Header
                real :: Data(:,:)
        end subroutine
end interface