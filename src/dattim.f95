SUBROUTINE DATTIM(DMY,HMS)
 
      character DMY*9, HMS*8

!  As written these are HPUX 7.0 non-standard calls.
      CALL DATE( DMY )
      CALL TIME( HMS )

RETURN
END
