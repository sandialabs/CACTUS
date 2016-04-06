Module fnames

    use pathseparator

    implicit none

    integer :: nargin, FNLength, status, DecInd
    character(80) :: InputFN, SFOutputFN, RevOutputFN, TSOutputFN, ELOutputFN, RegOutputFN, WakeOutputFN, FieldOutputFN, GPOutputFN, FSOutputFN, DSOutputFN, FNBase
    character(255) :: OutputPath,WakeElementOutputPath,WallOutputPath,FieldOutputPath,ProbeOutputPath

contains

    subroutine get_FNBase()
        logical :: back

        Call get_command_argument(1,InputFN,FNLength,status)

        back=.true.
        FNBase=InputFN((index(InputFN,'/',back)+1):len(InputFN))
        DecInd=index(FNBase,'.',back)

        if (DecInd > 1) then
            FNBase=FNBase(1:(DecInd-1))
        end if

    end subroutine get_FNBase

End Module fnames
