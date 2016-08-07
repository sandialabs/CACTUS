module vecutils

    ! module vec : functions for performing vector operations

    implicit none

    private
    public cross3, mag3, calcrotation3

contains

    pure function cross3(a,b)

        ! cross3() : returns the cross product a x b, with a(3) and b(3)
        !
        !   Inputs:
        !   =======
        !   a,b      : real(3) vectors
        !
        !   Outputs:
        !   ========
        !   cross(3) : real(3) vector

        real, intent(in) :: a(3), b(3)
        real :: cross3(3)

        cross3(1) = a(2) * b(3) - a(3) * b(2)
        cross3(2) = a(3) * b(1) - a(1) * b(3)
        cross3(3) = a(1) * b(2) - a(2) * b(1)

    end function cross3


    pure function mag3(a)

        ! mag3() : returns the magnitude of a(3)
        !
        !   Inputs:
        !   =======
        !   a        : real(3) vector
        !
        !   Outputs:
        !   ========
        !   mag      : real, length of a

        real, intent(in) :: a(3)
        real             :: mag3

        mag3 = sqrt(sum(a**2))

    end function mag3


    subroutine calcrotation3(R,VecI,VecO,reverse)

        real    :: R(3,3), VecI(3), VecO(3)
        integer :: reverse

        ! Apply rotation matrix R to VecI to get VecO
        ! Reverse: 1 to use the transpose of R

        if (reverse == 0) then
            VecO(1)=R(1,1)*VecI(1)+R(1,2)*VecI(2)+R(1,3)*VecI(3)
            VecO(2)=R(2,1)*VecI(1)+R(2,2)*VecI(2)+R(2,3)*VecI(3)
            VecO(3)=R(3,1)*VecI(1)+R(3,2)*VecI(2)+R(3,3)*VecI(3)
        else
            VecO(1)=R(1,1)*VecI(1)+R(2,1)*VecI(2)+R(3,1)*VecI(3)
            VecO(2)=R(1,2)*VecI(1)+R(2,2)*VecI(2)+R(3,2)*VecI(3)
            VecO(3)=R(1,3)*VecI(1)+R(2,3)*VecI(2)+R(3,3)*VecI(3)
        end if

    end subroutine calcrotation3


end module vecutils