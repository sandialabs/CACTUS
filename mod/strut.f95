MODULE strut

	! Strut geometry and loads outputs

    use util

    implicit none

    type StrutType
        integer :: NElem
        real :: TtoC  ! Strut thickness to chord ratio

        ! Strut element end locations
        real, allocatable :: MCx(:)
        real, allocatable :: MCy(:)
        real, allocatable :: MCz(:)
        ! Strut chord to radius at element ends
        real, allocatable :: CtoR(:)
        ! Strut element center locations
        real, allocatable :: PEx(:)
        real, allocatable :: PEy(:)
        real, allocatable :: PEz(:)
        ! Strut element spanwise vector
        real, allocatable :: sEx(:)
        real, allocatable :: sEy(:)
        real, allocatable :: sEz(:)
        ! Strut element chord to radius and norm. area
        real, allocatable :: ECtoR(:)
        real, allocatable :: EAreaR(:)
        ! Blade and element indicies to which the strut connects (for interference drag calc)
        ! For struts that are attached to the rotor shaft at one end (not to another blade), set the appropriate BInd and EInd values to zero.
        ! BIndS, EIndS : first strut element
        ! BIndE, EIndE : last strut element 
        integer :: BIndS
        integer :: EIndS
        integer :: BIndE
        integer :: EIndE

        real :: LR      ! Strut length to radius ratio
        ! Interference drag calc parameters
        real :: tcS      ! Thickness to chord of the blade element at the strut-blade junction (first strut element)
        real :: tcE      ! Thickness to chord of the blade element at the strut-blade junction (last strut element)

        ! Current flow quantities at each element center
        real, allocatable :: ReStrut(:) 
        real, allocatable :: u(:) ! u velocity over Uinf
        real, allocatable :: v(:) ! v velocity over Uinf
        real, allocatable :: w(:) ! w velocity over Uinf
        real, allocatable :: ur(:) ! velocity mag over Uinf

        ! Current strut element coeffs (normalized by strut element scale parameters)
        real, allocatable :: Cd0(:)

        ! Current total strut output (nomalized by machine scale parameters)
        real :: CP  ! Power coefficient due to this strut
        real :: CTR ! Torque coefficient due to this strut
        real :: CFx ! Fx coefficient due to this strut
        real :: CFy ! Fy coefficient due to this strut
        real :: CFz ! Fz coefficient due to this strut

    end type StrutType

    type(StrutType), allocatable :: Struts(:)   

    integer :: NStrut ! number of struts

    real, parameter :: ReCrit=3.0e5
    real :: Cdpar   ! Additional strut parasitic interference drag coefficient based on "chord area" (chord squared)

    ! Current sum of output over all struts (nomalized by machine scale parameters)
    real :: CP_S  ! Power coefficient due to all struts
    real :: CTR_S ! Torque coefficient due to all struts
    real :: CFx_S ! Fx coefficient due to all struts
    real :: CFy_S ! Fy coefficient due to all struts
    real :: CFz_S ! Fz coefficient due to all struts


CONTAINS


    SUBROUTINE strut_comp_cns(SInd,NElem)

        implicit none

        ! Constructor for the arrays for a strut component

        integer :: SInd, NElem

        Struts(SInd)%NElem=NElem
        allocate(Struts(SInd)%MCx(NElem+1))
        allocate(Struts(SInd)%MCy(NElem+1))
        allocate(Struts(SInd)%MCz(NElem+1))
        allocate(Struts(SInd)%CtoR(NElem+1))
        allocate(Struts(SInd)%PEx(NElem))
        allocate(Struts(SInd)%PEy(NElem))
        allocate(Struts(SInd)%PEz(NElem))
        allocate(Struts(SInd)%sEx(NElem))
        allocate(Struts(SInd)%sEy(NElem))
        allocate(Struts(SInd)%sEz(NElem))
        allocate(Struts(SInd)%ECtoR(NElem))
        allocate(Struts(SInd)%EAreaR(NElem))
        allocate(Struts(SInd)%ReStrut(NElem)) 
        allocate(Struts(SInd)%u(NElem))
        allocate(Struts(SInd)%v(NElem))
        allocate(Struts(SInd)%w(NElem))
        allocate(Struts(SInd)%ur(NElem))
        allocate(Struts(SInd)%Cd0(NElem))         

    End SUBROUTINE strut_comp_cns


    SUBROUTINE RotateStrut(SNum,delt,nrx,nry,nrz,px,py,pz)

        implicit none

        integer :: SNum, j, NElem
        real :: delt,nrx,nry,nrz,px,py,pz
        real :: vrx,vry,vrz

        ! Rotates data in strut arrays. Rotate element end geometry and recalculate element geometry.

        NElem=Struts(SNum)%NElem

        ! Strut end locations
        do j=1,NElem+1
            Call QuatRot(Struts(SNum)%MCx(j),Struts(SNum)%MCy(j),Struts(SNum)%MCz(j),delt,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
            Struts(SNum)%MCx(j)=vrx
            Struts(SNum)%MCy(j)=vry
            Struts(SNum)%MCz(j)=vrz
        end do

        ! Calc element geometry
        Call CalcSEGeom(SNum)


    End SUBROUTINE RotateStrut


    SUBROUTINE CalcSEGeom(SNum)

        implicit none

        integer :: SNum, NElem, j
        real :: sE(3)
        real :: sEM

        ! Calculates element geometry from element end geometry

        NElem=Struts(SNum)%NElem

        do j=1,NElem

            ! Element center locations
            Struts(SNum)%PEx(j)=0.5*(Struts(SNum)%MCx(j+1)+Struts(SNum)%MCx(j))
            Struts(SNum)%PEy(j)=0.5*(Struts(SNum)%MCy(j+1)+Struts(SNum)%MCy(j))
            Struts(SNum)%PEz(j)=0.5*(Struts(SNum)%MCz(j+1)+Struts(SNum)%MCz(j))

            ! Set spanwise vectors
            sE=(/Struts(SNum)%MCx(j+1)-Struts(SNum)%MCx(j),Struts(SNum)%MCy(j+1)-Struts(SNum)%MCy(j),Struts(SNum)%MCz(j+1)-Struts(SNum)%MCz(j)/)
            sEM=sqrt(dot_product(sE,sE))
            sE=sE/sEM
            Struts(SNum)%sEx(j)=sE(1)
            Struts(SNum)%sEy(j)=sE(2)
            Struts(SNum)%sEz(j)=sE(3)

		    ! Calc element area and chord
		    Struts(SNum)%ECtoR(j)=0.5*(Struts(SNum)%CtoR(j+1)+Struts(SNum)%CtoR(j))
		    Struts(SNum)%EAreaR(j)=sEM*Struts(SNum)%ECtoR(j)

        end do


    End SUBROUTINE CalcSEGeom


    SUBROUTINE StrutElemCoeffs(SInd,EInd)

        implicit none

        integer :: SInd, EInd
        real :: st, ReS, vs, AOS, pc
        real :: Cflam, Cdlam, Cfturb, Cdturb, Fblend

        ! Updates strut element coeffs for current flow states

        ! Calc sideslip angle
        vs=Struts(SInd)%u(EInd)*Struts(SInd)%sEx(EInd)+Struts(SInd)%v(EInd)*Struts(SInd)%sEy(EInd)+Struts(SInd)%w(EInd)*Struts(SInd)%sEz(EInd)
        AOS=asin(vs/Struts(SInd)%ur(EInd)) 

        ! Effective section thickness and Re referenced to nominal flow path length rather than chord (p/c = 1/cos(AOS))
        ! Allows estimation at high element sideslip (azimuthal struts).
        pc=1.0/abs(cos(AOS))
        pc=min(pc,Struts(SInd)%LR/Struts(SInd)%ECtoR(EInd)) ! Limit p to strut length
        st=Struts(SInd)%TtoC/pc
        ReS=Struts(SInd)%ReStrut(EInd)*pc

        ! Calculate strut element profile drag
        Cflam = 2.66 / SQRT(ReS)  ! Laminar friction drag coefficient
        Cdlam = 2.0 * Cflam * (1 + st) + st**2 ! Laminar drag coefficient
        Cfturb = 0.044 / ReS**(1.0/6.0) ! Turbulent friction drag coefficient
        Cdturb = 2.0 * Cfturb * (1.0 + 2.0*st + 60.0*st**4) ! Turbulent drag coefficient
        Fblend = 0.5 * (1.0 + TANH((LOG10(ReS)-LOG10(ReCrit))/0.2)) ! Blending function for transition between laminar and turbulent drag 
        Struts(SInd)%Cd0(EInd) = (1.0-Fblend) * Cdlam + Fblend * Cdturb ! Profile drag coefficient


    End SUBROUTINE StrutElemCoeffs

End MODULE strut
