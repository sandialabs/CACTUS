SUBROUTINE SGeomSetup() 

	use parameters

    use configr       
	use element
    use strut 
    use airfoil            
	use pidef      

    implicit none

    integer :: i, j, NElem, BIndS, EIndS, BIndE, EIndE, iSc
    real :: sx, sy, sz, VMag, LR

	! Sets up strut geometry arrays.        

	! Set initial geometry (zero theta)                                                            
	do i=1,NStrut  
        NElem=Struts(i)%NElem 
        BIndS=Struts(i)%BIndS
        EIndS=Struts(i)%EIndS
        BIndE=Struts(i)%BIndE
        EIndE=Struts(i)%EIndE 

        ! Get blade thickness to chord
        if (BIndS > 0) then
            iSc=Blades(BIndS)%iSect(EIndS)
            Struts(i)%tcS=tc(iSc) 
        else 
            Struts(i)%tcS=0.0
        end if
        if (BIndE > 0) then
            iSc=Blades(BIndE)%iSect(EIndE)
            Struts(i)%tcE=tc(iSc) 
        else 
            Struts(i)%tcE=0.0
        end if

        ! Init length sum
        LR=0.0
        do j=1,NElem 
            VMag=sqrt(Struts(i)%sEx(j)**2+Struts(i)%sEy(j)**2+Struts(i)%sEz(j)**2)
            ! Sum length
            LR=LR+VMag
        end do
        Struts(i)%LR=LR

        ! Calc element geometry
        Call CalcSEGeom(i)

    end do

    return                                                            
end SUBROUTINE SGeomSetup
