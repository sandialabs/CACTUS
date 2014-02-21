SUBROUTINE InputGeom(FN)

    use element
    use strut
    use varscale 
    use configr    

    implicit none

    integer, parameter :: MaxReadLine = 1000    
    character(MaxReadLine) :: FN    ! path to geometry input file 

    integer :: NElem, i
    character(MaxReadLine) :: ReadLine

    ! Read geometry data file

	! Format example:
	!	NBlade: 3
	!	NStrut: 3
	!	RotN: 0 0 1
	!	RotP: 0 0 0
	!	RefAR: 2.0
	!	RefR: 10.0
	!	Type: VAWT
	!	Blade 1:
	!	    NElem: 5
	!	    FlipN: 0
	!	    QCx: 0 0 0 0 0 0
	!	    QCy: 1 2 3 4 5 6
	!	    QCz: 1 1 1 1 1 1
	!	    tx: 1 1 1 1 1 1
	!	    ty: 0 0 0 0 0 0
	!	    tz: 0 0 0 0 0 0
	!	    CtoR: .1 .1 .1 .1 .1 .1
	!	    PEx: 0 0 0 0 0
	!	    PEy: 1 2 3 4 5
	!	    PEz: 1 1 1 1 1
	!	    tEx: 0 0 0 0 0
	!	    tEy: 1 2 3 4 5
	!	    tEz: 1 1 1 1 1
	!	    nEx: 0 0 0 0 0
	!	    nEy: 1 2 3 4 5
	!	    nEz: 1 1 1 1 1
	!	    sEx: 1 1 1 1 1
	!	    sEy: 0 0 0 0 0
	!	    sEz: 1 2 3 4 5
	!	    ECtoR: .1 .1 .1 .1 .1
	!	    EAreaR: .1 .1 .1 .1 .1
	!	    iSect: 1 1 1 1 1
	!	Blade 2:
	!	    ...
	!	Strut 1:
	!	    NElem: 5
	!	    TtoC: .15
	!	    MCx: 0 0 0 0 0 0
	!	    MCy: 1 2 3 4 5 6
	!	    MCz: 1 1 1 1 1 1
	!	    CtoR: .1 .1 .1 .1 .1 .1
	!	    PEx: 0 0 0 0 0
	!	    PEy: 1 2 3 4 5
	!	    PEz: 1 1 1 1 1
    !       sEx: 0 0 0 0 0
    !       sEy: 1 2 3 4 5
    !       sEz: 1 1 1 1 1
	!	    ECtoR: .1 .1 .1 .1 .1
	!	    EAreaR: .1 .1 .1 .1 .1
	!	    BIndS: 0
	!	    EIndS: 0
	!	    BIndE: 1
	!	    EIndE: 3
	!	Strut 2:
	!	    ...


    ! Open input file for this section
    open(15, file=FN)

    ! Read header data
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) nb

    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) NStrut

    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) RotX, RotY, RotZ

    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) RotPX, RotPY, RotPZ

    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) at

    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Rmax

    read(15,'(A)') ReadLine

    ! Allocate blades struct
    allocate(Blades(nb))

    ! Allocate struts struct
    if (NStrut>0) then
        allocate(Struts(NStrut))     
    end if

    ! Read blade data
    do i=1,nb

        read(15,'(A)') ReadLine

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) NElem

        ! Allocate arrays in blade structure
        Call blade_geom_cns(i,NElem)

        ! JCM: Currently (until blade/wake geometry module is restructured), NElem must be the same for all blades.
        nbe=NElem

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%FlipN

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%QCx(1:NElem+1)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%QCy(1:NElem+1)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%QCz(1:NElem+1)  

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tx(1:NElem+1)                      

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%ty(1:NElem+1) 

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tz(1:NElem+1) 

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%CtoR(1:NElem+1) 

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%PEx(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%PEy(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%PEz(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tEx(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tEy(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tEz(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%nEx(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%nEy(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%nEz(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%sEx(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%sEy(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%sEz(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%ECtoR(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%EAreaR(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%iSect(1:NElem) 

    end do

    ! Read strut data
    do i=1,NStrut

        read(15,'(A)') ReadLine

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) NElem

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%TtoC

        ! Allocate arrays in blade structure
        Call strut_comp_cns(i,NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%MCx(1:NElem+1)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%MCy(1:NElem+1)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%MCz(1:NElem+1)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%CtoR(1:NElem+1)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%PEx(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%PEy(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%PEz(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%sEx(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%sEy(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%sEz(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%ECtoR(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%EAreaR(1:NElem)

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%BIndS

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%EIndS

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%BIndE

        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Struts(i)%EIndE

    end do

    ! Close input file 
    close(15)

    Return
End SUBROUTINE InputGeom
