MODULE element

    ! Blade element geometry data

    use util

    implicit none

    ! Blade geometry input structure
    ! JCM: Currently used only to store input geometry file data and blade loads outputs. Should eventually 
    ! replace the arrays below that are used for the internal calculation, similar to the strut module. This will 
    ! require a change to the blade/element/wake iterators throughout the entire code...
    type BladeType
        ! blade geometry
        integer :: NElem
        integer :: FlipN

        real, allocatable :: QCx(:)
        real, allocatable :: QCy(:)
        real, allocatable :: QCz(:)
        real, allocatable :: tx(:)
        real, allocatable :: ty(:)
        real, allocatable :: tz(:)
        real, allocatable :: CtoR(:)
        real, allocatable :: PEx(:)
        real, allocatable :: PEy(:)
        real, allocatable :: PEz(:)
        real, allocatable :: tEx(:)
        real, allocatable :: tEy(:)
        real, allocatable :: tEz(:)
        real, allocatable :: nEx(:)
        real, allocatable :: nEy(:)
        real, allocatable :: nEz(:)
        real, allocatable :: sEx(:)
        real, allocatable :: sEy(:)
        real, allocatable :: sEz(:)
        real, allocatable :: ECtoR(:)
        real, allocatable :: EAreaR(:)
        integer, allocatable :: iSect(:)

        ! Current total blade output (nomalized by machine scale parameters)
        real :: CP  ! Power coefficient due to this blade
        real :: CTR ! Torque coefficient due to this blade
        real :: CFx ! Fx coefficient due to this blade
        real :: CFy ! Fy coefficient due to this blade
        real :: CFz ! Fz coefficient due to this blade
    end type BladeType

    type(BladeType), allocatable :: Blades(:)   ! Input blade geometry

    real, allocatable :: xBE(:)     ! X location for each blade segment end (quarter chord)     
    real, allocatable :: yBE(:)     ! Y location for each blade segment end (quarter chord)     
    real, allocatable :: zBE(:)     ! Z location for each blade segment end (quarter chord) 

    real, allocatable :: txBE(:)     ! Tangential X for each blade segment end
    real, allocatable :: tyBE(:)     ! Tangential Y for each blade segment end
    real, allocatable :: tzBE(:)     ! Tangential Z for each blade segment end

    real, allocatable :: CtoR(:)     ! Chord to radius ratio for each blade segment end

    real, allocatable :: xBC(:)             ! X location for each blade segment center (quarter chord)         
    real, allocatable :: yBC(:)             ! Y location for each blade segment center (quarter chord)         
    real, allocatable :: zBC(:)             ! Z location for each blade segment center (quarter chord)  

    real :: dSGeom                          ! Geometry discretization level used in vortex core calculation
    real :: CrRef                           ! Ref chord to radius ratio     
    real, allocatable :: nxBC(:)        ! Normal X for each blade segment
    real, allocatable :: nyBC(:)        ! Normal Y for each blade segment
    real, allocatable :: nzBC(:)        ! Normal Z for each blade segment 
    real, allocatable :: txBC(:)        ! Tangential X for each blade segment   
    real, allocatable :: tyBC(:)        ! Tangential Y for each blade segment   
    real, allocatable :: tzBC(:)        ! Tangential Z for each blade segment   
    real, allocatable :: sxBC(:)        ! Spanwise X for each blade segment     
    real, allocatable :: syBC(:)        ! Spanwise Y for each blade segment     
    real, allocatable :: szBC(:)        ! Spanwise Z for each blade segment 
    ! JCM: note CircSign could be made a function of blade not element with the new geometry spec, or eliminated as it is a function of FlipN
    real, allocatable :: CircSign(:)        ! Direction of segment circulation on wake grid at positive lift
    real, allocatable :: eArea(:)       ! Element area to radius ratio for each element
    real, allocatable :: eChord(:)      ! Element chord to radius ratio for each element
    integer, allocatable :: iSect(:)        ! Array of indicies of the section table to apply to each blade element

    ! Current sum of output over all blades (nomalized by machine scale parameters)
    real :: CP_B  ! Power coefficient due to all blades
    real :: CTR_B ! Torque coefficient due to all blades
    real :: CFx_B ! Fx coefficient due to all blades
    real :: CFy_B ! Fy coefficient due to all blades
    real :: CFz_B ! Fz coefficient due to all blades

CONTAINS


    SUBROUTINE blade_geom_cns(BInd,NElem)

        implicit none

        ! Constructor for the arrays in this module

        integer :: BInd, NElem

        Blades(BInd)%NElem=NElem
        allocate(Blades(BInd)%QCx(NElem+1))    
        allocate(Blades(BInd)%QCy(NElem+1)) 
        allocate(Blades(BInd)%QCz(NElem+1)) 
        allocate(Blades(BInd)%tx(NElem+1)) 
        allocate(Blades(BInd)%ty(NElem+1)) 
        allocate(Blades(BInd)%tz(NElem+1))
        allocate(Blades(BInd)%CtoR(NElem+1)) 
        allocate(Blades(BInd)%PEx(NElem))
        allocate(Blades(BInd)%PEy(NElem))
        allocate(Blades(BInd)%PEz(NElem))
        allocate(Blades(BInd)%tEx(NElem))
        allocate(Blades(BInd)%tEy(NElem))
        allocate(Blades(BInd)%tEz(NElem))
        allocate(Blades(BInd)%nEx(NElem))
        allocate(Blades(BInd)%nEy(NElem))
        allocate(Blades(BInd)%nEz(NElem))
        allocate(Blades(BInd)%sEx(NElem))
        allocate(Blades(BInd)%sEy(NElem))
        allocate(Blades(BInd)%sEz(NElem))
        allocate(Blades(BInd)%ECtoR(NElem))
        allocate(Blades(BInd)%EAreaR(NElem))
        allocate(Blades(BInd)%iSect(NElem))         

    End SUBROUTINE blade_geom_cns


    SUBROUTINE element_cns(MaxSegEnds,MaxSegEndPerBlade)

        implicit none

        ! Constructor for the arrays in this module

        integer :: MaxSegEnds,MaxSegEndPerBlade

        allocate(xBE(MaxSegEnds))
        allocate(yBE(MaxSegEnds))
        allocate(zBE(MaxSegEnds))
        allocate(txBE(MaxSegEnds))
        allocate(tyBE(MaxSegEnds))
        allocate(tzBE(MaxSegEnds))
        allocate(CtoR(MaxSegEnds))
        allocate(xBC(MaxSegEnds))
        allocate(yBC(MaxSegEnds))
        allocate(zBC(MaxSegEnds))              
        allocate(nxBC(MaxSegEnds))      
        allocate(nyBC(MaxSegEnds))      
        allocate(nzBC(MaxSegEnds))      
        allocate(txBC(MaxSegEnds))      
        allocate(tyBC(MaxSegEnds))      
        allocate(tzBC(MaxSegEnds))      
        allocate(sxBC(MaxSegEnds))      
        allocate(syBC(MaxSegEnds))      
        allocate(szBC(MaxSegEnds))
        allocate(CircSign(MaxSegEnds))                          
        allocate(eArea(MaxSegEnds)) 
        allocate(eChord(MaxSegEnds))            
        allocate(iSect(MaxSegEnds))              

    End SUBROUTINE element_cns


    SUBROUTINE RotateBlade(BNum,delt,nrx,nry,nrz,px,py,pz)

        implicit none

        integer :: BNum
        real :: delt,nrx,nry,nrz,px,py,pz 
        integer :: nbe, j, nei, nej
        real :: vrx,vry,vrz,VMag

        ! Rotates data in blade arrays. Rotate element end geometry and recalculate element geometry.

        ! JCM: Eventually, should just be able to loop through Blades(BNum) data structure
        ! While data is still held in arrays concatenated across blades, need to replicate
        ! nbe (stored in configr) from Blades(1).NElem
        nbe=Blades(1)%NElem

        nei=1+(BNum-1)*(nbe+1)

        do j=0,nbe   
            nej=nei+j ! element index 

            ! Blade end locations (quarter chord). xBE(MaxSegEnds)
            Call QuatRot(xBE(nej),yBE(nej),zBE(nej),delt,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
            xBE(nej)=vrx                                       
            yBE(nej)=vry                                                                                              
            zBE(nej)=vrz 

            ! Tangent vectors
            Call QuatRot(txBE(nej),tyBE(nej),tzBE(nej),delt,nrx,nry,nrz,px,py,pz,vrx,vry,vrz)
            VMag=sqrt(vrx**2+vry**2+vrz**2)
            txBE(nej)=vrx/VMag
            tyBE(nej)=vry/VMag
            tzBE(nej)=vrz/VMag

        end do

        ! Calc element geometry
        Call CalcBEGeom(BNum)


    End SUBROUTINE RotateBlade


    SUBROUTINE CalcBEGeom(BNum)

        implicit none

        integer :: BNum
        integer :: nbe, nei, FlipN, nej, j
        real :: sEM, tEM, nEM
        real :: PE(3), sE(3), tE(3), normE(3), P1(3), P2(3), P3(3), P4(3), V1(3), V2(3), V3(3), V4(3), A1(3), A2(3)

        ! Calculates element geometry from element end geometry

        ! JCM: Eventually, should just be able to loop through Blades(BNum) data structure
        ! While data is still held in arrays concatenated across blades, need to replicate
        ! nbe (stored in configr) from Blades(1).NElem
        nbe=Blades(1)%NElem

        FlipN=Blades(BNum)%FlipN

        nei=1+(BNum-1)*(nbe+1)

        do j=1,nbe
            nej=nei+j

            ! Element center locations
            xBC(nej)=(xBE(nej)+xBE(nej-1))/2.0
            yBC(nej)=(yBE(nej)+yBE(nej-1))/2.0
            zBC(nej)=(zBE(nej)+zBE(nej-1))/2.0

            ! Set spanwise and tangential vectors
            sE=-(/xBE(nej)-xBE(nej-1),yBE(nej)-yBE(nej-1),zBE(nej)-zBE(nej-1)/) ! nominal element spanwise direction set opposite to QC line
            sEM=sqrt(dot_product(sE,sE))
            sE=sE/sEM
            tE=(/txBE(nej)+txBE(nej-1),tyBE(nej)+tyBE(nej-1),tzBE(nej)+tzBE(nej-1)/)/2.0
            ! Force tE normal to sE
            tE=tE-dot_product(tE,sE)*sE
            tEM=sqrt(dot_product(tE,tE))
            tE=tE/tEM
            sxBC(nej)=sE(1)
            syBC(nej)=sE(2)
            szBC(nej)=sE(3)
            txBC(nej)=tE(1)
            tyBC(nej)=tE(2)
            tzBC(nej)=tE(3)

            ! Calc normal vector
            Call cross(sE(1),sE(2),sE(3),tE(1),tE(2),tE(3),normE(1),normE(2),normE(3))
            nEM=sqrt(dot_product(normE,normE))
            normE=normE/nEM
            nxBC(nej)=normE(1)
            nyBC(nej)=normE(2)
            nzBC(nej)=normE(3)

            ! Flip normal direction if requested
            CircSign(nej)=1.0
            if (FlipN .eq. 1) then
                nxBC(nej)=-nxBC(nej)
                nyBC(nej)=-nyBC(nej)
                nzBC(nej)=-nzBC(nej)
                sxBC(nej)=-sxBC(nej)
                syBC(nej)=-syBC(nej)
                szBC(nej)=-szBC(nej)
                CircSign(nej)=-1.0
            end if

            ! Calc element area and chord
            P1=(/xBE(nej-1)-0.25*CtoR(nej-1)*txBE(nej-1),yBE(nej-1)-0.25*CtoR(nej-1)*tyBE(nej-1),zBE(nej-1)-0.25*CtoR(nej-1)*tzBE(nej-1)/)
            P2=(/xBE(nej-1)+0.75*CtoR(nej-1)*txBE(nej-1),yBE(nej-1)+0.75*CtoR(nej-1)*tyBE(nej-1),zBE(nej-1)+0.75*CtoR(nej-1)*tzBE(nej-1)/)
            P3=(/xBE(nej)+0.75*CtoR(nej)*txBE(nej),yBE(nej)+0.75*CtoR(nej)*tyBE(nej),zBE(nej)+0.75*CtoR(nej)*tzBE(nej)/)
            P4=(/xBE(nej)-0.25*CtoR(nej)*txBE(nej),yBE(nej)-0.25*CtoR(nej)*tyBE(nej),zBE(nej)-0.25*CtoR(nej)*tzBE(nej)/)
            V1=P2-P1
            V2=P3-P2
            V3=P4-P3
            V4=P1-P4
            ! Calc quad area from two triangular facets
            Call cross(V1(1),V1(2),V1(3),V2(1),V2(2),V2(3),A1(1),A1(2),A1(3))
            A1=A1/2.0
            Call cross(V3(1),V3(2),V3(3),V4(1),V4(2),V4(3),A2(1),A2(2),A2(3))
            A2=A2/2.0
            eArea(nej)=sqrt(dot_product(A1,A1))+sqrt(dot_product(A2,A2))
            ! Calc average element chord from area and span
            eChord(nej)=eArea(nej)/sEM

        end do
    End SUBROUTINE CalcBEGeom


End MODULE element
