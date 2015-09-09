SUBROUTINE BGeomSetup()

    use parameters

    use configr
    use element
    use pidef
    use util

    implicit none

    integer i, j, nei, nej
    real dx, dy, dz, vMag, ds

    ! Sets up blade geometry arrays for each blade, starting from the root of the first blade, and continuing for each blade.

    ! Set reference cr
    CrRef=0.0
    do i=1,nb
        do j=1,nbe
            CrRef=max(CrRef,Blades(i)%ECtoR(j))
        end do
    end do

    ! Set representative geometry discretization level (for vortex core calculation)
    dSGeom=0.0
    do i=1,nb
        do j=1,nbe
            dx=Blades(i)%QCx(j+1)-Blades(i)%QCx(j)
            dy=Blades(i)%QCy(j+1)-Blades(i)%QCy(j)
            dz=Blades(i)%QCz(j+1)-Blades(i)%QCz(j)
            ds=sqrt(dx**2+dy**2+dz**2)
            dsGeom=max(dsGeom,ds)
        end do
    end do

    ! Set initial turbine geometry (zero theta)
    do i=1,nb

        nei=1+(i-1)*(nbe+1)

        ! Set blade end geometry
        do j=0,nbe
            nej=nei+j ! element index

            xBE(nej)=Blades(i)%QCx(j+1)
            yBE(nej)=Blades(i)%QCy(j+1)
            zBE(nej)=Blades(i)%QCz(j+1)

            txBE(nej)=Blades(i)%tx(j+1)
            tyBE(nej)=Blades(i)%ty(j+1)
            tzBE(nej)=Blades(i)%tz(j+1)

            CtoR(nej)=Blades(i)%CtoR(j+1)
        end do


        ! Calc element geometry

        ! JCM: currently, although these values are for each element, they are held in arrays sized for element ends, where the first value for each blade
        ! is simply ignored in bsload (where these values are used).

        ! JCM: These zeros are ignored...
        eArea(nei)=0.0
        eChord(nei)=0.0
        iSect(nei)=0.0
        xBC(nei)=0.0
        yBC(nei)=0.0
        zBC(nei)=0.0
        nxBC(nei)=0.0
        nyBC(nei)=0.0
        nzBC(nei)=0.0
        txBC(nei)=0.0
        tyBC(nei)=0.0
        tzBC(nei)=0.0
        sxBC(nei)=0.0
        syBC(nei)=0.0
        szBC(nei)=0.0
        Call CalcBEGeom(i)

        do j=1,nbe
            nej=nei+j
            iSect(nej)=Blades(i)%iSect(j)
        end do

    end do

    return
end SUBROUTINE BGeomSetup
