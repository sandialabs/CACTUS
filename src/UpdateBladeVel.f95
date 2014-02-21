SUBROUTINE UpdateBladeVel(IFLG)   

    use configr
    use blade
    use wallsoln  

    integer :: i,ygcErr
    real :: Point(3), dVel(3), dUdX                                                       

    ! Calculate the velocity induced on the blades by wake, wall, and freestream


    if (iflg .eq. 0) then
        ! re-initialize uiwake viwake wiwake as we are beginning a new time step
        uiwake(:)=0.0
        viwake(:)=0.0
        wiwake(:)=0.0
    end if


    do I=1,NE                                                      

        ! If flag is set, just recompute the velocity contiribution due to the bound vorticies on the blades.
        ! Otherwise, calculate all wake, wall and freestream induced velocity.                                                                              
        if (IFLG .eq. 0) then 

            !  Calculate freestream velocity at blade elements
            CALL CalcFreestream(X(NT,I),Y(NT,I),Z(NT,I),UFSB(I),VFSB(I),WFSB(I),ygcErr)

            ! Set freestream velocity of next shed wake elements to that calculated on the blade
            UFS(NT,I)=UFSB(I)
            VFS(NT,I)=VFSB(I)
            WFS(NT,I)=WFSB(I)


            USUM=0.0                                                          
            VSUM=0.0                                                          
            WSUM=0.0                                                                                                                                                         
            if (NT > 1) then                                         

                ! Calculate wake velocity at blade elements (excluding bound vorticity component)
                Call BladeIndVel(NT,ntTerm,NBE,NB,NE,X(NT,I),Y(NT,I),Z(NT,I),USUM,VSUM,WSUM,dUdX,2,0)

                ! Calculate wall induced velocities at blade locations   
                Point=[X(NT,I),Y(NT,I),Z(NT,I)]
                Call WallIndVel(Point,dVel)
                USUM=USUM+dVel(1)
                VSUM=VSUM+dVel(2)
                WSUM=WSUM+dVel(3)

            end if

            uiwake(I)=USUM                                                     
            viwake(I)=VSUM                                                     
            wiwake(I)=WSUM                                                     
        else                                                         
            USUM=uiwake(I)                                                     
            VSUM=viwake(I)                                                     
            WSUM=wiwake(I)                                                     
        end if

        ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO JUST THE BOUND VORTICIES ON THE BLADES ( GS(NT,:) ) 
        Call BladeIndVel(NT,ntTerm,NBE,NB,NE,X(NT,I),Y(NT,I),Z(NT,I),UP,VP,WP,dUdX,1,0)

        ! Set wake and wall velocities on blade
        UB(I)=USUM+UP
        VB(I)=VSUM+VP
        WB(I)=WSUM+WP

        ! Set induced velocity of next shed wake elements
        if (iut .eq. -2) then
            ! Fix wake velocities at freestream velocity
            U(NT,I)=0.0
            V(NT,I)=0.0
            W(NT,I)=0.0
        else
            U(NT,I)=UB(I)
            V(NT,I)=VB(I)
            W(NT,I)=WB(I)
        end if

    end do

    RETURN                                                            
END SUBROUTINE UpdateBladeVel
