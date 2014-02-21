SUBROUTINE BladeIndVel(NT,ntTerm,NBE,NB,NE,XP,YP,ZP,UP,VP,WP,DUDX,Mode,CalcDer)

    use blade

    integer :: nt, ntTerm, nbe, nb, ne, Mode, CalcDer
    real :: XP, YP, ZP, UP, VP, WP, DUDX        
    integer :: i, j, k, nei, nej, nt1 
    integer :: VFlag                                                      

    ! COMPUTE THE VORTEX INDUCED VELOCITY AT POINT XP,YP,ZP FROM THE BLADE SYSTEM  
    ! Mode: 0 -> Calc whole system
    !       1 -> Calc bound vorticity only
    !       2 -> Calc whole system minus bound vorticity 
    ! CalcDer: 1 to calc DUDX

    NT1=NT-1                                                                 
    UP=0.0                                                            
    VP=0.0                                                            
    WP=0.0     
    DUDX=0.0        

    if (Mode == 0) then
        ! Calc induced velocity from the whole blade system (bound, spanwise, and trailing wake vorticity)                                                            

        ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO TRAILING VORTICIES  ( GT(1:NT-1,:) )   
        ! ntTerm represents the furthest away wake elements that are to be considered. (Calculated using user input xstop) 
        VFlag=1         
        do i=1,ne
            do j=ntTerm,NT1
                Call VorIVel(VFlag,CalcDer,GT(j,i),X(j,i),Y(j,i),Z(j,i),X(j+1,i),Y(j+1,i),Z(j+1,i),XP,YP,ZP,UP,VP,WP,DUDX)                              
            end do
        end do

        ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO SPANWISE VORTICIES, including current bound vorticity ( GS(1:NT,:) )
        do i=1,nb
            nei=(i-1)*(nbe+1)
            do j=1,nbe
                nej=nei+j
                do k=ntTerm,NT
                    if (k==NT) then
                        VFlag=0
                    else
                        VFlag=2
                    end if
                    Call VorIVel(VFlag,CalcDer,GS(k,nej),X(k,nej),Y(k,nej),Z(k,nej),X(k,nej+1),Y(k,nej+1),Z(k,nej+1),XP,YP,ZP,UP,VP,WP,DUDX)                      
                end do
            end do
        end do

    else if (Mode == 1) then
        ! Calc induced velocity from bound vorticity only                                                                               

        VFlag=0                                                 
        do i=1,nb
            nei=(i-1)*(nbe+1)
            do j=1,nbe
                nej=nei+j
                k=NT
                Call VorIVel(VFlag,CalcDer,GS(k,nej),X(k,nej),Y(k,nej),Z(k,nej),X(k,nej+1),Y(k,nej+1),Z(k,nej+1),XP,YP,ZP,UP,VP,WP,DUDX)                      
            end do
        end do

    else
        ! Calc induced velocity from spanwise and trailing wake vorticity (exclude bound vorticity)

        ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO TRAILING VORTICIES  ( GT(1:NT-1,:) )   
        ! ntTerm represents the furthest away wake elements that are to be considered. (Calculated using user input xstop) 
        VFlag=1         
        do i=1,ne
            do j=ntTerm,NT1
                Call VorIVel(VFlag,CalcDer,GT(j,i),X(j,i),Y(j,i),Z(j,i),X(j+1,i),Y(j+1,i),Z(j+1,i),XP,YP,ZP,UP,VP,WP,DUDX)                              
            end do
        end do

        ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO SPANWISE VORTICIES, NOT including current bound vorticity ( GS(1:NT-1,:) )
        VFlag=2
        do i=1,nb
            nei=(i-1)*(nbe+1)
            do j=1,nbe
                nej=nei+j
                do k=ntTerm,NT1
                    Call VorIVel(VFlag,CalcDer,GS(k,nej),X(k,nej),Y(k,nej),Z(k,nej),X(k,nej+1),Y(k,nej+1),Z(k,nej+1),XP,YP,ZP,UP,VP,WP,DUDX)                      
                end do
            end do
        end do

    end if

    RETURN                                                            
END SUBROUTINE BladeIndVel
