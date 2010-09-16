SUBROUTINE pivel(NT,ntTerm,NBE,NB,NE,XP,YP,ZP,UP,VP,WP,IncGB) 

        use blade
                                                                       
        integer :: nt, ntTerm, nbe, nb, ne
        real :: XP, YP, ZP, UP, VP, WP        
        integer :: i, j, k, nei, nej, nt1                                                       
                                                                    
        ! COMPUTE THE VORTEX INDUCED VELOCITY AT POINT XP,YP,ZP              
                                                                      
        NT1=NT-1                                                                 
        UP=0.0                                                            
        VP=0.0                                                            
        WP=0.0                                                            
                                                                       
        ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO TRAILING VORTICIES  ( GT(1:NT-1,:) )   
        ! ntTerm represents the furthest away wake elements that are to be considered. (Calculated using user input xstop)          
        do i=1,ne
                do j=ntTerm,NT1
                        Call VorIVel(GT(j,i),X(j,i),Y(j,i),Z(j,i),X(j+1,i),Y(j+1,i),Z(j+1,i),XP,YP,ZP,UP,VP,WP)                              
                end do                                                       
        end do
                                                                       
                                                                       
        if (IncGB == 1) then                                                              
                ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO SPANWISE VORTICIES, including current bound vorticity ( GS(1:NT,:) )
                do i=1,nb
                        nei=(i-1)*(nbe+1)
                        do j=1,nbe
                                nej=nei+j
                                do k=ntTerm,NT
                                        Call VorIVel(GS(k,nej),X(k,nej),Y(k,nej),Z(k,nej),X(k,nej+1),Y(k,nej+1),Z(k,nej+1),XP,YP,ZP,UP,VP,WP)                      
                                end do 
                        end do                                                      
                end do 
        else
                ! CALCULATE THE VELOCITY CONTRIBUTIONS DUE TO SPANWISE VORTICIES, NOT including current bound vorticity ( GS(1:NT-1,:) )
                do i=1,nb
                        nei=(i-1)*(nbe+1)
                        do j=1,nbe
                                nej=nei+j
                                do k=ntTerm,NT1
                                        Call VorIVel(GS(k,nej),X(k,nej),Y(k,nej),Z(k,nej),X(k,nej+1),Y(k,nej+1),Z(k,nej+1),XP,YP,ZP,UP,VP,WP)                      
                                end do 
                        end do                                                      
                end do 
        end if
                                                                                                               
RETURN                                                            
END                                                               
