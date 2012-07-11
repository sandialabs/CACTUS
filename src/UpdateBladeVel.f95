SUBROUTINE UpdateBladeVel(IFLG)   
        
        use configr
        use blade
        use freestream 
        use wallsoln  

        integer :: i,j,k,m,ygcErr,VFlag
        real :: Point(3), dVel(3), dUdX                                                       
                                                                 
        ! Calculate the velocity induced on the blades by wake, wall, and freestream
                                                                       
        J=NT                                                              
        NT1=NT-1  

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
                        CALL CalcFreestream(Y(J,I),UFS(J,I),VFS(J,I),WFS(J,I),ygcErr)                                
                                                       
                        USUM=0.0                                                          
                        VSUM=0.0                                                          
                        WSUM=0.0                                                                                                                                                         
                        if (NT > 1) then                                         
                                                                       
                                ! Calculate wake velocity at blade elements (excluding bound vorticity component)
                                Call BladeIndVel(NT,ntTerm,NBE,NB,NE,X(J,I),Y(J,I),Z(J,I),USUM,VSUM,WSUM,dUdX,2,0)                                                  
  
                                ! Calculate wall induced velocities at blade locations   
                                Point=[X(J,I),Y(J,I),Z(J,I)]
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
                Call BladeIndVel(NT,ntTerm,NBE,NB,NE,X(J,I),Y(J,I),Z(J,I),UP,VP,WP,dUdX,1,0) 
                
                ! Set wake and wall velocities on blade
                U(J,I)=USUM+UP                                                       
                V(J,I)=VSUM+VP                                                       
                W(J,I)=WSUM+WP 
                                                                      
        end do   
                                                       
RETURN                                                            
END                                                               
