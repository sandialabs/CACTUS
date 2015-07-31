module quadsourcepanel

    use vecutils, only: mag3, cross3
    use geomutils, only: dist_point_to_segment_2d, is_point_in_quad
    use pidef, only: pi

    implicit none

    real, parameter :: rp_far_multiplier = 10.0    ! multiplier for far-field threshold (far-field if rp > rp_far_multplier*lc)
    real, parameter :: ztol              = 1.0e-7  ! threshold for distance from plane
    real, parameter :: edgetol           = 1.0e-7  ! threshold for distance from edge
    real, parameter :: eps               = 1.0e-7  ! small value to avoid taking log(0) in quadsourcevel_full
    
    private quadsourcevel_full, quadsourcevel_farfield, reverse_point_order
    public quadsourcevel

contains

    subroutine quadsourcevel(p,p1,p2,p3,p4,sigma,selfinfluence,calcder,vel,dudx,info)

        ! quadsourcevel() : compute the induced velocity of a quadrilateral source panel with points p1,p2,p3,p4
        !   on a point p. The panel strength is given by sigma.
        !   Requires that the panel be on the x-y plane (z=0).
        !   Expects that p1,p2,p3,p4 are given CW when the plane is viewed from behind (in the direction of the normal)
        !
        !   Inputs:
        !   =======
        !   p1,p2,p3,p4     : locations of panel nodes
        !   p               : location at which to compute induced velocity
        !   selfinfluence   : 1 to compute the velocity at the panel's center
        !   calcder         : 1 to compute the velocity dudx
        !
        !   Outputs:
        !   ========
        !   vel             : a vector of u,v,w velocities
        !   dudx            : derivative du/dx
        !   info            : integer containing information about the solution method
        !                     0 - self influence
        !                     1 - point is on plane near single edge
        !                     2 - point is on plane near two edges
        !                     3 - point is on plane outside quad (uv computed, w=0)
        !                     4 - point is on plane inside quad  (uv computed, w=sigma/2)
        !                     5 - mid-field (full solution)
        !                     6 - far-field (point approximation)
        !                     

        real, intent(in)     :: p(3), p1(3), p2(3), p3(3), p4(3), sigma
        integer, intent(in)  :: selfinfluence, calcder
        real, intent(out)    :: vel(3), dudx
        integer, intent(out) :: info
        real                 :: p1r(3), p2r(3), p3r(3), p4r(3)
        real                 :: pc(3), rp, lc, sz
        integer              :: nedgesnear
        real                 :: quad_points(4,3)

        ! reverse the points!
        !! in CACTUS the node points are given CW, while the equations expect CCW ordering!!
        call reverse_point_order(p1,p2,p3,p4, &
                                 p1r,p2r,p3r,p4r)

        ! check if z <= ztol for all points, which is required
        if (abs(p1r(3)) > ztol .or. abs(p2r(3)) > ztol .or. abs(p3r(3)) > ztol .or. abs(p4r(3)) > ztol) then
            write(*,*) "Error: Quadrilateral source panel does not sit on a panel-local xy-plane (z must be ~ 0 for all 4 corners)!"
            stop
        end if
        
        ! compute panel center coordinates
        pc = 0.25*(p1r+p2r+p3r+p4r)

        ! compute panel characteristic length  as average length of 4 sides
        lc = 0.25*(mag3(p2r-p1r) + mag3(p3r-p2r) + mag3(p4r-p3r) + mag3(p1r-p4r))

        ! compute |r| (distance from point to panel center)
        rp = mag3(p-pc)

        ! set value to zero
        nedgesnear = 0

        ! get sign of z
        sz = sign(1.0,p(3))

        if (selfinfluence == 1) then
            ! if computing induced velocity of panel at its own center
            vel = [0.0, 0.0, sz*sigma/2.0]
            
            if (calcder == 1) then
                !! todo: add in the derivative computation
                ! dudx = 
            else
                dudx = 0.0
            end if

            ! set info flag
            info = 0

        else
            ! near/mid/far field
            
            if (rp > rp_far_multiplier*lc) then
                ! far field (point solution)
                Call quadsourcevel_farfield(p,p1r,p2r,p3r,p4r,sigma,calcder,vel,dudx)

                ! set info flag
                info = 6

            else
                ! point is either in the near-field or on the panel

                if (abs(p(3)) <= ztol) then
                    ! if point is on the panel
                    Call quadsourcevel_full(p,p1r,p2r,p3r,p4r,sigma,calcder,vel,dudx)

                    ! finally, check if point is near edges
                    if (dist_point_to_segment_2d(p(1),p(2),p1r(1),p1r(2),p2r(1),p2r(2)) <= edgetol) nedgesnear = nedgesnear + 1
                    if (dist_point_to_segment_2d(p(1),p(2),p2r(1),p2r(2),p3r(1),p3r(2)) <= edgetol) nedgesnear = nedgesnear + 1
                    if (dist_point_to_segment_2d(p(1),p(2),p3r(1),p3r(2),p4r(1),p4r(2)) <= edgetol) nedgesnear = nedgesnear + 1
                    if (dist_point_to_segment_2d(p(1),p(2),p4r(1),p4r(2),p1r(1),p1r(2)) <= edgetol) nedgesnear = nedgesnear + 1

                    if (nedgesnear == 1) then
                        ! panel adjacent to 1 other panel, average the induced w-velocity with 2 other panels
                        vel(3) = 0.25*sigma

                        ! set info flag
                        info = 1

                    else if (nedgesnear == 2) then
                        ! point is 
                        ! panel adjacent to 3 other panels, average the induced w-velocity with 4 panels
                        vel(3) = 0.125*sigma

                        ! set info flag
                        info = 2

                    else if (nedgesnear > 2) then
                        ! point is near more than 2 edges
                        write(*,*) "Warning in quadsourcevel(): Point is near more than two panel edges"
                        write(*,*) " -- infinitesimally small panel?"

                        ! panel adjacent to 3+? other panels, average the induced w-velocity with 4 panels
                        vel(3) = 0.125*sigma

                        ! set info flag
                        info = 2

                    else
                        ! point is not near any edges

                        ! if the point isn't near any edges, set the w-velocity based on whether or not the point
                        ! is in inside/outside the quadrilateral

                        ! check if the point is in the quad
                        quad_points(1,:) = p1r(1:3)
                        quad_points(2,:) = p2r(1:3)
                        quad_points(3,:) = p3r(1:3)
                        quad_points(4,:) = p4r(1:3)

                        if (is_point_in_quad(quad_points, p)) then
                            ! if the point is inside the polygon, the w-velocity is 0.5*sigma
                            vel(3) = sz*sigma/2.0

                            ! set info flag
                            info = 4
                            
                        else
                            ! if the point is outside the polygon, the w-velocity is 0.0
                            vel(3) = 0.0

                            ! set info flag
                            info = 3

                        end if
                    end if

                else
                    ! point is neither far away nor near the panel -- use the full-field solution
                    Call quadsourcevel_full(p,p1r,p2r,p3r,p4r,sigma,calcder,vel,dudx)

                    ! set info flag
                    info = 5

                end if ! near/mid test
            end if ! far-field test
        end if ! selfinfluence test

    end subroutine quadsourcevel


    subroutine quadsourcevel_full(p,p1,p2,p3,p4,sigma,calcder,vel,dudx)
        
        ! quadsourcevel_full() : full form of the quadrilateral source panel induced velocity calculation

        real, intent(in)    :: p1(3), p2(3), p3(3), p4(3), p(3), sigma
        real, intent(out)   :: vel(3), dudx
        integer, intent(in) :: calcder

        real                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        real                :: xp,yp,zp
        
        ! intermediate values
        real                :: d12,d23,d34,d41
        real                :: m12,m23,m34,m41
        real                :: r1,r2,r3,r4
        real                :: e1,e2,e3,e4
        real                :: h1,h2,h3,h4

        ! repeated calculations (for speed savings)
        real                :: a1,a2,a3,a4
        real                :: w1,w2,w3,w4,w5,w6,w7,w8

        ! breakout components of nodes
        x1 = p1(1)
        y1 = p1(2)
        z1 = p1(3)

        x2 = p2(1)
        y2 = p2(2)
        z2 = p2(3)

        x3 = p3(1)
        y3 = p3(2)
        z3 = p3(3)

        x4 = p4(1)
        y4 = p4(2)
        z4 = p4(3)

        ! breakout components of point
        xp = p(1)
        yp = p(2)
        zp = p(3)

        ! compute intermediate values (panel geometry dependent)
        d12 = ((x2-x1)**2 + (y2-y1)**2)**0.5
        d23 = ((x3-x2)**2 + (y3-y2)**2)**0.5
        d34 = ((x4-x3)**2 + (y4-y3)**2)**0.5
        d41 = ((x1-x4)**2 + (y1-y4)**2)**0.5

        m12 = (y2-y1) / (x2-x1)
        m23 = (y3-y2) / (x3-x2)
        m34 = (y4-y3) / (x4-x3)
        m41 = (y1-y4) / (x1-x4)

        ! compute intermediate values (panel geometry and point dependent)
        r1 = ((xp-x1)**2 + (yp-y1)**2 + zp**2)**0.5
        r2 = ((xp-x2)**2 + (yp-y2)**2 + zp**2)**0.5
        r3 = ((xp-x3)**2 + (yp-y3)**2 + zp**2)**0.5
        r4 = ((xp-x4)**2 + (yp-y4)**2 + zp**2)**0.5 

        e1 = (xp-x1)**2 + zp**2
        e2 = (xp-x2)**2 + zp**2
        e3 = (xp-x3)**2 + zp**2
        e4 = (xp-x4)**2 + zp**2

        h1 = (xp-x1)*(yp-y1)
        h2 = (xp-x2)*(yp-y2)
        h3 = (xp-x3)*(yp-y3)
        h4 = (xp-x4)*(yp-y4)

        ! compute repeated values (for u,v computations)
        a1 = (r1+r2-d12)/(r1+r2+d12)
        a2 = (r2+r3-d23)/(r2+r3+d23)
        a3 = (r3+r4-d34)/(r3+r4+d34)
        a4 = (r4+r1-d41)/(r4+r1+d41)

        ! compute temporary values (for w computations)
        w1 = (m12*e1-h1)/(zp*r1)
        w2 = (m12*e2-h2)/(zp*r2)
        w3 = (m23*e2-h2)/(zp*r2)
        w4 = (m23*e3-h3)/(zp*r3)
        w5 = (m34*e3-h3)/(zp*r3)
        w6 = (m34*e4-h4)/(zp*r4)
        w7 = (m41*e4-h4)/(zp*r4)
        w8 = (m41*e1-h1)/(zp*r1)

        ! check for bad values of a,b,c,d, which can lead to 
        if (a1 <= 0) a1 = eps
        if (a2 <= 0) a2 = eps
        if (a3 <= 0) a3 = eps
        if (a4 <= 0) a4 = eps
        
        ! compute induced velocities
        vel(1) = sigma/(4*pi)*((y2-y1)/d12 *log(a1) + &
                               (y3-y2)/d23 *log(a2) + &
                               (y4-y3)/d34 *log(a3) + &
                               (y1-y4)/d41 *log(a4) )

        vel(2) = sigma/(4*pi)*((x1-x2)/d12 * log(a1) + &
                               (x2-x3)/d23 * log(a2) + &
                               (x3-x4)/d34 * log(a3) + &
                               (x4-x1)/d41 * log(a4) )

        vel(3) = sigma/(4*pi)*(atan(w1) - atan(w2) + &
                               atan(w3) - atan(w4) + &
                               atan(w5) - atan(w6) + &
                               atan(w7) - atan(w8) )

        if (calcder == 1) then
            !! need to add derivative!
            dudx = 0.0
        end if

    end subroutine quadsourcevel_full


    subroutine quadsourcevel_farfield(p,p1,p2,p3,p4,sigma,calcder,vel,dudx)

        ! quadsourcevel_farfield() : far-field (point-source) quadrilateral source panel induced velocity calculation

        real, intent(in)    :: p(3), p1(3), p2(3), p3(3), p4(3), sigma
        integer, intent(in) :: calcder
        real, intent(out)   :: vel(3), dudx

        real :: pc(3), a(3), b(3), c(3), d(3)
        real :: area
        real :: rp_squared
        real :: xp,yp,zp,xc,yc,zc

        ! compute panel center coordinates
        pc = 0.25*(p1+p2+p3+p4)

        ! breakout
        xp = p(1)
        yp = p(2)
        zp = p(3)

        xc = pc(1)
        yc = pc(2)
        zc = pc(3)

        ! compute |r|**2 (distance from point to panel center)
        rp_squared = sum((p-pc)**2)

        ! compute vectors of sides
        a = p2-p1
        b = p3-p2
        c = p4-p3
        d = p1-p4

        ! compute the quadrilateral area using cross product
        area = 0.5*mag3(cross3(b+c,a+b))

        ! compute induced velocities
        vel(1) = sigma/(4*pi) * area * (xp-xc)/rp_squared**(1.5)
        vel(2) = sigma/(4*pi) * area * (yp-yc)/rp_squared**(1.5)
        vel(3) = sigma/(4*pi) * area * (zp-zc)/rp_squared**(1.5)

        ! derivative w.r.t. x
        if (calcder == 1) then
            dudx = 0.0
        end if

    end subroutine quadsourcevel_farfield


    subroutine reverse_point_order(p1,p2,p3,p4,p1r,p2r,p3r,p4r)
        
        ! reverse_points() : reverses the ordering of four points given by p1,p2,p3,p4

        real, intent(in) :: p1(3), p2(3), p3(3), p4(3)
        real, intent(out) :: p1r(3), p2r(3), p3r(3), p4r(3)
        real :: temp1(3), temp2(3)

        temp1 = p1
        temp2 = p2
        
        p1r = p4
        p2r = p3
        p3r = temp2
        p4r = temp1

    end subroutine reverse_point_order


end module quadsourcepanel