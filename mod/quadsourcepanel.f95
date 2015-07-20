module quadsourcepanel

    use vecutils, only: mag3, cross3, dist_point_to_segment
    use pidef, only: pi

    implicit none

    real, parameter :: rp_far_multiplier = 10.0
    real, parameter :: ztol              = 1.0e-7
    real, parameter :: edgetol           = 1.0e-7

    private
    public quadsourcevel, quadsourcevel_full, quadsourcevel_farfield

contains

    subroutine quadsourcevel(p,p1,p2,p3,p4,sigma,selfinfluence,calcder,vel,dudx)

        ! quadsourcevel() : compute the induced velocity of a quadrilateral source panel with points p1,p2,p3,p4
        !   on a point p. The panel strength is given by sigma.

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

        real, intent(in) :: p(3), p1(3), p2(3), p3(3), p4(3), sigma
        real, intent(out) :: vel(3), dudx
        real :: pc(3), rp, lc, sz
        integer, intent(in) :: selfinfluence, calcder
        integer :: nedgesnear

        ! check if z <= ztol for all points, which is required
        if (p1(3) > ztol .or. p2(3) > ztol .or. p3(3) > ztol .or. p4(3) > ztol) then
            write(*,*) "Error: Quadrilateral source panel does not sit on a panel-local xy-plane (z must be ~ 0 for all 4 corners)!"
        end if
        
        ! compute panel center coordinates
        pc = 0.25*(p1+p2+p3+p4)

        ! compute panel characteristic length  as average length of 4 sides
        lc = 0.25*(mag3(p2-p1) + mag3(p3-p2) + mag3(p4-p3) + mag3(p1-p4))

        ! compute |r|**2 (distance from point to panel center)
        rp = mag3(p-pc)

        ! set value to zero
        nedgesnear = 0

        ! get sign of z
        sz = sign(1.0,p(3))

        if (selfinfluence == 1) then
            ! if computing induced velocity of panel at its own center
            vel = [0.0, 0.0, sz*sigma/2.0]
            
            if (calcder == 1) then
                ! add in the derivative computation
            else
                dudx = 0.0
            end if

        else
            ! near/mid/far field
            if (rp > rp_far_multiplier*lc) then
                ! far field (point solution)
                Call quadsourcevel_farfield(p,p1,p2,p3,p4,sigma,calcder,vel,dudx)

            else if (p(3) < 2*ztol) then
                ! if the z-coordinate of p is close to the surface, check for nearby edges
                if (dist_point_to_segment(p,p1,p2) <= edgetol) nedgesnear = nedgesnear + 1
                if (dist_point_to_segment(p,p2,p3) <= edgetol) nedgesnear = nedgesnear + 1
                if (dist_point_to_segment(p,p3,p4) <= edgetol) nedgesnear = nedgesnear + 1
                if (dist_point_to_segment(p,p4,p1) <= edgetol) nedgesnear = nedgesnear + 1
                
                ! handle the edge cases
                if (nedgesnear == 1) then
                    ! panel adjacent to 1 other panel, average the induced vel
                     vel = [0.0, 0.0, sz*sigma/4.0]

                else if (nedgesnear > 1) then
                    ! panel adjacent to 3 other panels, average the induced vel
                     vel = [0.0, 0.0, sz*sigma/8.0]

                else
                    ! point is close to the plane but not near an edge, just use the full solution
                    Call quadsourcevel_full(p,p1,p2,p3,p4,sigma,calcder,vel,dudx)
                end if

            else
                ! mid-field (full solution)
                Call quadsourcevel_full(p,p1,p2,p3,p4,sigma,calcder,vel,dudx)

            end if

        end if

    end subroutine quadsourcevel


    subroutine quadsourcevel_full(p,p1,p2,p3,p4,sigma,calcder,vel,dudx)
        
        ! quadsourcevel_full() : full form of the quadrilateral source panel induced velocity calculation

        real, intent(in)    :: p1(3), p2(3), p3(3), p4(3), p(3), sigma
        real, intent(out)   :: vel(3), dudx
        integer, intent(in) :: calcder

        real                :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
        real                :: xp,yp,zp
        
        real                :: d12,d23,d34,d41
        real                :: m12,m23,m34,m41
        real                :: r1,r2,r3,r4
        real                :: e1,e2,e3,e4
        real                :: h1,h2,h3,h4

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

        ! compute induced velocities
        vel(1) = sigma/(4*pi)*((y2-y1)/d12 *log((r1+r2-d12)/(r1+r2+d12)) + &
                          (y3-y2)/d23 *log((r2+r3-d23)/(r2+r3+d23)) + &
                          (y4-y3)/d34 *log((r3+r4-d34)/(r3+r4+d34)) + &
                          (y1-y4)/d41 *log((r4+r1-d41)/(r4+r1+d41)) )

        vel(2) = sigma/(4*pi)*((x1-x2)/d12 * log((r1+r2-d12)/(r1+r2+d12)) + &
                          (x2-x3)/d23 * log((r2+r3-d23)/(r2+r3+d23)) + &
                          (x3-x4)/d34 * log((r3+r4-d34)/(r3+r4+d34)) + &
                          (x4-x1)/d41 * log((r4+r1-d41)/(r4+r1+d41)) )

        vel(3) = sigma/(4*pi)*(atan((m12*e1-h1)/(zp*r1)) - atan((m12*e2-h2)/(zp*r2)) + &
                          atan((m23*e2-h2)/(zp*r2)) - atan((m23*e3-h3)/(zp*r3)) + &
                          atan((m34*e3-h3)/(zp*r3)) - atan((m34*e4-h4)/(zp*r4)) + &
                          atan((m41*e4-h4)/(zp*r4)) - atan((m41*e1-h1)/(zp*r1)))

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


end module quadsourcepanel
