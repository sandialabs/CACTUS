module geomutils

    ! module geomutils : functions for querying 2-D and 3-D geometries
    
    implicit none

    private
    public dist_point_to_segment_2d, is_point_in_quad
    
contains

     function dist_point_to_segment_2d(x,y,x1,y1,x2,y2)

        ! dist_point_to_segment_2d() : returns the shortest distance from a point p to a line segment.
        !   Based on the method of http://paulbourke.net/geometry/pointlineplane/
        !
        !   Inputs:
        !   =======
        !   x,y         : the coordinates of the point
        !
        !   Output:
        !   ======= 
        !   x1,y1,x2,y2 : the start and end points of the line segment
        !

        real, intent(in)  :: x,y,x1,y1,x2,y2
        real :: dist_point_to_segment_2d

        real :: u, dx, dy
        real :: x_closest, y_closest

        ! check if points are coincident
        if (x1==x2 .and. y1==y2) then
            write(*,*) 'Error in dist_point_to_segment_2d : Points are coincident!'
            stop
        end if
        
        dx = x2-x1
        dy = y2-y1

        u = ((x-x1)*dx + (y-y1)*dy) / abs(dx**2 - dy**2)

        if (u < 0) then
            x_closest = x1
            y_closest = y1
        else if (u > 1) then
            x_closest = x2
            y_closest = y2
        else
            x_closest = x1+u*dx
            y_closest = y1+u*dy
        end if

        dist_point_to_segment_2d = ((x_closest-x)**2 + (y_closest-y)**2)**0.5

    end function dist_point_to_segment_2d


    function is_point_in_quad(p_list, p)

        ! is_point_in_quad() : For a quadrilateral whose corners are given in p_list, returns
        !   True if the point p lies within the polygon and False if the point p lies outside of it.
        !   Based on the method of http://paulbourke.net/geometry/polygonmesh/
        !
        !   Inputs:
        !   =======
        !   p_list                     : a (4,3) array of points. expects that p(:,3) = 0 (the z-coordinate is zero)
        !
        !   Outputs:
        !   ========
        !   is_point_in_quad (logical) : .True. if point is inside the quadrilateral, .False. otherwise

        real, intent(in)   :: p_list(4,3)
        real, intent(in)   :: p(3)
        logical            :: is_point_in_quad
        
        integer, parameter :: nsides = 4
        real               :: p1(3), p2(3)
        real               :: xinters

        integer            :: i, counter

        counter = 0

        p1 = p_list(1,1:3)

        do i=1,nsides
            p2 = p_list(mod(i,nsides)+1,1:3)
            if (p(2) > min(p1(2),p2(2))) then
                if (p(2) <= max(p1(2),p2(2))) then
                    if (p(1) <= max(p1(1),p2(1))) then
                        if (p1(2) /= p2(2)) then
                            xinters = (p(2)-p1(2))*(p2(1)-p1(1))/(p2(2)-p1(2)) + p1(1)
                            if (p1(1) == p2(1) .or. p(1) <= xinters) then
                                counter = counter + 1
                            end if  
                        end if
                    end if
                end if
            end if
            p1 = p2
        end do

        if (mod(counter,2) == 0) then
            is_point_in_quad = .False.
        else
            is_point_in_quad = .True.
        end if

    end function is_point_in_quad


end module geomutils