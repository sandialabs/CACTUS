module vecutils

    implicit none
    
    private
    public cross3, mag3, dist_point_to_segment

contains

    pure function cross3(a,b)

        ! cross3() : returns the cross product a x b, with a(3) and b(3)

        real, intent(in) :: a(3), b(3)
        real :: cross3(3)
        
        cross3(1) = a(2) * b(3) - a(3) * b(2)
        cross3(2) = a(3) * b(1) - a(1) * b(3)
        cross3(3) = a(1) * b(2) - a(2) * b(1)

    end function cross3


    pure function mag3(a)

        ! mag3() : returns the magnitude of a(3)

        real, intent(in) :: a(3)
        real :: mag3

        mag3 = sqrt(sum(a**2))

    end function mag3


    pure function dist_point_to_segment(p,p0,p1)

        ! dist_point_to_segment() : returns the shortest distance from a point p to a line 
        !   segment with endpoints p0, p1
        !   based on the method http://geomalgorithms.com/a02-_lines.html

        real, intent(in) :: p(3),p0(3),p1(3)
        real :: dist_point_to_segment
        real :: v(3), w(3), c1, c2, b, pb(3)

        v = p1 - p0
        w = p - p0

        c1 = dot_product(w,v)
        if (c1 <= 0) then
            dist_point_to_segment = mag3(p-p0)
            return
        end if
            
        c2 = dot_product(v,v)
        if (c2 <= 0) then
            dist_point_to_segment = mag3(p-p1)
            return
        end if

        b = c1/c2
        pb = p0 + b*v
        dist_point_to_segment = mag3(p-pb)

    end function dist_point_to_segment


end module vecutils