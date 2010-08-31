subroutine CalcIndVel(NT,ntTerm,NBE,NB,NE,Point,Vel)

	use wallsoln 
	use wakeloc
	
	! Calculate wall and wake induced velocity (including bound vorticity component)
	
        real :: Point(3), dVel(3), Vel(3)
        integer :: nt, ntTerm, nbe, nb, ne
	
	! Calc wake induced velocity at wake locations                                                  
	CALL PIVEL(NT,ntTerm,NBE,NB,NE,Point(1),Point(2),Point(3),Vel(1),Vel(2),Vel(3),1)  

	! Calculate wall induced velocities at wake locations  
	Call WallIndVel(Point,dVel)
	Vel=Vel+dVel

return
end

