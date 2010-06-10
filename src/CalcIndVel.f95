subroutine CalcIndVel(NT,ntTerm,NBE,NB,NE,Point,Vel)

	use wallsoln 
	use wakeloc
	
	! Calculate wall and wake induced velocity (including bound vorticity component)
	
        real :: Point(3,1), dVel(3,1), Vel(3,1)
        integer :: nt, ntTerm, nbe, nb, ne
	
	! Calc wake induced velocity at wake locations                                                  
	CALL PIVEL(NT,ntTerm,NBE,NB,NE,Point(1,1),Point(2,1),Point(3,1),Vel(1,1),Vel(2,1),Vel(3,1),1)  

	! Calculate wall induced velocities at wake locations  
	Call WallIndVel(Point,dVel)
	Vel=Vel+dVel

return
end

