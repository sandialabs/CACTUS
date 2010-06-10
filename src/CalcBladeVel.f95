SUBROUTINE CalcBladeVel(wx,wy,wz,rx,ry,rz,uBlade,vBlade,wBlade)

	real wx,wy,wz,rx,ry,rz,uBlade,vBlade,wBlade
	
	! Blade rotation velocity (w x r) 
	CALL cross(wx,wy,wz,rx,ry,rz,uBlade,vBlade,wBlade)

Return
End