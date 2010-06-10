SUBROUTINE cross(ax,ay,az,bx,by,bz,cx,cy,cz) 
	
	real ax,ay,az,bx,by,bz,cx,cy,cz	

	cx = ay*bz - az*by
	cy = az*bx - ax*bz
	cz = ax*by - ay*bx

End
