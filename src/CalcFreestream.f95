SUBROUTINE CalcFreestream(yElem,u,v,w,ygcErr)

	use shear
	
	real u, v, w
	real yElem
	integer ygcErr
	
	! Freestream velocity (with ground shear model)    
	! At y/R = 0, u/Uinf = 0. At y/R = yref, u/Uinf = 1
	! slex = 0 : Constant freestream
	! slex = 1/2 : Laminar shear layer (approx)
	! slex = 1/7 : Turbulent shear layer (approx)   
	if ((yElem+ygc) <= 0.0) then                                        
		ygcErr=1                                                           
		u=0.0
		v=0.0
		w=0.0   
	else
		u=((yElem+ygc)/yref)**slex 
		v=0.0
		w=0.0                                                        
	end if                         

Return
End