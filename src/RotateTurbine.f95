SUBROUTINE RotateTurbine(delt) 
	
        use configr       
	use element
	use blade
	
	! Updates current turbine geometry
                                                                 
        ! Rotate turbine axis of rotation and origin if necessary...                                                         
                                                                 
        ! Rotate blades
        do i=1,nb
            Call RotateBlade(i,delt,RotX,RotY,RotZ,RotPX,RotPY,RotPZ)
        end do                                                        
                                                                                                          
        ! Struts etc...                                                                                                  
	               
return                                                            
end                                                               
