SUBROUTINE bgeom(nGeom) 
	
        use configr       
	use element
	use blade
	
	! JCM: Gets blade element geometry for each blade starting from the bottom of the first blade                                                         
	                                                                                                    
	do i=1,nb                                                      
		                                            
		nei=1+(i-1)*(nbe+1)                                               
	                                                                              
		! Blade element end locations (quarter chord). 
		do j=0,nbe   
			nej=nei+j ! element index
			x(nt,nej)=xBE(nGeom,nej)                                           
			y(nt,nej)=yBE(nGeom,nej)                                                                                            
			z(nt,nej)=zBE(nGeom,nej)   
		end do  
	end do                                                          
	               
return                                                            
end                                                               
