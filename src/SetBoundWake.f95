SUBROUTINE SetBoundWake() 

    use configr       
	use element
	use blade

    ! Set new wake element positions
	do i=1,nb                                                      

		nei=1+(i-1)*(nbe+1)                                               

        ! Blade element end locations (quarter chord).
		do j=0,nbe   
			nej=nei+j ! element index
			x(nt,nej)=xBE(nej)                                           
			y(nt,nej)=yBE(nej)                                                                                            
			z(nt,nej)=zBE(nej)   
		end do
	end do


    return                                                            
end SUBROUTINE SetBoundWake
