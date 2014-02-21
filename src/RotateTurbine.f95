SUBROUTINE RotateTurbine

    use configr       
	use element
	use blade
    use strut       

	! Updates current turbine geometry

    ! Rotate turbine axis of rotation and origin if necessary...                                                         

    ! Rotate blades
    do i=1,nb
        Call RotateBlade(i,delt,RotX,RotY,RotZ,RotPX,RotPY,RotPZ)
    end do

    ! Rotate struts
    do i=1,NStrut
        Call RotateStrut(i,delt,RotX,RotY,RotZ,RotPX,RotPY,RotPZ)
    end do

    return                                                            
end SUBROUTINE RotateTurbine
