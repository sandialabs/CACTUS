MODULE regtest

	! Regression test outputs

    integer :: RegTFlag             ! Set to 1 to perform a single iteration regression test, 0 for normal operation  
    integer :: Reg_TS               ! timestep number
	integer :: Reg_NLIter	        ! non linear convergence loop iteration
    integer :: Reg_ElemNum          ! element number
    integer :: Reg_DFL              ! BV lift dynamic stall flag
    integer :: Reg_LBC              ! LB dynamic stall model logic checksum

    real :: Reg_CPOut               ! cp after one time step
    real :: Reg_ElemAOA             ! AOA on each element
    real :: Reg_ElemCirc            ! Circulation on each element
    real :: Reg_dElemCirc           ! delta circulation between NL iterations
    real :: Reg_MaxWS               ! max wall source value
    real :: Reg_MaxWVM              ! max velocity mag in first wake shed

End MODULE regtest
