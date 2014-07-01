SUBROUTINE WriteRegTOutput(Flag)
    
    use regtest
    
    integer :: Flag
    if (Flag == 0) then
        write(7,*) 'Timestep NLIter ElemNum ElemAOA ElemCirc dElemCirc BVDynFlagL LBCheck'
    else if (Flag == 1) then
        ! Write data during non-linear iteration 
        write(7,'(3I3,3E13.5,2I3)') Reg_TS, Reg_NLIter, Reg_ElemNum, Reg_ElemAOA, Reg_ElemCirc, Reg_dElemCirc, Reg_DFL, Reg_LBC
    else if (Flag == 2) then
        ! Write cp average over timesteps
        write(7,*) 'CPave'
        write(7,'(E13.5)') Reg_CPOut
        ! Write max calculated wall source strength
        write(7,*) 'MaxWS'
        write(7,'(E13.5)') Reg_MaxWS
        ! Write max calculated velocity in first wake shed
        write(7,*) 'MaxWakeV'
        write(7,'(E13.5)') Reg_MaxWVM
    end if
	
    Return 
End SUBROUTINE WriteRegTOutput
