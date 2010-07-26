SUBROUTINE WriteRegTOutput(Flag)
        
        use regtest
        
        integer :: Flag
        if (Flag == 0) then
                write(7,*) 'Timestep NLIter ElemNum DynFlag ElemAOA ElemCirc dElemCirc'
        else if (Flag == 1) then
                ! Write data during non-linear iteration on first time step
                write(7,'5I3,3E13.5') Reg_TS, Reg_NLIter, Reg_ElemNum, Reg_DFL, Reg_DFD, Reg_ElemAOA, Reg_ElemCirc, Reg_dElemCirc
        else if (Flag == 2) then
                ! Write cp average over timesteps
                write(7,*) 'CPave'
                write(7,'E13.5') Reg_CPOut
                ! Write max calculated wall source strength
                write(7,*) 'MaxWS'
                write(7,'E13.5') Reg_MaxWS
                ! Write max calculated velocity in first wake shed
                write(7,*) 'MaxWakeV'
                write(7,'E13.5') Reg_MaxWVM
        end if
	
Return 
End