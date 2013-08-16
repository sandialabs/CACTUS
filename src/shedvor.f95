SUBROUTINE shedvor()   

    use configr
    use blade

    ! CALCULATE THE STRENGTHS OF THE VORTICIES TO BE SHED                

    NTp1=NT+1                                                          
    do I=1,NB                                                      
        NEI=1+(I-1)*(NBE+1)                                               
        NEJ=NEI+NBE                                                       

		! CALCULATE THE STRENGTHS FOR THE SPANWISE VORTICIES               

        NEJ1=NEJ-1                                                        
        do J=1,NBE                                                     
            K=NEI+(J-1)  
			! Set the initial values of the bound vorticity for the next time step GS(NT+1,Element)
			! to the current bound vorticity values. Create shed spanwise vorticity at the current
			! bound vortex location GS(NT,Element).                                                  
            GS(NTp1,K)=GB(K)                                                    
            GS(NT,K)=OGB(K)-GB(K)                                            
            OGB(K)=GB(K)                                                      
		end do

		! Note that the trailing vorticity exists only in the wake (behind the blade) at any given time step.
		! While iterating the non linear system for the bound vorticity, this trailing
		! vorticity (as well as the spanwise shed vorticity) is left constant...
		! Here we set the strength for the new trailing vorticity, originating at the current
		! bound vorticity location, to be included in the velocity calculations in the next time step.

		! GET STRENGTHS OF TRAILING VORTICIES FOR BLADE ENDS

        GT(NT,NEI)=GB(NEI)                                               
        GT(NT,NEJ)=-GB(NEJ1)                                             
        KSTRT=NEI+1                                                       

		! CALCULATE STRENGTHS OF TRAILING VORTICIES FOR THE BLADE SEGMENTS                                                         

        do K=KSTRT,NEJ1                                                
            K1=K-1                                                            
            GT(NT,K)=GB(K)-GB(K1)                                            
		end do
	end do

    RETURN                                                            
END SUBROUTINE shedvor
