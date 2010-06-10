SUBROUTINE WriteFinalOutput()

	use time
	use test
        
	! Final time                  
	TimeF = secnds(t0)                                                                                                                    
	                            
	call listit (1,6)       
	call listit (1,12)       
	if (ilxtp .gt. 0) write (6,615)
	if (iuxtp .gt. 0) write (6,618)
	
Return
618  FORMAT ('AT LEAST ONE BLADE REYNOLDS NUMBER WAS ABOVE TABLE LIMIT.  UPPER LIMIT USED.')
615  FORMAT ('AT LEAST ONE BLADE REYNOLDS NUMBER WAS BELOW TABLE LIMIT.  LOWER LIMIT USED.')   
End