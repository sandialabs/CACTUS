SUBROUTINE WriteFinalOutput()
        
	use time
        use output    
        use airfoil, only : ilxtp, iuxtp   
        
        implicit none
        
        include 'csvwrite.inc'
                                                                                                                          
	! Check error flags and write notifications to stdout                                
	if (ilxtp .gt. 0) write (6,615)
	if (iuxtp .gt. 0) write (6,618)

        ! Write revolution average data csv file
        Call csvwrite(9,Output_RevHead,Output_RevData,Output_RevRow)
        
        ! Write timestep data csv file
        Call csvwrite(10,Output_TSHead,Output_TSData,Output_TSRow)
        
        ! Write element loads data csv file
        if (Output_ELFlag == 1) then
                Call csvwrite(11,Output_ELHead,Output_ELData,Output_ELRow)
        end if

Return
618  FORMAT ('AT LEAST ONE BLADE REYNOLDS NUMBER WAS ABOVE TABLE LIMIT.  UPPER LIMIT USED.')
615  FORMAT ('AT LEAST ONE BLADE REYNOLDS NUMBER WAS BELOW TABLE LIMIT.  LOWER LIMIT USED.')  
End