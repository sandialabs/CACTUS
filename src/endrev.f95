SUBROUTINE endrev(cpave)    

	use configr
        use output      
	use time
		
        real cpave        
	
        real kpave, trave, dtime, etime										
	
	time2 = secnds(t0)                                                
	
	dtime=time2-time1                                           
	etime=time2                                          
	time1=time2  
                                                                                                             
        ! Calc average power over last revolution
        cpave=0.0
        trave=0.0
        do i=(Output_TSRow-nti+1),Output_TSRow
                ! Average torque and power coefficients based on freestream speed and Rmax
                trave=trave+Output_TSData(i,3)/nti
                cpave=cpave+Output_TSData(i,4)/nti        
        end do
        ! Torque in ft-lbs
        torque=trave*torquec
        ! Power coefficient based on tip speed
        kpave=cpave/ut**3
        ! Power in kW                                                      
        power=kpave*powerc                                                      
                                                             
        ! Set revolution average output 
        Output_RevRow=irev                                                    
        Output_RevData(irev,1)=irev
        Output_RevData(irev,2)=cpave
        Output_RevData(irev,3)=kpave
        Output_RevData(irev,4)=power
        Output_RevData(irev,5)=torque
        Output_RevData(irev,6)=dtime
        Output_RevData(irev,7)=etime                                                   							
                                          
Return                                                              
End                                                               
