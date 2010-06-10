SUBROUTINE endrev(cpsum)    

	use configr
	use time
		
	real cpsum
	real nfptf
	integer nfpt										
	
	time2 = secnds(t0)                                                
	
	dtime(irev)=time2-time1                                           
	etime(irev)=time2                                          
	time1=time2                                                       
	cp(irev)=cpsum/nti                                                
	kp(irev)=cp(irev)/ut**3                                           
	write (6,620) cp(irev),kp(irev),irev 
	power=kp(irev)*powerc                                                                                   
								
	! Linear regression curve fit to cp as function of 1/rev     
	! Use last 3 points for the fit...
	nfpt=3
	                                  
	if (irev < (FitStartRev+nfpt-1)) then
		cpf(irev)=0.0                                                     
		kpf(irev)=0.0   
	else                            
		sumx=0.0                                                          
		sumy=0.0                                                          
		sumxy=0.0                                                         
		sumx2=0.0                                                                                                
		do i=(irev-nfpt+1),irev                                                
			x=1.0/float(i)                                                    
			sumx=sumx+x                                                       
			sumy=sumy+cp(i)                                                   
			sumxy=sumxy+cp(i)*x                                               
			sumx2=sumx2+x*x                                                   
		end do 
		nfptf=float(nfpt)                                                         
		b=(sumxy-sumx*sumy/nfptf)/(sumx2-sumx**2/nfptf)                       
		cpf(irev)=sumy/nfptf-b*sumx/nfptf                                     
		kpf(irev)=cpf(irev)/ut**3                                         
	end if                                                     
                                               
Return 
620 FORMAT(10X,'AVERAGE ROTOR CP =',E13.5,', KP =',E13.6,', FOR REVOLUTION NUMBER ',I2)                                                                
End                                                               
