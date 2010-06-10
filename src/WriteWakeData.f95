SUBROUTINE WriteWakeData(tind,DelT,WakeLineInd)   

	! JCM: Print out wake data arrays for viewing in Matlab 

	use wakedata
	use wakeloc
	use configr  
	use wallsoln 
	
	integer :: tind, WakeLineInd(4)
	integer :: wCount
	integer :: tCount, tCountMax, jCount, kCount, yGCErr, nCount
	real :: DelT
	real :: ThetaOut
	real :: VelOut(3,1), IndVel(3,1)
	real, allocatable :: NVelOut(:)
	
	! allocate output
	allocate(NVelOut(NumWP))
	
	if (NT .eq. 1) then
		! Initialize on first run
	
		! Open output files for WriteWakeData
		OPEN(13, FILE='VelOut.dat')   
		! Note: file ID 15 used in input.f95 
		OPEN(16, FILE='VorLine1.dat')  
		OPEN(17, FILE='VorLine2.dat')  
		OPEN(18, FILE='VorLine3.dat') 
		OPEN(19, FILE='VorLine4.dat')  
	
		! Initialize wake deficit to zeros on first run...
		WakeDefAve(:,:)=0  
	end if
	
	! Write wake positions for a couple of wake lines
	! Comma delimited write...
	! TimeStep, Theta, X(1:NT), Y(1:NT), Z(1:NT)
	ThetaOut=(NT-1)*DelT
	tCountMax=NT
	write(16,'I5,",",F10.7,",",$') tCountMax,ThetaOut ! write with commas and no carriage return...
	write(17,'I5,",",F10.7,",",$') tCountMax,ThetaOut ! write with commas and no carriage return...
	write(18,'I5,",",F10.7,",",$') tCountMax,ThetaOut ! write with commas and no carriage return...
	write(19,'I5,",",F10.7,",",$') tCountMax,ThetaOut ! write with commas and no carriage return...
	
	do tCount=1,tCountMax
		write(16,'(F10.7,",",$)') X(tCount,WakeLineInd(1)) ! write with a comma and no carriage return...
		write(17,'(F10.7,",",$)') X(tCount,WakeLineInd(2)) ! write with a comma and no carriage return...
		write(18,'(F10.7,",",$)') X(tCount,WakeLineInd(3)) ! write with a comma and no carriage return...
		write(19,'(F10.7,",",$)') X(tCount,WakeLineInd(4)) ! write with a comma and no carriage return...
	end do
	
	do tCount=1,tCountMax
		write(16,'(F10.7,",",$)') Y(tCount,WakeLineInd(1)) ! write with a comma and no carriage return...
		write(17,'(F10.7,",",$)') Y(tCount,WakeLineInd(2)) ! write with a comma and no carriage return...
		write(18,'(F10.7,",",$)') Y(tCount,WakeLineInd(3)) ! write with a comma and no carriage return...
		write(19,'(F10.7,",",$)') Y(tCount,WakeLineInd(4)) ! write with a comma and no carriage return...
	end do
	
	do tCount=1,tCountMax
		if (tCount .ne. tCountMax) then
			write(16,'(F10.7,",",$)') Z(tCount,WakeLineInd(1))  ! write with a comma and no carriage return...
			write(17,'(F10.7,",",$)') Z(tCount,WakeLineInd(2))  ! write with a comma and no carriage return...
			write(18,'(F10.7,",",$)') Z(tCount,WakeLineInd(3))  ! write with a comma and no carriage return...
			write(19,'(F10.7,",",$)') Z(tCount,WakeLineInd(4))  ! write with a comma and no carriage return...
		else 
			write(16,'F10.7') Z(tCount,WakeLineInd(1)) ! write with a carriage return (default)
			write(17,'F10.7') Z(tCount,WakeLineInd(2)) ! write with a carriage return (default)
			write(18,'F10.7') Z(tCount,WakeLineInd(3)) ! write with a carriage return (default)
			write(19,'F10.7') Z(tCount,WakeLineInd(4)) ! write with a carriage return (default)
		end if
	end do
		
	if (irev .eq. nr) then
	
		if (tind .eq. NTI) then
			
			! Calc wall normal velocities
			if (GPFlag == 1) then
			
				do wCount=1,NumWP
				
					! Freestream
					Call CalcFreestream(WCPoints(2,wCount),VelOut(1,1),VelOut(2,1),VelOut(3,1),ygcErr) 
					
					! Calculate wall and wake induced velocities at wake locations
					
					Call CalcIndVel(nt,ntTerm,nbe,nb,ne,WCPoints(1:3,wCount),IndVel)
					VelOut=VelOut+IndVel
					
					NVelOut(wCount)=sum(WZVec(1:3,wCount)*VelOut(1:3,1))
	
				end do
	
				! Comma delimited write...
				nCount=1
				do jCount=1,NumWPx
					do kCount=1,NumWPx
						if (kCount .ne. NumWPx) then
							write(13,'(F10.7,",",$)') NVelOut(nCount)  ! write with a comma and no carriage return...
						else 
							write(13,'F10.7') NVelOut(nCount) ! write with a carriage return (default)
						end if
						nCount=nCount+1
					end do
				end do
			
			end if 
		
		end if
	
	end if	
		
Return
End	    
