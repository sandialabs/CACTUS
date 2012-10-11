SUBROUTINE InputGeom(FN)
        
        use element
        use varscale 
        use configr    
        
        implicit none
               
        integer, parameter :: MaxReadLine = 1000    
        character(MaxReadLine) :: FN    ! path to geometry input file 
        
        integer :: NElem, i
        character(MaxReadLine) :: ReadLine
                                                                                                                                                                                                                                     
        ! Read geometry data file
        
        ! Format Example:
!         NBlade: 3
!         RotN: 0 0 1
!         RotP: 0 0 0
!         RefAR: 2.0
!         RefR: 10.0
!         Type: VAWT_Par
!         Blade 1:
!             NElem: 5
!             QCx: 0 0 0 0 0 0
!             QCy: 1 2 3 4 5 6
!             QCz: 1 1 1 1 1 1
!             nx: 0 0 0 0 0 0
!             ny: 0 0 0 0 0 0
!             nz: -1 -1 -1 -1 -1 -1
!             tx: 1 1 1 1 1 1
!             ty: 0 0 0 0 0 0
!             tz: 0 0 0 0 0 0
!             CtoR: .1 .1 .1 .1 .1 .1
!             AreaR: .1 .1 .1 .1 .1
!             iSect: 1 1 1 1 1
!         Blade 2:
!             ...  
               
        ! Open input file for this section
        open(15, file=FN)
        
        ! Read header data
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) nb
        
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) RotX, RotY, RotZ
        
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) RotPX, RotPY, RotPZ
        
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) at
        
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Rmax
        
        read(15,'(A)') ReadLine
        
        ! Allocate blades struct
        allocate(Blades(nb))
        
        ! Read blade data
        do i=1,nb
        
            read(15,'(A)') ReadLine
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) NElem
            
            ! Allocate arrays in blade structure
            Call blade_geom_cns(i,NElem)
            
            ! JCM: Currently (until blade/wake geometry module is restructured), NElem must be the same for all blades.
            nbe=NElem
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%QCx(1:NElem+1)
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%QCy(1:NElem+1)
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%QCz(1:NElem+1)  
                      
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%nx(1:NElem+1)                      
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%ny(1:NElem+1) 
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%nz(1:NElem+1) 
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tx(1:NElem+1)                      
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%ty(1:NElem+1) 
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%tz(1:NElem+1) 
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%CtoR(1:NElem+1) 
            
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%AreaR(1:NElem) 
       
            read(15,'(A)') ReadLine
            read(ReadLine(index(ReadLine,':')+1:),*) Blades(i)%iSect(1:NElem) 
        
        end do 
        
        ! Close input file 
        close(15)
                               
Return
End