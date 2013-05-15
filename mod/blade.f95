MODULE blade
        
        ! Blade and blade wake data
        
        ! Blade loading data
        real, allocatable :: GB(:)              ! Bound vorticity
        real, allocatable :: OGB(:)             ! Old bound vorticity (previous time step)
        real, allocatable :: GB_Raw(:)          ! Raw (pre-filter) bound vorticity
        real, allocatable :: AOA(:)             ! AOA on blade elements    
        real, allocatable :: AOA_Last(:)        ! Last value of AOA on blade elements
         
        real, allocatable :: UIWake(:)          ! Velocity induced at blade from wake
        real, allocatable :: VIWake(:)          ! Velocity induced at blade from wake
        real, allocatable :: WIWake(:)          ! Velocity induced at blade from wake
        
        ! Wake circulation data
        real, allocatable :: GT(:,:)            ! Trailing wake (streamwise) vorticity
        real, allocatable :: GS(:,:)            ! Shed wake (spanwise) vorticity
        
        ! Wake location data 
        real, allocatable :: X(:,:)             ! X position
        real, allocatable :: Y(:,:)             ! Y position
        real, allocatable :: Z(:,:)             ! Z position
        
        ! Wake lattice point velocity data
        real, allocatable :: U(:,:)             ! Lattice point x velocity for each wake point
        real, allocatable :: V(:,:)             ! Lattice point y velocity for each wake point
        real, allocatable :: W(:,:)             ! Lattice point z velocity for each wake point
        
        ! Old wake lattice point velocities
        real, allocatable :: UO(:,:)            ! Last lattice point x velocity for each wake point
        real, allocatable :: VO(:,:)            ! Last lattice point y velocity for each wake point
        real, allocatable :: WO(:,:)            ! Last lattice point z velocity for each wake point
        
        ! Timestep filter
        integer :: TSFilFlag            ! 1 to enable timestep filtering, 0 for no filtering (default) 
        integer :: ntsf                 ! Number of timesteps over which the bound vorticity is filtered smooth (if TSFilFlag = 1)
        real    :: KTF                
                
        CONTAINS

        SUBROUTINE blade_cns(MaxWakeNodes, MaxSegEnds)

                ! Constructor for the arrays in this module

                integer :: MaxWakeNodes, MaxSegEnds
        
                allocate(GB(MaxSegEnds))            
                allocate(OGB(MaxSegEnds))
                allocate(GB_Raw(MaxSegEnds))
                allocate(AOA(MaxSegEnds))   
                allocate(AOA_Last(MaxSegEnds))   
                allocate(UIWake(MaxSegEnds))
                allocate(VIWake(MaxSegEnds))
                allocate(WIWake(MaxSegEnds))
                allocate(GT(MaxWakeNodes,MaxSegEnds))
                allocate(GS(MaxWakeNodes+1,MaxSegEnds))         ! needs extra spanwise station for shedvor
                allocate(X(MaxWakeNodes,MaxSegEnds))
                allocate(Y(MaxWakeNodes,MaxSegEnds))
                allocate(Z(MaxWakeNodes,MaxSegEnds))
                allocate(U(MaxWakeNodes,MaxSegEnds))
                allocate(V(MaxWakeNodes,MaxSegEnds))
                allocate(W(MaxWakeNodes,MaxSegEnds))
                allocate(UO(MaxWakeNodes,MaxSegEnds))
                allocate(VO(MaxWakeNodes,MaxSegEnds))
                allocate(WO(MaxWakeNodes,MaxSegEnds))               
                
        End SUBROUTINE
        
        SUBROUTINE UpdateAOALast(ne)
        
                integer :: ne
                integer :: k
        
                ! Save last AOA values for each element
        
                do k=1,ne                                                                                                             
                        AOA_Last(k)=AOA(k)                                                         
                end do 
        
        End SUBROUTINE
        
        SUBROUTINE UpdateTSFilter(ne)
        
                integer :: ne
                integer :: k
        
                ! Update filtered bound vorticity (filtered smooth over approx. ntsf timesteps using a first order discrete filter)
        
                do k=1,ne                                                                                                             
                        GB(k)=KTF*GB_Raw(k) + (1.0-KTF)*GB(k)                                                        
                end do 
        
        End SUBROUTINE        

End