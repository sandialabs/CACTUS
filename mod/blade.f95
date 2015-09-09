MODULE blade

    ! Blade data

    ! Blade loading data
    real, allocatable :: GB(:)              ! Bound vorticity
    real, allocatable :: OGB(:)             ! Old bound vorticity (previous time step)
    real, allocatable :: GB_Raw(:)          ! Raw (pre-filter) bound vorticity
    real, allocatable :: AOA(:)             ! AOA on blade elements
    real, allocatable :: AOA_Last(:)        ! Last value of AOA on blade elements

    real, allocatable :: UIWake(:)          ! Velocity induced at blade from wake
    real, allocatable :: VIWake(:)          ! Velocity induced at blade from wake
    real, allocatable :: WIWake(:)          ! Velocity induced at blade from wake

    ! Induced velocity on blade lattice points
    real, allocatable :: UB(:)             ! Lattice point x velocity for each element end lattice point
    real, allocatable :: VB(:)             ! Lattice point y velocity for each element end lattice point
    real, allocatable :: WB(:)             ! Lattice point z velocity for each element end lattice point

    ! Freestream velocity on blade lattice points
    real, allocatable :: UFSB(:)             ! Lattice point x velocity for each element end lattice point
    real, allocatable :: VFSB(:)             ! Lattice point y velocity for each element end lattice point
    real, allocatable :: WFSB(:)             ! Lattice point z velocity for each element end lattice point

    ! Timestep filter
    integer :: TSFilFlag            ! 1 to enable timestep filtering, 0 for no filtering (default)
    integer :: ntsf                 ! Number of timesteps over which the bound vorticity is filtered smooth (if TSFilFlag = 1)
    real    :: KTF

CONTAINS

    SUBROUTINE blade_cns(MaxSegEnds)

        ! Constructor - allocates memory for arrays containing element induced velocities on all blades

        integer :: MaxSegEnds

        allocate(GB(MaxSegEnds))
        allocate(OGB(MaxSegEnds))
        allocate(GB_Raw(MaxSegEnds))
        allocate(AOA(MaxSegEnds))
        allocate(AOA_Last(MaxSegEnds))
        allocate(UIWake(MaxSegEnds))
        allocate(VIWake(MaxSegEnds))
        allocate(WIWake(MaxSegEnds))
        allocate(UB(MaxSegEnds))
        allocate(VB(MaxSegEnds))
        allocate(WB(MaxSegEnds))
        allocate(UFSB(MaxSegEnds))
        allocate(VFSB(MaxSegEnds))
        allocate(WFSB(MaxSegEnds))


    End SUBROUTINE blade_cns

    SUBROUTINE UpdateAOALast(ne)

        integer :: ne
        integer :: k

        ! Save last AOA values for each element

        do k=1,ne
            AOA_Last(k)=AOA(k)
        end do

    End SUBROUTINE UpdateAOALast

    SUBROUTINE UpdateTSFilter(ne)

        integer :: ne
        integer :: k

        ! Update filtered bound vorticity (filtered smooth over approx. ntsf timesteps using a first order discrete filter)

        do k=1,ne
            GB(k)=KTF*GB_Raw(k) + (1.0-KTF)*GB(k)
        end do

    End SUBROUTINE UpdateTSFilter

End MODULE blade
