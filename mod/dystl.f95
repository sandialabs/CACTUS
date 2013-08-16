MODULE dystl

    ! Dynamic stall data

    ! Module private data
    private :: pi, conrad, condeg

    ! Pi definition
    real :: pi
    real :: conrad
    real :: condeg

    integer :: DSFlag                               ! 0 for no dynamic stall, 1 for BV model, 2 for LB model


    ! Modified Boeing-Vertol (Gormont) model      
    real :: K1Pos                                   ! lagged AOA magnitude tweak for CL increasing
    real :: K1Neg                                   ! lagged AOA magnitude tweak for CL decreasing
    integer, allocatable :: BV_DynamicFlagL(:)      ! Dynamic stall active flags
    integer, allocatable :: BV_DynamicFlagD(:)      ! Dynamic stall active flags

    ! Additional BV diagnostic output     
    real :: BV_alpha, BV_adotnorm, BV_alrefL, BV_alrefD     


    ! Leishman-Beddoes model

    ! Time step normalized as ds = 2*U*dt/c where U is the relative velocity and c is the chord (note: at the same t, s is different on each element)
    real, allocatable :: ds(:)                      

    ! Temporally updated state variables (held for each element)
    real, allocatable :: dp(:)
    real, allocatable :: dF(:)
    real, allocatable :: dCNv(:)
    real, allocatable :: sLEv(:) 
    integer, allocatable :: LESepState(:)

    ! Other states needed at the program level and discrete lagged values used in the state update process (held for each element)
    real, allocatable :: CLRef(:) 
    real, allocatable :: CLRefLE(:) 
    real, allocatable :: CLCritP(:)
    real, allocatable :: CLCritN(:)
    integer, allocatable :: CLRateFlag(:)
    real, allocatable :: Fstat(:)
    real, allocatable :: F(:)
    real, allocatable :: cv(:)
    real, allocatable :: dcv(:)      
    real, allocatable :: CLRef_Last(:) 
    real, allocatable :: CLRefLE_Last(:) 
    real, allocatable :: Fstat_Last(:)
    real, allocatable :: cv_Last(:)

    ! Additional LB diagnostic output
    integer, parameter :: NLBL = 9
    integer, allocatable :: LB_LogicOutputs(:,:) 
    integer :: Logic_W(NLBL) = [1,1,1,1,1,3,2,1,1]             ! Logic weights for logic checksum  


CONTAINS

    SUBROUTINE dystl_cns(MaxAirfoilSect, MaxReVals, MaxSegEnds)

        ! Constructor for the arrays in this module

        integer :: MaxAirfoilSect, MaxReVals, MaxSegEnds

        ! Pi definition
        pi = 4.0*atan(1.0)
        conrad = pi/180.0                                                  
        condeg = 180.0/pi         

        ! Boeing-Vertol
        allocate(BV_DynamicFlagL(MaxSegEnds))
        allocate(BV_DynamicFlagD(MaxSegEnds))

        ! Leishman-Beddoes
        allocate(ds(MaxSegEnds))
        allocate(dp(MaxSegEnds))
        allocate(dF(MaxSegEnds))
        allocate(dCNv(MaxSegEnds))
        allocate(sLEv(MaxSegEnds))
        allocate(LESepState(MaxSegEnds))

        allocate(CLRef(MaxSegEnds))
        allocate(CLRefLE(MaxSegEnds))
        allocate(CLCritP(MaxSegEnds))
        allocate(CLCritN(MaxSegEnds))
        allocate(CLRateFlag(MaxSegEnds))
        allocate(Fstat(MaxSegEnds))
        allocate(F(MaxSegEnds))
        allocate(cv(MaxSegEnds))
        allocate(dcv(MaxSegEnds))      
        allocate(CLRef_Last(MaxSegEnds))
        allocate(CLRefLE_Last(MaxSegEnds)) 
        allocate(Fstat_Last(MaxSegEnds))
        allocate(cv_Last(MaxSegEnds))

        allocate(LB_LogicOutputs(MaxSegEnds,NLBL))

    End SUBROUTINE dystl_cns

    SUBROUTINE dystl_init_LB()

        ! Initialize LB model
        dp(:) = 0.0
        dF(:) = 0.0
        dCNv(:) = 0.0
        LESepState(:) = 0
        sLEv(:) = 0.0
        CLRef_Last(:) = 0.0
        CLRefLE_Last(:) = 0.0
        Fstat_Last(:) = 1.0
        cv_Last(:) = 0.0

        LB_LogicOutputs(:,:)=0

    End SUBROUTINE dystl_init_LB

    SUBROUTINE dystl_init_BV()

        BV_DynamicFlagL(:)=0
        BV_DynamicFlagD(:)=0

    End SUBROUTINE dystl_init_BV

    SUBROUTINE LB_EvalIdealCL(AOA,AOA0,CLa,RefFlag,CLID)

        ! AOA inputs in radians
        ! AOA0 is zero lift AOA
        ! RefFlag defines whether to output reference CL or ideal CL
        ! CLa is reference lift slope (per radian) to be used for reference CL (ideal CLa is 2*pi)

        real :: AOA, AOA0, CLa
        integer :: RefFlag
        real :: CLID

        real :: IDS, aID, d1, CLaI, aIDc

        aID=AOA-AOA0
        Call Force180(aID)

        ! reflect function across axis for abs(aID)>90
        IDS=1
        if (aID>pi/2.0) then
            aID=pi-aID
            IDS=-1.0
        else if (aID<-pi/2.0) then
            aID=-pi-aID
            IDS=-1.0
        end if

        ! If RefFlag is 1, output reference CL, otherwise round off ideal CL at high AOA
        if (RefFlag==1) then
            CLID=IDS*CLa*aID
        else
            ! round off the ideal CL after cutoff AOA
            aIDc=30.0*conrad
            d1=1.8
            CLaI=2.0*pi
            if (abs(aID)<aIDc) then
                CLID=IDS*CLaI*aID
            else
                CLID=IDS*(CLaI*(aIDc-1.0/d1*sin(d1*aIDc))+CLaI/d1*sin(d1*aID))
            end if
        end if

    End SUBROUTINE LB_EvalIdealCL

    SUBROUTINE Force180(a)

        real :: a

        ! alpha in radians
        if (a>pi) then
            a=a-2.0*pi
        else if (a<-pi) then
            a=a+2.0*pi
        end if

    End SUBROUTINE Force180

    SUBROUTINE LB_UpdateStates(nb,nbe)

        ! Update states for the LB model 
        ! Note dynstall should be included eventually in an expanded blade module, at which point it would have 
        ! access to the geometry info it needs...

        integer :: nb, nbe      

        integer :: i, nei, j, nElem, IsBE

        ! Set model parameters. All of these are potentially a function of Mach
        ! and are set to low mach values...
        Tp=1.7                  ! time constant on LE pressure response to change in CL
        TfRef=3.0               ! time constant on TE separation point travel
        TvRef=6.0               ! time constant on LE vortex lift indicial function
        TvL=11.0                ! Characteristic LE vortex travel time

        ! Update states for each blade element
        do i=1,nb                                                                                                     
            nei=1+(i-1)*(nbe+1)                                               
            do j=1,nbe  

                ! Blade element is referenced by its upper end location                                           
                nElem=nei+j                                                         

                IsBE=0
                if (j==1 .OR. j==nbe) then
                    IsBE=1
                end if

                ! Dynamic stall model not active for blade end elements
                if (IsBE==0) then

                    ! Eval LE separation state
                    ! Note: CLCrit is critical (ideal) CL value for LE separation. This is
                    ! approximately equal to the CL that would exist at the angle of attack at
                    ! max CL if the CL curve had remained linear
                    if (LESepState(nElem)==0 .AND. (CLRefLE(nElem)>CLCritP(nElem) .OR. CLRefLE(nElem)<CLCritN(nElem))) then
                        ! In LE separation state
                        LESepState(nElem)=1
                        sLEv(nElem)=0 ! reset leading edge vortex time counter

                        ! Set logic state flags (for model diagnosis output)
                        LB_LogicOutputs(nElem,5)=1
                    else if (LESepState(nElem)==1 .AND. (CLRefLE(nElem)<CLCritP(nElem) .AND. CLRefLE(nElem)>CLCritN(nElem))) then
                        ! Out of LE separation state
                        LESepState(nElem)=0
                        sLEv(nElem)=0 ! reset leading edge vortex time counter

                        ! Set logic state flags (for model diagnosis output)
                        LB_LogicOutputs(nElem,5)=2
                    end if

                    ! Set time constants based on LE separation state and TE separation point
                    ! location. Different time constants for abs(CL) increasing vs. decreasing
                    if (LESepState(nElem)==1) then
                        if (sLEv(nElem)<TvL) then
                            if (CLRateFlag(nElem)>0) then
                                !Tf=3.0*TfRef ! original
                                Tf=4.0*TfRef
                                Tv=TvRef

                                ! Set logic state flags (for model diagnosis output)
                                LB_LogicOutputs(nElem,8)=1
                            else
                                Tf=1.0/2.0*TfRef
                                Tv=1.0/2.0*TvRef

                                ! Set logic state flags (for model diagnosis output)
                                LB_LogicOutputs(nElem,8)=2
                            end if

                            ! Set logic state flags (for model diagnosis output)
                            LB_LogicOutputs(nElem,7)=1
                        else if (sLEv(nElem)<2.0*TvL) then
                            if (CLRateFlag(nElem)>0) then
                                ! orig 
                                !Tf=1.0/3.0*TfRef
                                !Tv=1.0/4.0*TvRef
                                Tf=2.0*TfRef
                                Tv=TvRef

                                ! Set logic state flags (for model diagnosis output)
                                LB_LogicOutputs(nElem,8)=3
                            else
                                Tf=1.0/2.0*TfRef
                                Tv=1.0/2.0*TvRef

                                ! Set logic state flags (for model diagnosis output)
                                LB_LogicOutputs(nElem,8)=4                                                        
                            end if

                            ! Set logic state flags (for model diagnosis output)
                            LB_LogicOutputs(nElem,7)=2                                                
                        else
                            ! orig
                            !Tf=4.0*TfRef
                            !Tv=0.9*TvRef
                            Tf=TfRef
                            Tv=TvRef

                            ! Set logic state flags (for model diagnosis output)
                            LB_LogicOutputs(nElem,7)=3
                        end if

                        ! Set logic state flags (for model diagnosis output)
                        LB_LogicOutputs(nElem,6)=1
                    else
                        if (F(nElem)>0.7) then
                            Tf=TfRef

                            ! Set logic state flags (for model diagnosis output)
                            LB_LogicOutputs(nElem,7)=4
                        else
                            Tf=2*TfRef

                            ! Set logic state flags (for model diagnosis output)
                            LB_LogicOutputs(nElem,7)=5
                        end if
                        Tv=TvRef

                        ! Set logic state flags (for model diagnosis output)
                        LB_LogicOutputs(nElem,6)=2
                    end if

                    ! update LE vortex time counter if in LE separation state
                    if (LESepState(nElem)==1) then
                        sLEv(nElem)=sLEv(nElem)+ds(nElem)

                        ! Set logic state flags (for model diagnosis output)
                        LB_LogicOutputs(nElem,9)=1
                    end if

                    ! Update states, first order lag equations, exponential recursion form (midpoint rule version)
                    dp(nElem)=dp(nElem)*exp(-ds(nElem)/Tp)+(CLRef(nElem)-CLRef_Last(nElem))*exp(-ds(nElem)/(2*Tp))
                    dF(nElem)=dF(nElem)*exp(-ds(nElem)/Tf)+(Fstat(nElem)-Fstat_Last(nElem))*exp(-ds(nElem)/(2*Tf))
                    dCNv(nElem)=dCNv(nElem)*exp(-ds(nElem)/Tv)+dcv(nElem)*exp(-ds(nElem)/(2*Tv))

                    ! update lagged values
                    CLRef_Last(nElem)=CLRef(nElem)
                    CLRefLE_Last(nElem)=CLRefLE(nElem)
                    Fstat_Last(nElem)=Fstat(nElem)
                    cv_Last(nElem)=cv(nElem)

                end if
            end do

        end do

    End SUBROUTINE LB_UpdateStates

    SUBROUTINE LB_LogicChecksum(nElem,LBCheck)

        integer :: LBCheck, nElem
        integer :: Loop

        ! Calculates a checksum for the logic states in the LB model

        LBCheck=0
        do Loop=1,NLBL
            LBCheck=LBCheck+LB_LogicOutputs(nElem,Loop)*Logic_W(Loop)
        end do

    End SUBROUTINE LB_LogicChecksum

End MODULE dystl
