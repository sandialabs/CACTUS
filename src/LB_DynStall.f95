SUBROUTINE LB_DynStall(nElem,CLstat,CDstat,alphaL,alpha5,umach,Re,SectInd,CL,CD)

    use airfoil
    use dystl
    use pidef

    implicit none

    real :: CLstat, CDstat, alphaL, alpha5, umach, Re, CL, CD
    integer :: SectInd, nElem

    real :: AOA0, CLID, Trans, dCLRefLE, dAOARefLE, AOARefLE, CLstatF, C, C1, CLIDF, CLRatio, CLsep, CLF, dCDF, KD, CLa, NOF, dCLv, dCDv, acut

    ! Leishman-Beddoes dynamic stall model (incompressible reduction).

    ! Airfoil data
    AOA0=alzer(SectInd)
    Call CalcLBStallAOALim(Re,SectInd,CLa,CLCritP(nElem),CLCritN(nElem))

    ! Model constants
    KD=.1                  ! TE separation drag factor

    ! Evaluate ideal CL curve at current AOA
    Call LB_EvalIdealCL(alphaL,AOA0,CLa,1,CLRef(nElem))
    Call LB_EvalIdealCL(alphaL,AOA0,CLa,0,CLID)

    ! calc lagged ideal CL for comparison with critical LE separation CL
    Trans=(cos(alphaL-AOA0))**2 ! fair effect to zero at 90 deg. AOA...
    dCLRefLE=Trans*dp(nElem)  ! dp is lagged CLRef change
    dAOARefLE=dCLRefLE/CLa

    ! define reference LE CL and AOA
    CLRefLE(nElem)=CLRef(nElem)-dCLRefLE
    if (CLRefLE(nElem)*(CLRefLE(nElem)-CLRefLE_Last(nElem)) > 0) then
        CLRateFlag(nElem)=1
    else
        CLRateFlag(nElem)=0
    end if
    AOARefLE=alphaL-dAOARefLE
    Call Force180(AOARefLE)

    ! calc effective static TE separation point using effective LE AOA
    Call intp(Re,AOARefLE*condeg,CLstatF,C,C1,SectInd)
    Call LB_EvalIdealCL(AOARefLE,AOA0,CLa,0,CLIDF)
    if (abs(CLIDF)<0.001) then
        CLRatio=999
    else
        CLRatio=CLstatF/CLIDF;
    end if

    if (CLRatio > 0.25) then
        Fstat(nElem)=min((sqrt(4.0*CLRatio)-1.0)**2,1.0)

        ! Test logic
        LB_LogicOutputs(nElem,1)=1
    else
        Fstat(nElem)=0

        ! Test logic
        LB_LogicOutputs(nElem,1)=2
    end if
    ! calc lagged Fstat to represent dynamic TE separation point
    F(nElem)=Fstat(nElem)-dF(nElem)
    ! force limits on lagged F (needed due to discretization error...)
    F(nElem)=min(max(F(nElem),0.0),1.0)

    ! Calc dynamic CL due to TE separation as fairing between fully attached and fully separated predictions from the Kirchoff approximation at current AOA
    if (abs(CLID)<0.001) then
        CLRatio=999
    else
        CLRatio=CLstat/CLID
    end if

    if (CLRatio > 1.0) then
        CLID=CLstat

        ! Test logic
        LB_LogicOutputs(nElem,2)=1
    end if

    if (CLRatio > 0.25) then
        CLsep=CLID/4.0

        ! Test logic
        LB_LogicOutputs(nElem,3)=1
    else
        CLsep=CLstat

        ! Test logic
        LB_LogicOutputs(nElem,3)=2
    end if
    CLF=CLsep+CLID*0.25*(F(nElem)+2.0*sqrt(F(nElem)))
    dCDF=KD*(CLstat-CLF)*sign(1.0,CLstat)

    ! LE vortex lift component, dCNv is a lagged change in the added normal force due
    ! to LE vortex shedding. Assumed to affect lift coeff as an added circulation...
    dCLv=dCNv(nElem)*cos(alpha5)
    dCDv=dCNv(nElem)*sin(alpha5)
    ! vortex feed is given by the rate at which lift (circulation) is being shed due to dynamic separation. Lift component due to separation is defined by the
    ! difference between the ideal lift and the lift including dynamic separation effects.
    cv(nElem)=CLID-CLF
    dcv(nElem)=cv(nElem)-cv_Last(nElem)
    ! If the sign of dcv is opposite the reference LE CL, set to zero to disallow negative vorticity from shedding from the leading edge. Also, limit the model 
    ! at AOA>acut or if the magnitude of the reference CL is decreasing...
    acut=50.0*conrad
    if (sign(1.0,dcv(nElem)*CLRefLE(nElem))<0 .OR. abs(alphaL-AOA0)>acut .OR. CLRateFlag(nElem)<0) then
        dcv=0.0

        ! Test logic
        LB_LogicOutputs(nElem,4)=1
    end if

    ! Total lift and drag
    CL=CLF+dCLv
    CD=CDstat+dCDF+dCDv

    Return
End SUBROUTINE LB_DynStall
