subroutine EndRev()

    use util
    use configr
    use output
    use time
    use fnames

    Call cpu_time(time2)
!$  time2 = omp_get_wtime()

    dtime=time2-time1
    etime=time2-t0
    time1=etime
!$  time1=omp_get_wtime()

    ! Calc average power over last revolution
    CPAve=CPSum/nti
    CTRAve=CTRSum/nti
    CFxAve=CFxSum/nti
    CFyAve=CFySum/nti
    CFzAve=CFzSum/nti
    ! Torque in ft-lbs
    TorqueAve=CTRAve*TorqueC
    ! Power coefficient based on tip speed
    KPAve=CPAve/ut**3
    ! Power in kW
    PowerAve=KPave*PowerC

    ! Set revolution average output
    RevOutData(1,1)=irev
    RevOutData(1,2)=CPAve
    RevOutData(1,3)=KPAve
    RevOutData(1,4)=CTRAve
    RevOutData(1,5)=CFxAve
    RevOutData(1,6)=CFyAve
    RevOutData(1,7)=CFzAve
    RevOutData(1,8)=PowerAve
    RevOutData(1,9)=TorqueAve

    ! Write to revolution average data csv file
    OPEN(9, FILE=RevOutputFN, POSITION='append')
    Call csvwrite(9,RevOutHead,RevOutData,0,1)
    CLOSE(9)

    ! Reset rev average sums
    CPSum=0.0
    CTRSum=0.0
    CFxSum=0.0
    CFySum=0.0
    CFzSum=0.0

    return
end subroutine EndRev
