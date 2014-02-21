SUBROUTINE EndRev()    

    use util
    use configr
    use output
    use time


    Call cpu_time(time2)

    dtime=time2-time1
    etime=time2-t0
    time1=etime

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
    Output_RevData(1,1)=irev
    Output_RevData(1,2)=CPAve
    Output_RevData(1,3)=KPAve
    Output_RevData(1,4)=CTRAve
    Output_RevData(1,5)=CFxAve
    Output_RevData(1,6)=CFyAve
    Output_RevData(1,7)=CFzAve
    Output_RevData(1,8)=PowerAve
    Output_RevData(1,9)=TorqueAve

    ! Write to revolution average data csv file
    Call csvwrite(9,Output_RevHead,Output_RevData,0,1)

    ! Reset rev average sums
    CPSum=0.0
    CTRSum=0.0
    CFxSum=0.0
    CFySum=0.0
    CFzSum=0.0

    Return                                                              
End SUBROUTINE EndRev
