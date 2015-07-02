subroutine VorIVel(VFlag,CalcDer,G,X1,Y1,Z1,X2,Y2,Z2,XP,YP,ZP,UP,VP,WP,DUDX)      

    use vortex

    real :: DVFDX, RRAT
    integer :: VFlag       ! 0 for bound vortex element, 1 for trailing wake element, 2 for spanwise wake element (sets core radius)
    integer :: CalcDer     ! 1 to calc dudx

    AX = X2-X1 
    AY = Y2-Y1 
    AZ = Z2-Z1
    A2 = AX**2 + AY**2 + AZ**2

    BX = X2-XP                                                                                             
    BY = Y2-YP                                                 
    BZ = Z2-ZP                                        
    B = SQRT(BX**2 + BY**2 + BZ**2)                                   

    CX = X1-XP                                                     
    CY = Y1-YP                                                     
    CZ = Z1-ZP                                                                                                                        
    C = SQRT(CX**2 + CY**2 + CZ**2) 

    if (B == 0.0) then 
        B = 1.0
    end if
    if (C == 0.0) then
        C = 1.0
    end if

    ADB = (AX*BX + AY*BY + AZ*BZ)
    ADBDB = ADB/B
    ADC = (AX*CX + AY*CY + AZ*CZ)
    ADCDC = ADC/C

    CCAX = CY*AZ - AY*CZ                                              
    CCAY = AX*CZ - CX*AZ                                              
    CCAZ = CX*AY - AX*CY                                              

    CCAV = CCAX**2 + CCAY**2 + CCAZ**2                                

    ! Select the vortex core radius for bound vortex, trailing
    ! or spanwise wake elements
    if (VFlag == 1) then
        VRAD2=VRAD2_T
    else if (VFlag == 2) then
        VRAD2=VRAD2_S
    else
        VRAD2=VRAD2_B
    end if

    if (CCAV >= vCutOffRad) then
        if (ivtxcor > 0 .and. CCAV < A2*VRAD2) then
            if (ivtxcor == 1) then
                ! Constant velocity (approx) core
                CCAV = A2*VRAD2 ! limit denominator to value at core radius
                VF   = (ADBDB-ADCDC)*G/(12.56637*CCAV)
                if (CalcDer == 1) then
                    DVFDX = G/12.56637*((AX/C-CX*ADC/C**3-AX/B+BX*ADB/B**3)/CCAV)
                end if
            else if (ivtxcor == 2) then
                ! Linear velocity (approx) core
                RRAT=sqrt(CCAV/(A2*VRAD2)) ! radius to core radius ratio
                CCAV = A2*VRAD2 ! limit denominator to value at core radius
                VF   = RRAT*(ADBDB-ADCDC)*G/(12.56637*CCAV)
                if (CalcDer == 1) then
                    DVFDX = G/12.56637*((AX/C-CX*ADC/C**3-AX/B+BX*ADB/B**3)/CCAV)
                end if
            end if
        else
            VF   = (ADBDB-ADCDC)*G/(12.56637*CCAV)
            if (CalcDer == 1) then
                DVFDX = G/12.56637*((AX/C-CX*ADC/C**3-AX/B+BX*ADB/B**3)/CCAV-2*(ADBDB-ADCDC)/CCAV**2*(AZ*CCAY-AY*CCAZ))
            end if
        end if
    else
        VF   = 0.0
        DVFDX = 0.0
    end if

    UP = UP+CCAX*VF                                                 
    VP = VP+CCAY*VF                                                 
    WP = WP+CCAZ*VF
    if (CalcDer == 1) then
        DUDX = DUDX+CCAX*DVFDX   
    end if

    RETURN                                                            
END subroutine VorIVel
