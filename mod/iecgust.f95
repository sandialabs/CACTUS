MODULE iecgust

	! Parameters defining IEC Gust

    integer :: Igust              ! Flag for IEC Gust
    real :: gustamp		! Gust Amplitude (m/s)
    real :: gusttime		! Gust Timescale (sec)
    real :: gustA                 ! Non-dimensional gust amplitude (U'/Uinf)
    real :: gustT                 ! Non-dimensional gust timescale (T*Uinf/Rmax)      
    real :: gustX0                ! Starting position of the gust upstream (# rotor radii)

End MODULE iecgust
