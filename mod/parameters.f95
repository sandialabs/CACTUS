module parameters

    ! Sizes for arrays

    integer :: MaxBlades, MaxSegPerBlade, MaxSegEndPerBlade, MaxSegEnds, MaxSeg     ! Maximums for blade geometry arrays
    integer :: MaxStruts                                                            ! Maximums for strut arrays
    integer :: MaxAirfoilSect, MaxReVals, MaxAOAVals                ! Maximums for airfoil coeff data arrays
    integer :: MaxRevs                              ! Max revolutions
    integer :: MaxTimeStepPerRev                            ! Max time increments in a revolution
    integer :: MaxTimeSteps                                                         ! Max total time increments (MaxRevs * MaxTimeStepPerRev)
    integer :: MaxWakeNodes                             ! Max number of points in wake for an element (MaxRevs * MaxTimeStepPerRev)
    integer :: MaxNLIters                               ! Max iterations to perform in converging non linear component of blade element system

end module parameters
