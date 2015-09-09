module vortex

    ! Vortex core model parameters for trailing and spanwise wake elements

    integer :: ivtxcor              ! 0 for no core model, 1 for constant velocity core model, 2 for linear velocity core model
    real :: vRad_B          ! bound vortex core radius
    real :: vRad2_B         ! bound vortex core radius squared
    real :: vRad_T      ! trailing wake core radius
    real :: vRad2_T     ! trailing wake core radius squared
    real :: vRad_S          ! spanwise wake core radius
    real :: vRad2_S         ! spanwise wake core radius squared
    real :: vCutOffRad       ! cut-off radius (for bound/trailing/spanwise cores)


end module vortex
