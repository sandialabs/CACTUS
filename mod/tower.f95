MODULE tower

  implicit none

  ! user inputs
  real :: tower_D, tower_x, tower_ybot, tower_ytop
  integer :: Itower, tower_Npts

  real :: tower_CD, theta0
  real, parameter :: Delta_wake = 0.289, W0_wake = 1.75
  real, parameter :: S_wake = 0.083, alpha_wake = 0.693

  real, allocatable, dimension(:) :: tower_y, tower_Vx

CONTAINS

  SUBROUTINE setup_tower

    implicit none
    
    integer :: i

    Allocate(tower_y(tower_Npts), tower_Vx(tower_Npts))
    

    If (tower_Npts .LE. 1) Then
       Write(*,'(''Inadequate number of tower points.'')')
    End If

    Do i = 1,tower_Npts
       tower_y(i) = tower_ybot + (REAL(i)-1) / &
            (REAL(tower_Npts) - 1) * (tower_ytop-tower_ybot)
       tower_Vx(i) = 0.0
    End Do

    theta0 = 0.5 * tower_CD * tower_D

  END SUBROUTINE setup_tower


  FUNCTION interp_tower_Vx(y)

    implicit none

    real :: interp_tower_Vx, y
    integer :: i

    interp_tower_Vx = 0.0
    Do i = 1,tower_Npts-1
       If ((y .GE. tower_y(i)) .AND. (y .LE. tower_y(i+1))) Then
          interp_tower_Vx = tower_Vx(i) + (y - tower_y(i)) &
               / (tower_y(i+1) - tower_y(i)) * (tower_Vx(i+1) - &
               tower_Vx(i))
          Exit
       End If
    End Do

    Return

  END FUNCTION interp_tower_Vx


  FUNCTION wake_defect_velocity(x,y,z)

    implicit none
    
    real :: x, y, z, wake_defect_velocity
    real :: W, xi, Vx_tow, u0, f_wake

    If (x .LE. tower_x) Then
       wake_defect_velocity = 0.0
       Return
    End If

    If ((y .LT. tower_ybot) .OR. (y .GT. tower_ytop)) Then
       wake_defect_velocity = 0.0
       Return
    End If


    ! Calculate the wake half-width at this x location
    W = Delta_wake * sqrt((x-tower_x)*theta0)
    xi = z / W ! non-dimensional transverse coordinate
    Vx_tow = interp_tower_Vx(y)  ! x-dir fluid velocity at tower
    !Write(20,'(F20.12)') Vx_tow
    u0 = W0_wake * sqrt(theta0/(x-tower_x )) * Vx_tow ! max velocity
    !  defect at this x location
    f_wake = exp(-alpha_wake * xi * xi) ! non-dimensional velocity defect    
    wake_defect_velocity = u0 * f_wake ! dimensional velocity defect

    if (wake_defect_velocity .GT. 0.9 * Vx_tow) then
       wake_defect_velocity = 0.9 * Vx_tow
    end if
    
    Return

  END FUNCTION wake_defect_velocity

END MODULE tower
