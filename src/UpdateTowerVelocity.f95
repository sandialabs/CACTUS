SUBROUTINE UpdateTowerVelocity()
  
  use tower
  use configr

  implicit none

  real :: uFSt,vFSt,wFSt,utow,vtow,wtow
  integer :: i,ygcerr


  Do i = 1,tower_Npts

     ! Freestream velocity at tower location
     Call CalcFreestream(tower_x,tower_y(i),0.0,uFSt,vFSt,wFSt,ygcerr)  

     ! Induced velocity at tower location
     Call CalcIndVel(NT,ntTerm,NBE,NB,NE,tower_x,tower_y(i),0.0,utow,vtow,wtow)

     tower_Vx(i) = uFSt + utow

  End Do


  Return

End SUBROUTINE UpdateTowerVelocity
