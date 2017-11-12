MODULE mdle_taper

CONTAINS

SUBROUTINE bt_exp_create(taper,nb,F)
IMPLICIT NONE
  REAL :: F
  INTEGER :: nb
  REAL, DIMENSION(nb) :: taper
  INTEGER :: i

  DO i=1,nb
    taper(i) = exp( -(F*(REAL(nb) - REAL(i) ))**2 )
!    PRINT*, taper(i)
  END DO

END SUBROUTINE bt_exp_create


SUBROUTINE bt_apply(pp,Nx,Nz,nb,taper)
IMPLICIT NONE
  INTEGER :: NX,NZ,nb
  REAL, DIMENSION(nb) :: taper
  REAL, DIMENSION(NX,NZ) :: pp
  INTEGER :: i

  DO i=1,NZ
     pp(1:nb,i) = pp(1:nb,i) * taper
     pp(NX:NX-nb+1:-1,i) =  pp(NX:NX-nb+1:-1,i) * taper
  END DO

  DO i=1,NX
    pp(i,1:nb:1) = pp(i,1:nb:1) * taper
    pp(i,NZ:NZ-nb+1:-1) =  pp(i,NZ:NZ-nb+1:-1) * taper  
  END DO

END SUBROUTINE bt_apply


SUBROUTINE bt_apply_left(pp,Nx,Nz,nb,taper)
IMPLICIT NONE
  INTEGER :: NX,NZ,nb
  REAL, DIMENSION(nb) :: taper
  REAL, DIMENSION(NX,NZ) :: pp
  INTEGER :: i

  DO i=1,NZ
     pp(1:nb,i) = pp(1:nb,i) * taper
  END DO

END SUBROUTINE bt_apply_left


SUBROUTINE bt_apply_right(pp,Nx,Nz,nb,taper)
IMPLICIT NONE
  INTEGER :: NX,NZ,nb
  REAL, DIMENSION(nb) :: taper
  REAL, DIMENSION(NX,NZ) :: pp
  INTEGER :: i

  DO i=1,NZ
     pp(NX:NX-nb+1:-1,i) =  pp(NX:NX-nb+1:-1,i) * taper
  END DO

END SUBROUTINE bt_apply_right


SUBROUTINE bt_apply_top(pp,Nx,Nz,nb,taper)
IMPLICIT NONE
  INTEGER :: NX,NZ,nb
  REAL, DIMENSION(nb) :: taper
  REAL, DIMENSION(NX,NZ) :: pp
  INTEGER :: i

  DO i=1,NX
    pp(i,1:nb:1) = pp(i,1:nb:1) * taper
  END DO

END SUBROUTINE bt_apply_top


SUBROUTINE bt_apply_bottom(pp,Nx,Nz,nb,taper)
IMPLICIT NONE
  INTEGER :: NX,NZ,nb
  REAL, DIMENSION(nb) :: taper
  REAL, DIMENSION(NX,NZ) :: pp
  INTEGER :: i

  DO i=1,NX
    pp(i,NZ:NZ-nb+1:-1) =  pp(i,NZ:NZ-nb+1:-1) * taper  
  END DO

END SUBROUTINE bt_apply_bottom


END MODULE mdle_taper



