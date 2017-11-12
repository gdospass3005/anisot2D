MODULE mdle_source

IMPLICIT NONE

CONTAINS
SUBROUTINE source_init(trs,dt,ns,fpeak,tdelay)
  IMPLICIT NONE
  REAL, INTENT(INOUT), DIMENSION(ns) :: trs
  INTEGER, INTENT(IN) :: ns
  REAL, INTENT(IN) :: dt, fpeak
  REAL, INTENT(OUT) :: tdelay
  INTEGER :: i
  REAL :: t,pi=3.141592653589793238462643383279502884197

  ! fpeak = fcut/(3*SQRT(pi))
  tdelay = 4.0 /(3.*fpeak*SQRT(2.*pi))

  do i=1,ns
     t=(i-1)*dt - tdelay
     trs(i) = (1 - 2*pi*(pi*fpeak*(t-tdelay))**2)*exp(-pi*(pi*fpeak*(t-tdelay))**2)
  end do

END SUBROUTINE source_init

END MODULE
