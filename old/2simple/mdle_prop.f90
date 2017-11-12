MODULE mdle_prop

IMPLICIT NONE

CONTAINS 
  SUBROUTINE prop_wave(nxp,acp,pvp,pVel,maxVel,&
       Nxm,Nx,Nzm,Nz,dx,dz,dt,t,lx,lz)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nxm, Nx, Nzm, Nz, lx, lz
    REAL, INTENT(IN) :: dx, dz, dt, maxVel, t
    REAL, INTENT(IN),  DIMENSION(:,:) :: acp,pVel
    REAL, INTENT(INOUT),  DIMENSION(:,:) :: pvp
    REAL, INTENT(OUT), DIMENSION(:,:) :: nxp
    REAL :: sc1,sc2,tfast
    INTEGER :: ix,iz

    sc1=(dt * dt)/(dx * dx)
    sc2=0.25*sc1
    tfast=(maxVel*t)**2

    DO ix=Nxm+1,Nx-1
       DO iz=Nzm+1,Nz-1
          IF (( ((ix-lx)*dx)**2+((iz-lz)*dz)**2 ) <= tfast ) THEN
             nxp(ix,iz) = sc1 * &
                  (acp(ix+1,iz) + acp(ix-1,iz) + acp(ix,iz+1) + acp(ix,iz-1)&
                  - 4.0*acp(ix,iz))

             nxp(ix,iz) = pVel(ix,iz)*nxp(ix,iz) +&
                  2.0 *acp(ix,iz)  - pvp(ix,iz)
          END IF

       END DO
    END DO


  END SUBROUTINE prop_wave

END MODULE
