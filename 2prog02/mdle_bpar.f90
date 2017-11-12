MODULE mdle_bpar

IMPLICIT NONE

CONTAINS 
  SUBROUTINE bpar_left(nxp,acp,pvp,oVel,Nzm,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nzm, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(1:4,Nzm:Nz) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(1:4,Nzm:Nz) :: pvp
    REAL, INTENT(INOUT),  DIMENSION(1:2,Nzm:Nz) :: nxp
    INTEGER :: iz,ix

    ! Borda paraxial esquerda
    DO ix=1,2
    DO iz=Nzm,Nz
       nxp(ix,iz) = (oVel(ix,iz)) * (dt/dx) *&
         (acp(ix+1,iz)-acp(ix,iz)-(pvp(ix+2,iz)-pvp(ix+1,iz))) +&
         acp(ix,iz) + acp(ix+1,iz) - pvp(ix+1,iz)
     END DO
    END DO
  END SUBROUTINE bpar_left


  SUBROUTINE bpar_right(nxp,acp,pvp,oVel,Nx,Nzm,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nzm, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(Nx-3:Nx,Nzm:Nz) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(Nx-3:Nx,Nzm:Nz) :: pvp
    REAL, INTENT(INOUT),  DIMENSION(Nx-1:Nx,Nzm:Nz) :: nxp
    INTEGER :: iz,ix

    ! Borda paraxial direita
    DO ix=Nx-1,Nx
    DO iz=Nzm,Nz
       nxp(ix,iz) = -(oVel(ix,iz)) * (dt/dx) *&
         (acp(ix,iz)-acp(ix-1,iz)-(pvp(ix-1,iz)-pvp(ix-2,iz))) +&
         acp(ix,iz) + acp(ix-1,iz) - pvp(ix-1,iz)
    END DO
    END DO
  END SUBROUTINE bpar_right

    
  SUBROUTINE bpar_bottom(nxp,acp,pvp,oVel,Nx,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(1:Nx,Nz-3:Nz) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(1:Nx,Nz-3:Nz) :: pvp
    REAL, INTENT(INOUT),  DIMENSION(1:Nx,Nz-1:Nz) :: nxp
    INTEGER :: ix,iz

    ! Borda paraxial de baixo
    DO ix=2,(Nx-1)
    DO iz=Nz-1,Nz
       nxp(ix,iz) = -(oVel(ix,iz)) * (dt/dx) *&
         (acp(ix,iz)-acp(ix,iz-1)-(pvp(ix,iz-1)-pvp(ix,iz-2))) +&
         acp(ix,iz) + acp(ix,iz-1) - pvp(ix,iz-1)
    END DO
    END DO
  END SUBROUTINE bpar_bottom


  SUBROUTINE bpar_top(nxp,acp,pvp,oVel,Nx,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(1:Nx,1:4) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(1:Nx,1:4) :: pvp
    REAL, INTENT(INOUT),  DIMENSION(1:Nx,1:2) :: nxp
    INTEGER :: ix,iz
    
    ! Borda paraxial de cima
    DO ix=2,(Nx-1)
       DO iz=1,2
          nxp(ix,iz) =  (oVel(ix,iz)) * (dt/dx) *&
        (acp(ix,iz+1)-acp(ix,iz)-(pvp(ix,iz+2)-pvp(ix,iz+1))) +&
        acp(ix,iz) + acp(ix,iz+1) - pvp(ix,iz+1)
       END DO
    END DO
  END SUBROUTINE bpar_top


END MODULE



