MODULE mdle_bpar

IMPLICIT NONE

CONTAINS 
  SUBROUTINE bpar_left(nxp,acp,pvp,oVel,Nx,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(:,:) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(:,:) :: pvp
    REAL, INTENT(OUT), DIMENSION(:,:) :: nxp
    INTEGER :: ix,iz

     ! Borda paraxial esquerda
    ix=1
    DO iz=2,(Nz-1)
       nxp(ix,iz) = SQRT(oVel(ix,iz)) * (dt/dx) *&
         (acp(ix+1,iz)-acp(ix,iz)-(pvp(ix+2,iz)-pvp(ix+1,iz))) +&
         acp(ix,iz) + acp(ix+1,iz) - pvp(ix+1,iz)
    END DO
  END SUBROUTINE bpar_left

    
  SUBROUTINE bpar_right(nxp,acp,pvp,oVel,Nx,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(:,:) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(:,:) :: pvp
    REAL, INTENT(OUT), DIMENSION(:,:) :: nxp
    INTEGER :: ix,iz

    ! Borda paraxial direita
    ix=Nx
    DO iz=2,(Nz-1)
       nxp(ix,iz) = -SQRT(oVel(ix,iz)) * (dt/dx) *&
         (acp(ix,iz)-acp(ix-1,iz)-(pvp(ix-1,iz)-pvp(ix-2,iz))) +&
         acp(ix,iz) + acp(ix-1,iz) - pvp(ix-1,iz)
    END DO
  END SUBROUTINE bpar_right

    
  SUBROUTINE bpar_top(nxp,acp,pvp,oVel,Nx,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(:,:) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(:,:) :: pvp
    REAL, INTENT(OUT), DIMENSION(:,:) :: nxp
    INTEGER :: ix,iz
    
    ! Borda paraxial de cima
    iz=1
    DO ix=2,(Nx-1)
    
      nxp(ix,iz) = SQRT(oVel(ix,iz)) * (dt/dx) *&
        (acp(ix,iz+1)-acp(ix,iz)-(pvp(ix,iz+2)-pvp(ix,iz+1))) +&
        acp(ix,iz) + acp(ix,iz+1) - pvp(ix,iz+1)
    !   nxp(ix,iz) = 0 ! Free surface
    END DO
  END SUBROUTINE bpar_top

    
  SUBROUTINE bpar_bottom(nxp,acp,pvp,oVel,Nx,Nz,dx,dz,dt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx, Nz
    REAL, INTENT(IN) :: dx, dz, dt
    REAL, INTENT(IN),  DIMENSION(:,:) :: acp,oVel
    REAL, INTENT(INOUT),  DIMENSION(:,:) :: pvp
    REAL, INTENT(OUT), DIMENSION(:,:) :: nxp
    INTEGER :: ix,iz

    ! Borda paraxial de baixo
    iz=Nz
    DO ix=2,(Nx-1)
       nxp(ix,iz) = -SQRT(oVel(ix,iz)) * (dt/dx) *&
         (acp(ix,iz)-acp(ix,iz-1)-(pvp(ix,iz-1)-pvp(ix,iz-2))) +&
         acp(ix,iz) + acp(ix,iz-1) - pvp(ix,iz-1)
    END DO
  END SUBROUTINE bpar_bottom
END MODULE



