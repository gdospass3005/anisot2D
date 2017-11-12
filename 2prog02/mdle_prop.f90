MODULE mdle_prop

CONTAINS

  SUBROUTINE prop_anisot_2D(Vx,Vz,Sxx,Sxz,Szz,&
       C11,C13,C33,C44,rhox,rhoz,dx,dt,Nx,Nz)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: Vx,Vz,Sxx,Sxz,Szz
    REAL, DIMENSION(:,:), INTENT(IN)    :: C11,C13,C33,C44
    REAL, DIMENSION(:,:), INTENT(IN)    :: rhox,rhoz
    REAL,    INTENT(IN) :: dx,dt
    INTEGER, INTENT(IN) :: Nx,Nz

    REAL :: aux1,aux2,aux3,aux4,aux
    INTEGER :: i,j
    REAL :: fr1,fr2,ldx

    fr1 = -1./24; fr2=9./8; ldx=1./dx

    ! Space derivatives of stresses => velocities
    DO j=3,Nz-2
       DO i=3,Nx-2

          aux1 = (fr1*(Sxx(i+1,j)-Sxx(i-2,j)) + &
               fr2*(Sxx(i,j)-Sxx(i-1,j))) * ldx
          Vx(i,j) = Vx(i,j) + (dt * rhox(i,j) * aux1)

          aux1 = (fr1*(Sxz(i,j+1)-Sxz(i,j-2)) + &
               fr2*(Sxz(i,j)-Sxz(i,j-1))) * ldx
          Vx(i,j) = Vx(i,j) + (dt * rhox(i,j) * aux1)

          aux1 = (fr1*(Sxz(i+2,j)-Sxz(i-1,j)) + &
               fr2*(Sxz(i+1,j)-Sxz(i,j))) * ldx
          Vz(i,j) = Vz(i,j) + (dt * rhox(i,j) * aux1)

          aux1 = (fr1*(Szz(i,j+2)-Szz(i,j-1)) + &
               fr2*(Szz(i,j+1)-Szz(i,j))) * ldx
          Vz(i,j) = Vz(i,j) + (dt * rhox(i,j) * aux1)

       END DO
    END DO

    ! Space derivatives of velocities => stresses
    DO j=3,Nz-2
       DO i=3,Nx-2

          aux1 = (fr1*(Vx(i+2,j)-Vx(i-1,j)) + &
               fr2*(Vx(i+1,j)-Vx(i,j))) * ldx
          Sxx(i,j) = Sxx(i,j) + (dt * C11(i,j) * aux1)
          Szz(i,j) = Szz(i,j) + (dt * C13(i,j) * aux1)

          aux1 = (fr1*(Vz(i,j+1)-Vz(i,j-2)) + &
               fr2*(Vz(i,j)-Vz(i,j-1))) * ldx
          Sxx(i,j) = Sxx(i,j) + (dt * C13(i,j) * aux1)
          Szz(i,j) = Szz(i,j) + (dt * C33(i,j) * aux1)

          aux1 = (fr1*(Vz(i+1,j)-Vz(i-2,j)) + &
               fr2*(Vz(i,j)-Vz(i-1,j))) * ldx
          Sxz(i,j) = Sxz(i,j) + (dt * C44(i,j) * aux1)

          aux1 = (fr1*(Vx(i,j+2)-Vx(i,j-1)) + &
               fr2*(Vx(i,j+1)-Vx(i,j))) * ldx
          Sxz(i,j) = Sxz(i,j) + (dt * C44(i,j) * aux1)

       END DO
    END DO

  END SUBROUTINE prop_anisot_2D

END MODULE mdle_prop
