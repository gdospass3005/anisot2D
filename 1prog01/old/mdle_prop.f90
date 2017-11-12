MODULE mdle_prop

CONTAINS
  SUBROUTINE prop_anisot(Vx,Vz,Sxx,Sxz,Szz,&
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

    ! Space derivatives of velocities => stresses
    DO j=3,Nz-2
       DO i=3,Nx-2

          aux1 = fr1*(Vx(i+2,j)-Vx(i-1,j)) + &
               fr2*(Vx(i+1,j)-Vx(i,j)); aux3=aux1
          aux1 = C11(i,j)*aux1/dx
          aux2 = fr1*(Vz(i,j+1)-Vz(i,j-2)) + &
               fr2*(Vz(i,j)-Vz(i,j-1)); aux4=aux2
          aux2 = C13(i,j)*aux2*ldx
          aux = dt*(aux1+aux2)
          Sxx(i,j) = aux + Sxx(i,j)

          aux1 = C13(i,j)*aux3*ldx
          aux2 = C33(i,j)*aux4*ldx
          aux = dt*(aux1+aux2)
          Szz(i,j) = aux + Szz(i,j)

          aux1 = fr1*(Vz(i+1,j)-Vz(i-2,j)) + &
               fr2*(Vz(i,j)-Vz(i-1,j))
          aux1 = C44(i,j)*aux1/dx
          aux2 = fr1*(Vx(i,j+2)-Vx(i,j-1)) + &
               fr2*(Vx(i,j+1)-Vx(i,j))
          aux2 = C44(i,j)*aux2*ldx
          aux = dt*(aux1+aux2)
          Sxz(i,j) = aux + Sxz(i,j)

       END DO
    END DO

    ! Space derivatives of stresses => velocities
    DO j=3,Nz-2
       DO i=3,Nx-2

          aux1 = fr1*(Sxx(i+1,j)-Sxx(i-2,j)) + &
               fr2*(Sxx(i,j)-Sxx(i-1,j))
          aux1 = aux1/dx
          aux2 = fr1*(Sxz(i,j+1)-Sxz(i,j-2)) + &
               fr2*(Sxz(i,j)-Sxz(i,j-1))
          aux2 = aux2*ldx
          aux = dt*rhox(i,j)*(aux1+aux2)
          Vx(i,j) = aux + Vx(i,j)

          aux1 = fr1*(Sxz(i+2,j)-Sxz(i-1,j)) + &
               fr2*(Sxz(i+1,j)-Sxz(i,j))
          aux1 = aux1*ldx
          aux2 = fr1*(Szz(i,j+2)-Szz(i,j-1)) + &
               fr2*(Szz(i,j+1)-Szz(i,j))
          aux2 = aux2*ldx
          aux = dt*rhoz(i,j)*(aux1+aux2)
          Vz(i,j) = aux + Vz(i,j)

       END DO
    END DO


  END SUBROUTINE prop_anisot

END MODULE mdle_prop
