MODULE mdle_model2D

CONTAINS

  SUBROUTINE modele2D(Nx,Nz,nmedia,dx,dz,x,z,value1,value2,value3,&
       value4,value5,C11,C13,C33,C44,rhox)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx,Nz,nmedia
    REAL, INTENT(IN) :: dx,dz
    REAL, DIMENSION(4,nmedia), INTENT(IN) :: x,z
    REAL, DIMENSION(nmedia), INTENT(IN) :: value1,value2,value3,&
       value4,value5
    REAL, DIMENSION(Nx,Nz), INTENT(OUT) :: C11,C13,C33,C44,rhox
    INTEGER :: imedia
    INTEGER :: ix,iz
    REAL    :: xpos,zpos
    REAL    :: xallow1,xallow2
    REAL    :: zallow1,zallow2
    REAL    :: alpha1,alpha2
    REAL    :: beta1,beta2

    C11=0.0
    DO imedia=1,nmedia
       DO ix=1,Nx
          DO iz=1,Nz
             xpos = (ix-1)*dx
             zpos = (iz-1)*dz
             alpha1 = (z(4,imedia)-z(1,imedia))/(x(4,imedia)-x(1,imedia))
             beta1 = 1/((z(2,imedia)-z(1,imedia))/(x(2,imedia)-x(1,imedia)))
             alpha2 = (z(2,imedia)-z(3,imedia))/(x(2,imedia)-x(3,imedia))
             beta2 = 1/((z(4,imedia)-z(3,imedia))/(x(4,imedia)-x(3,imedia)))
             xallow1 = x(1,imedia) + (zpos-z(1,imedia))*beta1
             xallow2 = x(3,imedia) + (zpos-z(3,imedia))*beta2
             zallow1 = z(1,imedia) + (xpos-x(1,imedia))*alpha1
             zallow2 = z(3,imedia) + (xpos-x(3,imedia))*alpha2
             IF ((xpos >= xallow1).and.(xpos <= xallow2).and.&
                  (zpos >= zallow1).and.(zpos <= zallow2))  &
                  C11(ix,iz) = value1(imedia)
          END DO
       END DO
    END DO

    C13=0.0
    DO imedia=1,nmedia
       DO ix=1,Nx
          DO iz=1,Nz
             xpos = (ix-1)*dx
             zpos = (iz-1)*dz
             alpha1 = (z(4,imedia)-z(1,imedia))/(x(4,imedia)-x(1,imedia))
             beta1 = 1/((z(2,imedia)-z(1,imedia))/(x(2,imedia)-x(1,imedia)))
             alpha2 = (z(2,imedia)-z(3,imedia))/(x(2,imedia)-x(3,imedia))
             beta2 = 1/((z(4,imedia)-z(3,imedia))/(x(4,imedia)-x(3,imedia)))
             xallow1 = x(1,imedia) + (zpos-z(1,imedia))*beta1
             xallow2 = x(3,imedia) + (zpos-z(3,imedia))*beta2
             zallow1 = z(1,imedia) + (xpos-x(1,imedia))*alpha1
             zallow2 = z(3,imedia) + (xpos-x(3,imedia))*alpha2
             IF ((xpos >= xallow1).and.(xpos <= xallow2).and.&
                  (zpos >= zallow1).and.(zpos <= zallow2))  &
                  C13(ix,iz) = value2(imedia)
          END DO
       END DO
    END DO

    C33=0.0
    DO imedia=1,nmedia
       DO ix=1,Nx
          DO iz=1,Nz
             xpos = (ix-1)*dx
             zpos = (iz-1)*dz
             alpha1 = (z(4,imedia)-z(1,imedia))/(x(4,imedia)-x(1,imedia))
             beta1 = 1/((z(2,imedia)-z(1,imedia))/(x(2,imedia)-x(1,imedia)))
             alpha2 = (z(2,imedia)-z(3,imedia))/(x(2,imedia)-x(3,imedia))
             beta2 = 1/((z(4,imedia)-z(3,imedia))/(x(4,imedia)-x(3,imedia)))
             xallow1 = x(1,imedia) + (zpos-z(1,imedia))*beta1
             xallow2 = x(3,imedia) + (zpos-z(3,imedia))*beta2
             zallow1 = z(1,imedia) + (xpos-x(1,imedia))*alpha1
             zallow2 = z(3,imedia) + (xpos-x(3,imedia))*alpha2
             IF ((xpos >= xallow1).and.(xpos <= xallow2).and.&
                  (zpos >= zallow1).and.(zpos <= zallow2))  &
                  C33(ix,iz) = value3(imedia)
          END DO
       END DO
    END DO

    C44=0.0
    DO imedia=1,nmedia
       DO ix=1,Nx
          DO iz=1,Nz
             xpos = (ix-1)*dx
             zpos = (iz-1)*dz
             alpha1 = (z(4,imedia)-z(1,imedia))/(x(4,imedia)-x(1,imedia))
             beta1 = 1/((z(2,imedia)-z(1,imedia))/(x(2,imedia)-x(1,imedia)))
             alpha2 = (z(2,imedia)-z(3,imedia))/(x(2,imedia)-x(3,imedia))
             beta2 = 1/((z(4,imedia)-z(3,imedia))/(x(4,imedia)-x(3,imedia)))
             xallow1 = x(1,imedia) + (zpos-z(1,imedia))*beta1
             xallow2 = x(3,imedia) + (zpos-z(3,imedia))*beta2
             zallow1 = z(1,imedia) + (xpos-x(1,imedia))*alpha1
             zallow2 = z(3,imedia) + (xpos-x(3,imedia))*alpha2
             IF ((xpos >= xallow1).and.(xpos <= xallow2).and.&
                  (zpos >= zallow1).and.(zpos <= zallow2))  &
                  C44(ix,iz) = value4(imedia)
          END DO
       END DO
    END DO

    rhox=0.0
    DO imedia=1,nmedia
       DO ix=1,Nx
          DO iz=1,Nz
             xpos = (ix-1)*dx
             zpos = (iz-1)*dz
             alpha1 = (z(4,imedia)-z(1,imedia))/(x(4,imedia)-x(1,imedia))
             beta1 = 1/((z(2,imedia)-z(1,imedia))/(x(2,imedia)-x(1,imedia)))
             alpha2 = (z(2,imedia)-z(3,imedia))/(x(2,imedia)-x(3,imedia))
             beta2 = 1/((z(4,imedia)-z(3,imedia))/(x(4,imedia)-x(3,imedia)))
             xallow1 = x(1,imedia) + (zpos-z(1,imedia))*beta1
             xallow2 = x(3,imedia) + (zpos-z(3,imedia))*beta2
             zallow1 = z(1,imedia) + (xpos-x(1,imedia))*alpha1
             zallow2 = z(3,imedia) + (xpos-x(3,imedia))*alpha2
             IF ((xpos >= xallow1).and.(xpos <= xallow2).and.&
                  (zpos >= zallow1).and.(zpos <= zallow2))  &
                  rhox(ix,iz) = value5(imedia)
          END DO
       END DO
    END DO


  END SUBROUTINE modele2D

END MODULE mdle_model2D




