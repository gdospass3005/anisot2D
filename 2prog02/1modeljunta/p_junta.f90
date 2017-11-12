PROGRAM p_junta

  IMPLICIT NONE

  REAL, DIMENSION(:,:), ALLOCATABLE :: C11
  REAL, DIMENSION(:,:), ALLOCATABLE :: C13
  REAL, DIMENSION(:,:), ALLOCATABLE :: C33
  REAL, DIMENSION(:,:), ALLOCATABLE :: C44
  REAL, DIMENSION(:,:), ALLOCATABLE :: rhox
  REAL, DIMENSION(:,:), ALLOCATABLE :: aux
  CHARACTER(LEN=20) :: outfile1='C11.ad'
  CHARACTER(LEN=20) :: outfile2='C13.ad'
  CHARACTER(LEN=20) :: outfile3='C33.ad'
  CHARACTER(LEN=20) :: outfile4='C44.ad'
  CHARACTER(LEN=20) :: outfile5='rhox.ad'
  CHARACTER(LEN=20) :: infile='data.dat'
  CHARACTER(LEN=20) :: layerno
  REAL, DIMENSION(:), ALLOCATABLE :: value1
  REAL, DIMENSION(:), ALLOCATABLE :: value2
  REAL, DIMENSION(:), ALLOCATABLE :: value3
  REAL, DIMENSION(:), ALLOCATABLE :: value4
  REAL, DIMENSION(:), ALLOCATABLE :: value5
  INTEGER :: i,j,Nx,Nz,nmedia
  REAL :: dx,dz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT

  OPEN(25,FILE=infile,STATUS='UNKNOWN')
  READ(25,'(t10,i10)') Nx
  READ(25,'(t10,i10)') Nz
  READ(25,'(t10,f10.4)') dx
  READ(25,'(t10,f10.4)') dz
  READ(25,'(t10,i10)') nmedia
  ALLOCATE(value1(nmedia))
  ALLOCATE(value2(nmedia))
  ALLOCATE(value3(nmedia))
  ALLOCATE(value4(nmedia))
  ALLOCATE(value5(nmedia))

  DO i=1,nmedia
     READ(25,*) layerno
     READ(25,'(t10,f10.4,f10.4,f10.4,f10.4,f10.4)') value1(i),&
          value3(i),value4(i),value2(i),value5(i)
  END DO

  PRINT*,'PROGRAMA GERADOR DE PAINEIS PARA A MODELAGEM SISMICA'

  ALLOCATE(C11(Nx,Nz))
  ALLOCATE(C13(Nx,Nz))
  ALLOCATE(C33(Nx,Nz))
  ALLOCATE(C44(Nx,Nz))
  ALLOCATE(rhox(Nx,Nz))

  C11=0.;C13=0.;C33=0.;C44=0.;rhox=0.

  ALLOCATE(aux(Nx,Nz))

  OPEN(30,FILE="m_saltwater.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(1)
           C13(i,j) = value2(1)
           C33(i,j) = value3(1)
           C44(i,j) = value4(1)
           rhox(i,j) = value5(1)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_shale1.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(2)
           C13(i,j) = value2(2)
           C33(i,j) = value3(2)
           C44(i,j) = value4(2)
           rhox(i,j) = value5(2)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_shale2.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(3)
           C13(i,j) = value2(3)
           C33(i,j) = value3(3)
           C44(i,j) = value4(3)
           rhox(i,j) = value5(3)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_shale3.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(4)
           C13(i,j) = value2(4)
           C33(i,j) = value3(4)
           C44(i,j) = value4(4)
           rhox(i,j) = value5(4)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_bottomsand.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(5)
           C13(i,j) = value2(5)
           C33(i,j) = value3(5)
           C44(i,j) = value4(5)
           rhox(i,j) = value5(5)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_brine_sand.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(6)
           C13(i,j) = value2(6)
           C33(i,j) = value3(6)
           C44(i,j) = value4(6)
           rhox(i,j) = value5(6)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_gas_sand.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(7)
           C13(i,j) = value2(7)
           C33(i,j) = value3(7)
           C44(i,j) = value4(7)
           rhox(i,j) = value5(7)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_oil_sand.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(8)
           C13(i,j) = value2(8)
           C33(i,j) = value3(8)
           C44(i,j) = value4(8)
           rhox(i,j) = value5(8)
        END IF
     END DO
  END DO

  OPEN(30,FILE="m_salt.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx
     READ(30, REC=i) aux(i,:)
  END DO
  CLOSE(30)
  DO i=1,Nx
     DO j=1,Nz
        IF (aux(i,j) == 1.) THEN
           C11(i,j) = value1(9)
           C13(i,j) = value2(9)
           C33(i,j) = value3(9)
           C44(i,j) = value4(9)
           rhox(i,j) = value5(9)
        END IF
     END DO
  END DO

  OPEN(35,FILE=outfile1,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO I=1,Nx
     WRITE(35, REC=i) C11(i,:) 
  END DO
  OPEN(36,FILE=outfile2,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO I=1,Nx
     WRITE(36, REC=i) C13(i,:) 
  END DO
  OPEN(37,FILE=outfile3,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO I=1,Nx
     WRITE(37, REC=i) C33(i,:) 
  END DO
  OPEN(38,FILE=outfile4,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO I=1,Nx
     WRITE(38, REC=i) C44(i,:) 
  END DO
  OPEN(39,FILE=outfile5,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  DO I=1,Nx
     WRITE(39, REC=i) rhox(i,:) 
  END DO


END PROGRAM p_junta
