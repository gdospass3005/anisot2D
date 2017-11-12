PROGRAM model2D

  USE mdle_model2D
  IMPLICIT NONE

  INTEGER :: Nx,Nz,nmedia
  REAL    :: dx,dz
  REAL, DIMENSION(:,:), ALLOCATABLE :: x,z
  REAL, DIMENSION(:,:), ALLOCATABLE :: C11
  REAL, DIMENSION(:,:), ALLOCATABLE :: C13
  REAL, DIMENSION(:,:), ALLOCATABLE :: C33
  REAL, DIMENSION(:,:), ALLOCATABLE :: C44
  REAL, DIMENSION(:,:), ALLOCATABLE :: rhox
  REAL, DIMENSION(:), ALLOCATABLE :: value1
  REAL, DIMENSION(:), ALLOCATABLE :: value2
  REAL, DIMENSION(:), ALLOCATABLE :: value3
  REAL, DIMENSION(:), ALLOCATABLE :: value4
  REAL, DIMENSION(:), ALLOCATABLE :: value5
  CHARACTER(LEN=20) :: infile='data.dat'
  CHARACTER(LEN=20) :: outfile1='C11.ad'
  CHARACTER(LEN=20) :: outfile2='C13.ad'
  CHARACTER(LEN=20) :: outfile3='C33.ad'
  CHARACTER(LEN=20) :: outfile4='C44.ad'
  CHARACTER(LEN=20) :: outfile5='rhox.ad'
  INTEGER :: i
  CHARACTER(LEN=20) :: layerno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! INPUT

  OPEN(25,FILE=infile,STATUS='UNKNOWN')
  READ(25,'(t10,i10)') Nx
  READ(25,'(t10,i10)') Nz
  READ(25,'(t10,f10.4)') dx
  READ(25,'(t10,f10.4)') dz
  READ(25,'(t10,i10)') nmedia
  ALLOCATE(x(4,nmedia),z(4,nmedia))
  ALLOCATE(value1(nmedia))
  ALLOCATE(value2(nmedia))
  ALLOCATE(value3(nmedia))
  ALLOCATE(value4(nmedia))
  ALLOCATE(value5(nmedia))
  DO i=1,nmedia
     READ(25,*) layerno
     READ(25,'(t10,f10.4,f10.4,f10.4,f10.4,f10.4)') value1(i),&
          value2(i),value3(i),value4(i),value5(i)
     READ(25,'(t10,f10.4,f10.4)') x(1,i),z(1,i)
     READ(25,'(t10,f10.4,f10.4)') x(2,i),z(2,i)
     READ(25,'(t10,f10.4,f10.4)') x(3,i),z(3,i)
     READ(25,'(t10,f10.4,f10.4)') x(4,i),z(4,i)
     PRINT*, layerno
     PRINT*, value1(i),value2(i),value3(i),value4(i),value5(i)
     PRINT*, x(1,i),z(1,i)
     PRINT*, x(2,i),z(2,i)
     PRINT*, x(3,i),z(3,i)
     PRINT*, x(4,i),z(4,i)
  END DO

  PRINT*,'PROGRAMA GERADOR DE PAINEIS PARA A MODELAGEM POROELASTICA'

  ALLOCATE(C11(Nx,Nz))
  ALLOCATE(C13(Nx,Nz))
  ALLOCATE(C33(Nx,Nz))
  ALLOCATE(C44(Nx,Nz))
  ALLOCATE(rhox(Nx,Nz))
  CALL modele2D(Nx,Nz,nmedia,dx,dz,x,z,value1,value2,value3,&
       value4,value5,C11,C13,C33,C44,rhox)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! OUTPUT

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

END PROGRAM model2D




