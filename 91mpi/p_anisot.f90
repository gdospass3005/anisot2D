PROGRAM anisot

  USE mdle_io_utils                                                  
  USE mdle_source
  USE mdle_prop
  USE mdle_taper

  IMPLICIT NONE

  INCLUDE 'mpif.h'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: rank, size, ierr, comm, np
  INTEGER :: proot,status(MPI_STATUS_SIZE)
  INTEGER :: source, dest, tag
  INTEGER, DIMENSION(:), ALLOCATABLE :: div_nz
  INTEGER :: sNxNz, beg_iz, end_iz, rws, iz1, iz2, msz
  REAL, DIMENSION(:,:), ALLOCATABLE  :: save_P, messages, aux

  REAL, DIMENSION(:,:), ALLOCATABLE :: Vx,Vz,Sxx,Sxz,Szz

  INTEGER :: Nx,Nz
  REAL    :: dx,dz,dt,fpeak
  INTEGER :: itmax
  INTEGER :: lx, lz
  REAL, DIMENSION(:,:), ALLOCATABLE :: C11,C13,C33,C44
  REAL, DIMENSION(:,:), ALLOCATABLE :: rhox,rhoz

  REAL, DIMENSION(:),   ALLOCATABLE :: trs
  REAL    :: tdelay

  REAL, DIMENSION(:), ALLOCATABLE :: taper  
  REAL :: F
  INTEGER :: nb,fsf

  REAL    :: t
  INTEGER :: it,i,j
  !INTEGER :: aux

  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_Vx,csg_Vz,csg_P,csg_S
  INTEGER :: nphones,npmin_x,npmin_z,dnp_x,dnp_z
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER :: isnap,nsnaps,snapmin,dsnap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL MPI_Init(ierr)

  comm = MPI_COMM_WORLD
  CALL MPI_COMM_RANK(comm, rank, ierr)
  CALL MPI_COMM_SIZE(comm, size,   ierr)

  PRINT*, 'Alo, processo ', rank,' de ', size

  proot=0

  IF (rank==proot) THEN

     PRINT*,'anisot - Modelling seismic waves in anisotropic media'
     PRINT*,''

     ! Input
     CALL inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,&
          nphones,npmin_x,npmin_z,dnp_x,dnp_z,nsnaps,snapmin,dsnap,&
          nb,F,fsf)

     OPEN(26,FILE="parallel.dat",STATUS='UNKNOWN')
     READ(26,'(t10,i10)') np

  END IF

  CALL MPI_BCAST(np,      1,mpi_integer,  proot,comm,ierr)
  IF(np.ne.size) THEN
     stop 'Corrigir arquivo parallel.dat para np=size '
  END IF

  ALLOCATE(div_nz(0:np))

  IF (rank==proot) THEN

     DO i=0,np
        READ(26,'(t10,i10)') div_nz(i)
     END DO

  END IF

  CALL MPI_BCAST(Nx,      1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(Nz,      1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(itmax,   1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(dx,      1,mpi_real,     proot,comm,ierr)
  CALL MPI_BCAST(dt,      1,mpi_real,     proot,comm,ierr)
  CALL MPI_BCAST(fpeak,   1,mpi_real,     proot,comm,ierr)
  CALL MPI_BCAST(lx,      1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(lz,      1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(nsnaps,  1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(snapmin, 1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(dsnap,   1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(nphones, 1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(npmin_x, 1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(npmin_z, 1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(dnp_x,   1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(dnp_z,   1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(nb,      1,mpi_integer,  proot,comm,ierr)
  CALL MPI_BCAST(F,       1,mpi_real,     proot,comm,ierr)
  CALL MPI_BCAST(div_nz, np+1,mpi_integer,proot,comm,ierr)
  dz=dx

  beg_iz=div_nz(rank) - 3
  end_iz=div_nz(rank+1)


  ALLOCATE(C11(Nx,beg_iz:end_iz),C13(Nx,beg_iz:end_iz),&
       C33(Nx,beg_iz:end_iz),C44(Nx,beg_iz:end_iz))
  ALLOCATE(rhox(Nx,beg_iz:end_iz),rhoz(Nx,beg_iz:end_iz))

  CALL inputmodelmpi(Nx,Nz,beg_iz,end_iz,C11,C13,C33,C44,rhox,rank,comm,np,div_nz)

  rhoz=rhox
  ! Interpolando valores para malha escalonada

  do j=beg_iz,end_iz-1
     do i=1,Nx-1
        rhoz(i,j) = 0.25*(rhox(i,j) &
             + rhox(i+1,j) + rhox(i+1,j+1) + rhox(i,j+1))
     end do
  end do
  do j=beg_iz,end_iz
     do i=1,Nx-1
        C11(i,j) = 0.5*(C11(i,j) + C11(i+1,j))
        C13(i,j) = 0.5*(C13(i,j) + C13(i+1,j))
        C33(i,j) = 0.5*(C33(i,j) + C33(i+1,j))
     end do
  end do
  do j=beg_iz,end_iz-1
     do i=1,Nx
        C44(i,j) = 0.5*(C44(i,j) + C44(i,j+1))
     end do
  end do

  C11 = C11*1.e10
  C13 = C13*1.e10
  C33 = C33*1.e10
  C44 = C44*1.e10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE(trs(itmax))
  IF ((lz>=beg_iz).and.(lz<end_iz)) THEN
     CALL source_init(trs,dt,itmax,fpeak,tdelay)
     !  PRINT*,'time delay=',tdelay
  END IF

  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,F)

  IF (rank==proot) THEN

     ALLOCATE(csg_Vx(itmax,nphones),csg_Vz(itmax,nphones))
     ALLOCATE(csg_P(itmax,nphones),csg_S(itmax,nphones))

     n1=31; n2=32; n3=33; n4=34; n5=35; n6=36; n7=37; n8=38
     OPEN(n1,FILE="csg_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
     OPEN(n2,FILE="csg_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)

     OPEN(n3,FILE="snapshots_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
     OPEN(n4,FILE="snapshots_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

     OPEN(n5,FILE="csg_P.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
     OPEN(n6,FILE="snapshots_P.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

     OPEN(n7,FILE="csg_S.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
     OPEN(n8,FILE="snapshots_S.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
          ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

     ALLOCATE(save_P(Nx,Nz))

  END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(Vx(Nx,beg_iz:end_iz),Vz(Nx,beg_iz:end_iz),&
       Sxx(Nx,beg_iz:end_iz),Sxz(Nx,beg_iz:end_iz),&
       Szz(Nx,beg_iz:end_iz))

  Sxx=0.;  Sxz=0.;  Szz=0.;  Vx=0.;  Vz=0.
  rhox=1./rhox;  rhoz=1./rhoz

  isnap=0
  ALLOCATE(messages(Nx,10))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Beggining time loop
  DO it=1,itmax
     t=it*dt
     IF ((lz>=beg_iz).and.(lz<end_iz)) THEN

        ! Insert source function
        Sxx(lx,lz) = Sxx(lx,lz) + trs(it)
        Szz(lx,lz) = Szz(lx,lz) + trs(it)
      !  PRINT*,'it',it,'/',itmax,' Vx(lx+2,lz+2)=',Vx(lx+2,lz+2)

        OPEN(77,FILE="status.txt",STATUS='UNKNOWN',ACTION='WRITE')
        WRITE(77,*) &
             'anisot - Modeling seismic waves in anisotropic media'
        WRITE(77,*) &
             'Iteracao',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)
        CLOSE(77)

     END IF

     ! Do one time step ( MODIFICAR SUBROTINA PARA FAIXAS )
     CALL prop_anisot_2D(Vx,Vz,Sxx,Sxz,Szz,&
          C11,C13,C33,C44,rhox,rhoz,dx,dt,1,Nx,beg_iz,end_iz)

     IF (rank==0) THEN
        iz1=beg_iz
        iz2=end_iz
     ELSE
        iz1=div_nz(rank)
        iz2=div_nz(rank+1)
     END IF

     !!!!!!!!!!/* boundary taper (Cerjan) for Vx
     CALL bt_apply_left(Vx(1:nb,iz1:iz2),Nx,iz1,iz2,nb,taper)
     CALL bt_apply_right(Vx(Nx-nb+1:Nx,iz1:iz2),Nx,iz1,iz2,nb,taper)
     IF (rank==(np-1)) THEN
        CALL bt_apply_bottom(Vx(:,Nz-nb+1:Nz),Nx,Nz,nb,taper)
     END IF
     IF (rank==0) THEN
        IF (fsf.ne.1) THEN  !(free surface: fsf=1)
           CALL bt_apply_top(Vx(:,1:nb),Nx,nb,taper)
        END IF
     END IF
     !!!!!!!!!!/* boundary taper (Cerjan) for Vz
     CALL bt_apply_left(Vz(1:nb,iz1:iz2),Nx,iz1,iz2,nb,taper)
     CALL bt_apply_right(Vz(Nx-nb+1:Nx,iz1:iz2),Nx,iz1,iz2,nb,taper)
     IF (rank==(np-1)) THEN
        CALL bt_apply_bottom(Vz(:,Nz-nb+1:Nz),Nx,Nz,nb,taper)
     END IF
     IF (rank==0) THEN
        IF (fsf.ne.1) THEN  !(free surface: fsf=1)
           CALL bt_apply_top(Vz(:,1:nb),Nx,nb,taper)
        END IF
     END IF
     !!!!!!!!!!/* boundary taper (Cerjan) for Sxx
     CALL bt_apply_left(Sxx(1:nb,iz1:iz2),Nx,iz1,iz2,nb,taper)
     CALL bt_apply_right(Sxx(Nx-nb+1:Nx,iz1:iz2),Nx,iz1,iz2,nb,taper)
     IF (rank==(np-1)) THEN
        CALL bt_apply_bottom(Sxx(:,Nz-nb+1:Nz),Nx,Nz,nb,taper)
     END IF
     IF (rank==0) THEN
        IF (fsf.ne.1) THEN  !(free surface: fsf=1)
           CALL bt_apply_top(Sxx(:,1:nb),Nx,nb,taper)
        END IF
     END IF
     !!!!!!!!!!/* boundary taper (Cerjan) for Szz
     CALL bt_apply_left(Szz(1:nb,iz1:iz2),Nx,iz1,iz2,nb,taper)
     CALL bt_apply_right(Szz(Nx-nb+1:Nx,iz1:iz2),Nx,iz1,iz2,nb,taper)
     IF (rank==(np-1)) THEN
        CALL bt_apply_bottom(Szz(:,Nz-nb+1:Nz),Nx,Nz,nb,taper)
     END IF
     IF (rank==0) THEN
        IF (fsf.ne.1) THEN  !(free surface: fsf=1)
           CALL bt_apply_top(Szz(:,1:nb),Nx,nb,taper)
        END IF
     END IF
     !!!!!!!!!!/* boundary taper (Cerjan) for Sxz
     CALL bt_apply_left(Sxz(1:nb,iz1:iz2),Nx,iz1,iz2,nb,taper)
     CALL bt_apply_right(Sxz(Nx-nb+1:Nx,iz1:iz2),Nx,iz1,iz2,nb,taper)
     IF (rank==(np-1)) THEN
        CALL bt_apply_bottom(Sxz(:,Nz-nb+1:Nz),Nx,Nz,nb,taper)
     END IF
     IF (rank==0) THEN
        IF (fsf.ne.1) THEN  !(free surface: fsf=1)
           CALL bt_apply_top(Sxz(:,1:nb),Nx,nb,taper)
        END IF
     END IF

     IF (rank.ne.(np-1)) THEN
        source=rank
        dest=rank+1
        tag=1000*(rank+1)
        messages(1:Nx,1) = Vx(1:Nx,end_iz-3)
        messages(1:Nx,2) = Vx(1:Nx,end_iz-2)
        messages(1:Nx,3) = Vz(1:Nx,end_iz-3)
        messages(1:Nx,4) = Vz(1:Nx,end_iz-2)
        messages(1:Nx,5) =Sxx(1:Nx,end_iz-3)
        messages(1:Nx,6) =Sxx(1:Nx,end_iz-2)
        messages(1:Nx,7) =Sxz(1:Nx,end_iz-3)
        messages(1:Nx,8) =Sxz(1:Nx,end_iz-2)
        messages(1:Nx,9) =Szz(1:Nx,end_iz-3)
        messages(1:Nx,10)=Szz(1:Nx,end_iz-2)
        CALL MPI_SEND(messages,10*Nx,MPI_REAL,&
             dest,tag,comm,ierr)
        source=rank+1
        tag=1000*(rank+1) + 10
        CALL MPI_RECV(messages,10*Nx,MPI_REAL,&
             source,tag,comm,status,ierr)
        Vx(1:Nx,end_iz-1) = messages(1:Nx,1)
        Vx(1:Nx,end_iz)   = messages(1:Nx,2)
        Vz(1:Nx,end_iz-1) = messages(1:Nx,3)
        Vz(1:Nx,end_iz)   = messages(1:Nx,4)
        Sxx(1:Nx,end_iz-1)= messages(1:Nx,5)
        Sxx(1:Nx,end_iz)  = messages(1:Nx,6)
        Sxz(1:Nx,end_iz-1)= messages(1:Nx,7)
        Sxz(1:Nx,end_iz)  = messages(1:Nx,8)
        Szz(1:Nx,end_iz-1)= messages(1:Nx,9)
        Szz(1:Nx,end_iz)  = messages(1:Nx,10)
     END IF

     IF (rank.ne.0) THEN
        source=rank-1
        tag=1000*rank
        CALL MPI_RECV(messages,10*Nx,MPI_REAL,&
             source,tag,comm,status,ierr)
        Vx(1:Nx,beg_iz)   = messages(1:Nx,1)
        Vx(1:Nx,beg_iz+1) = messages(1:Nx,2)
        Vz(1:Nx,beg_iz)   = messages(1:Nx,3)
        Vz(1:Nx,beg_iz+1) = messages(1:Nx,4)
        Sxx(1:Nx,beg_iz)  = messages(1:Nx,5)
        Sxx(1:Nx,beg_iz+1)= messages(1:Nx,6)
        Sxz(1:Nx,beg_iz)  = messages(1:Nx,7)
        Sxz(1:Nx,beg_iz+1)= messages(1:Nx,8)
        Szz(1:Nx,beg_iz)  = messages(1:Nx,9)
        Szz(1:Nx,beg_iz+1)= messages(1:Nx,10)
        source=rank
        dest=rank-1
        tag=1000*rank + 10
        messages(1:Nx,1) = Vx(1:Nx,beg_iz+2)
        messages(1:Nx,2) = Vx(1:Nx,beg_iz+3)
        messages(1:Nx,3) = Vz(1:Nx,beg_iz+2)
        messages(1:Nx,4) = Vz(1:Nx,beg_iz+3)
        messages(1:Nx,5) =Sxx(1:Nx,beg_iz+2)
        messages(1:Nx,6) =Sxx(1:Nx,beg_iz+3)
        messages(1:Nx,7) =Sxz(1:Nx,beg_iz+2)
        messages(1:Nx,8) =Sxz(1:Nx,beg_iz+3)
        messages(1:Nx,9) =Szz(1:Nx,beg_iz+2)
        messages(1:Nx,10)=Szz(1:Nx,beg_iz+3)
        CALL MPI_SEND(messages,10*Nx,MPI_REAL,&
             dest,tag,comm,ierr)
     END IF

     !     PRINT*,'it',it,'/',itmax,'rank',rank
     !/* Output
     IF ((npmin_z>=beg_iz).and.(npmin_z<end_iz)) THEN 
          CALL save_shotgather(csg_Vx,&
               Vx(:,npmin_z),Nx,nphones,npmin_x,npmin_z,dnp_x,&
               n1,it,itmax)
          CALL save_shotgather(csg_Vz,&
               Vx(:,npmin_z),Nx,nphones,npmin_x,npmin_z,dnp_x,&
               n2,it,itmax)
     END IF
     CALL save_snapshots(Vx,&
        isnap,nsnaps,snapmin,dsnap,n3,Nx,Nz,it,itmax,&
        rank,np,beg_iz,end_iz,div_nz,comm,status,save_P)
     CALL save_snapshots(Vz,&
        isnap,nsnaps,snapmin,dsnap,n4,Nx,Nz,it,itmax,&
        rank,np,beg_iz,end_iz,div_nz,comm,status,save_P)

     ! Output
     !     CALL save_shotgather_n_snapshots(csg_Vx,csg_Vz,&
     !          csg_P,csg_S,Vx,Vz,&
     !          Sxx,Sxz,Szz,nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
     !          isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,n7,n8,&
     !          Nx,Nz,it,itmax)

     ! Boundary taper
     !     CALL bt_apply_multiple(Vx,Vz,Sxx,Sxz,Szz,Nx,Nz,nb,taper,fsf)

  END DO

  PRINT*,'Successful run.'

  CALL MPI_Finalize(ierr)

END PROGRAM anisot
