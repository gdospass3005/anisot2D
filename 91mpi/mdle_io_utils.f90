MODULE mdle_io_utils

CONTAINS

  SUBROUTINE save_shotgather(csg_ux,&
       next_ux,Nx,nphones,npmin_x,np_z,dnp_x,&
       n1,it,itmax)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: csg_ux
    INTEGER, INTENT(IN)    :: Nx,nphones
    INTEGER, INTENT(IN)    :: npmin_x,np_z,dnp_x
    REAL, DIMENSION(1:Nx,np_z:np_z), INTENT(IN)    :: next_ux
    INTEGER, INTENT(IN)    :: n1
    INTEGER, INTENT(IN)    :: it,itmax
    INTEGER :: i,j,ix,iz

    iz = np_z
    do i=1,nphones
       ix = (i-1)*dnp_x + npmin_x
       csg_ux(it,i) = next_ux(ix,iz)
    end do

    if (it == itmax) then  
       PRINT*,'Saving CSG'
       do i=1,nphones
          WRITE(n1, REC=i) csg_ux(:,i)
       end do
    end if

  END SUBROUTINE save_shotgather


  SUBROUTINE save_snapshots(next_P,&
       isnap,nsnaps,snapmin,dsnap,n3,Nx,Nz,it,itmax,&
       rank,np,beg_iz,end_iz,div_nz,comm,status,save_P)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN)    :: comm,status(MPI_STATUS_SIZE)
    INTEGER, INTENT(INOUT) :: isnap
    INTEGER, INTENT(IN)    :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(IN)    :: n3
    INTEGER, INTENT(IN)    :: Nx,Nz,it,itmax,rank,np,beg_iz,end_iz
    REAL,    DIMENSION(Nx,beg_iz:end_iz), INTENT(IN) :: next_P
    INTEGER, DIMENSION(0:np), INTENT(IN) :: div_nz
    INTEGER :: source,dest,tag, ierr
    INTEGER :: i,j,ix,iz,iz1,iz2,msz
    REAL, DIMENSION(:,:), INTENT(INOUT) :: save_P

    if (it == (isnap * dsnap) + snapmin) then
       isnap=isnap+1
       if (isnap <= nsnaps) then

          IF (rank.ne.0) THEN
             dest=0
             tag=1111*rank
             iz1=div_nz(rank)+1
             iz2=div_nz(rank+1)
             msz = iz2-iz1+1
             CALL MPI_SEND(next_P(:,iz1:iz2),Nx*msz,MPI_REAL,&
                  dest,tag,comm,ierr)
          END IF

          IF (rank==0) THEN
             PRINT*, 'Saving snapshot ',isnap,'/',nsnaps
             iz1=beg_iz
             iz2=end_iz
             save_P(:,iz1:iz2)=next_P(:,iz1:iz2)

             do i=1,np-1
                source=i
                tag=1111*i
                iz1=div_nz(i)+1
                iz2=div_nz(i+1)
                msz = iz2-iz1+1
                CALL MPI_RECV(save_P(:,iz1:iz2),Nx*msz,MPI_REAL,&
                     source,tag,comm,status,ierr)
             end do
             do j=1,Nx
                WRITE(n3, REC=((isnap -1)*Nx + j)) save_P(j,:)
             end do
          END IF

       end if
    end if

  END SUBROUTINE save_snapshots


  SUBROUTINE inputmodelmpi(Nx,Nz,beg_iz,end_iz,C11,C13,C33,C44,rhox,&
       rank,comm,np,div_nz)
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER, INTENT(IN) :: Nx,Nz,beg_iz,end_iz
    INTEGER, INTENT(IN) :: rank,comm,np
    INTEGER, INTENT(IN), DIMENSION(0:np) :: div_nz
    REAL,    INTENT(OUT), DIMENSION(1:Nx,beg_iz:end_iz) :: C11,C13,C33,C44,rhox
    INTEGER :: proot=0,source,dest,tag,iz1,iz2,msz,ierr,status(MPI_STATUS_SIZE)
    INTEGER :: i
    REAL,    DIMENSION(:,:), ALLOCATABLE :: aux

    IF (rank==proot) THEN
       ALLOCATE(aux(Nx,Nz))
       OPEN(34,FILE="C11.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
            ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
       DO i=1,Nx
          READ(34, REC=i) aux(i,:)
       END DO
       CLOSE(34)

       C11 = aux(:,beg_iz:end_iz)

       DO i=1,np-1
          source=0
          tag=100000*i
          iz1=div_nz(i)-3
          iz2=div_nz(i+1)
          msz = iz2-iz1+1
          dest=i
          CALL MPI_SEND(aux(:,iz1:iz2),Nx*msz,MPI_REAL,&
               dest,tag,comm,ierr)
       END DO
    !   DEALLOCATE(aux)

    END IF

    IF (rank.ne.proot) THEN
       source=0
       tag=100000*rank
       iz1=beg_iz
       iz2=end_iz
       msz = iz2-iz1+1
       CALL MPI_RECV(C11(:,iz1:iz2),Nx*msz,MPI_REAL,&
            source,tag,comm,status,ierr)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (rank==proot) THEN
     !  ALLOCATE(aux(Nx,Nz))

       OPEN(34,FILE="C33.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
            ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
       DO i=1,Nx
          READ(34, REC=i) aux(i,:)
       END DO
       CLOSE(34)

       C33 = aux(:,beg_iz:end_iz)

       DO i=1,np-1
          source=0
          tag=100033*i
          iz1=div_nz(i)-3
          iz2=div_nz(i+1)
          msz = iz2-iz1+1
          dest=i
          CALL MPI_SEND(aux(:,iz1:iz2),Nx*msz,MPI_REAL,&
               dest,tag,comm,ierr)
       END DO
      ! DEALLOCATE(aux)

    END IF

    IF (rank.ne.proot) THEN
       source=0
       tag=100033*rank
       iz1=beg_iz
       iz2=end_iz
       msz = iz2-iz1+1
       CALL MPI_RECV(C33(:,iz1:iz2),Nx*msz,MPI_REAL,&
            source,tag,comm,status,ierr)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (rank==proot) THEN
      ! ALLOCATE(aux(Nx,Nz))

       OPEN(34,FILE="C13.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
            ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
       DO i=1,Nx
          READ(34, REC=i) aux(i,:)
       END DO
       CLOSE(34)

       C13 = aux(:,beg_iz:end_iz)

       DO i=1,np-1
          source=0
          tag=100013*i
          iz1=div_nz(i)-3
          iz2=div_nz(i+1)
          msz = iz2-iz1+1
          dest=i
          CALL MPI_SEND(aux(:,iz1:iz2),Nx*msz,MPI_REAL,&
               dest,tag,comm,ierr)
       END DO
      ! DEALLOCATE(aux)

    END IF

    IF (rank.ne.proot) THEN
       source=0
       tag=100013*rank
       iz1=beg_iz
       iz2=end_iz
       msz = iz2-iz1+1
       CALL MPI_RECV(C13(:,iz1:iz2),Nx*msz,MPI_REAL,&
            source,tag,comm,status,ierr)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (rank==proot) THEN
      ! ALLOCATE(aux(Nx,Nz))

       OPEN(34,FILE="C44.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
            ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
       DO i=1,Nx
          READ(34, REC=i) aux(i,:)
       END DO
       CLOSE(34)

       C44 = aux(:,beg_iz:end_iz)

       DO i=1,np-1
          source=0
          tag=100044*i
          iz1=div_nz(i)-3
          iz2=div_nz(i+1)
          msz = iz2-iz1+1
          dest=i
          CALL MPI_SEND(aux(:,iz1:iz2),Nx*msz,MPI_REAL,&
               dest,tag,comm,ierr)
       END DO
      ! DEALLOCATE(aux)

    END IF

    IF (rank.ne.proot) THEN
       source=0
       tag=100044*rank
       iz1=beg_iz
       iz2=end_iz
       msz = iz2-iz1+1
       CALL MPI_RECV(C44(:,iz1:iz2),Nx*msz,MPI_REAL,&
            source,tag,comm,status,ierr)
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (rank==proot) THEN
      ! ALLOCATE(aux(Nx,Nz))

       OPEN(34,FILE="rhox.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
            ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
       DO i=1,Nx
          READ(34, REC=i) aux(i,:)
       END DO
       CLOSE(34)

       rhox = aux(:,beg_iz:end_iz)

       DO i=1,np-1
          source=0
          tag=100222*i
          iz1=div_nz(i)-3
          iz2=div_nz(i+1)
          msz = iz2-iz1+1
          dest=i
          CALL MPI_SEND(aux(:,iz1:iz2),Nx*msz,MPI_REAL,&
               dest,tag,comm,ierr)
       END DO
       DEALLOCATE(aux)

    END IF

    IF (rank.ne.proot) THEN
       source=0
       tag=100222*rank
       iz1=beg_iz
       iz2=end_iz
       msz = iz2-iz1+1
       CALL MPI_RECV(rhox(:,iz1:iz2),Nx*msz,MPI_REAL,&
            source,tag,comm,status,ierr)
    END IF

  END SUBROUTINE inputmodelmpi


  SUBROUTINE inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,nsnaps,snapmin,dsnap,&
       nb,F,fsf)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: Nx,Nz,itmax,lx,lz
    REAL,    INTENT(OUT) :: dx,dt,fpeak
    INTEGER, INTENT(OUT) :: nphones,npmin_x,npmin_z,dnp_x,dnp_z
    INTEGER, INTENT(OUT) :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(OUT) :: nb,fsf
    REAL,    INTENT(OUT) :: F
 
    OPEN(30,FILE="data.dat",STATUS='UNKNOWN',ACTION='READ')
    READ(30,'(t10,i10)') Nx
    READ(30,'(t10,i10)') Nz
    READ(30,'(t10,f10.4)') dx
    READ(30,'(t10,f10.8)') dt
    READ(30,'(t10,f10.4)') fpeak
    READ(30,'(t10,i10)') itmax
    READ(30,'(t10,i10)') lx
    READ(30,'(t10,i10)') lz
    READ(30,'(t10,i10)') nb
    READ(30,'(t10,f10.4)') F
    READ(30,'(t10,i10)') fsf
    READ(30,'(t10,i10)') nphones
    READ(30,'(t10,i10)') npmin_x
    READ(30,'(t10,i10)') npmin_z
    READ(30,'(t10,i10)') dnp_x
    READ(30,'(t10,i10)') dnp_z
    READ(30,'(t10,i10)') nsnaps
    READ(30,'(t10,i10)') snapmin
    READ(30,'(t10,i10)') dsnap
  END SUBROUTINE inputdata


  SUBROUTINE save_shotgather_n_snapshots(csg_ux,csg_uz,csg_P,csg_S,&
       next_ux,next_uz,Sxx,Sxz,Szz,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
       isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,n7,n8,&
       Nx,Nz,it,itmax)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: csg_ux,csg_uz,csg_P,csg_S
    REAL, DIMENSION(:,:), INTENT(IN)    :: next_ux,next_uz,&
         Sxx,Sxz,Szz
    REAL, DIMENSION(Nx) :: aux1
    REAL, DIMENSION(Nz) :: aux2
    INTEGER, INTENT(IN)    :: nphones
    INTEGER, INTENT(IN)    :: npmin_x,npmin_z,dnp_x,dnp_z
    INTEGER, INTENT(INOUT) :: isnap
    INTEGER, INTENT(IN)    :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(IN)    :: n1,n2,n3,n4,n5,n6,n7,n8
    INTEGER, INTENT(IN)    :: Nx,Nz,it,itmax
    INTEGER :: i,ix,iz

    do i=1,nphones
       ix = (i-1)*dnp_x + npmin_x
       iz = (i-1)*dnp_z + npmin_z
       csg_ux(it,i) = next_ux(ix,iz)
       csg_uz(it,i) = next_uz(ix,iz)
       csg_P(it,i) = Sxx(ix,iz) + Szz(ix,iz)
       csg_S(it,i) = Sxz(ix,iz)
    end do

    if (it == itmax) then  
       PRINT*,'Saving CSG'
       do i=1,nphones
          WRITE(n1, REC=i) csg_ux(:,i)
          WRITE(n2, REC=i) csg_uz(:,i)
          WRITE(n5, REC=i) csg_P(:,i)
          WRITE(n7, REC=i) csg_S(:,i)
       end do
    end if

    if (it == (isnap * dsnap) + snapmin) then
       isnap=isnap+1
       if (isnap <= nsnaps) then
          PRINT*, 'Saving snapshot ',isnap,'/',nsnaps
          do i=1,Nx
             WRITE(n3, REC=((isnap -1)*Nx + i)) next_ux(i,:)
             WRITE(n4, REC=((isnap -1)*Nx + i)) next_uz(i,:)
             WRITE(n6, REC=((isnap -1)*Nx + i)) Sxx(i,:)+Szz(i,:)
             WRITE(n8, REC=((isnap -1)*Nx + i)) Sxz(i,:)
          end do
       end if
    end if

  END SUBROUTINE save_shotgather_n_snapshots

END MODULE mdle_io_utils


