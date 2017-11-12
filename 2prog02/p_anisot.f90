PROGRAM anisot

  USE mdle_io_utils                                                  
  USE mdle_source
  USE mdle_prop
  USE mdle_taper
  USE mdle_bpar

  IMPLICIT NONE
  REAL, DIMENSION(:,:), ALLOCATABLE :: Vx,Vz,Sxx,Sxz,Szz
  REAL, DIMENSION(:,:), ALLOCATABLE  :: act_Sxx_left,act_Sxx_right
  REAL, DIMENSION(:,:), ALLOCATABLE  :: act_Sxx_top,act_Sxx_bottom
  REAL, DIMENSION(:,:), ALLOCATABLE  :: old_Sxx_left,old_Sxx_right
  REAL, DIMENSION(:,:), ALLOCATABLE  :: old_Sxx_top,old_Sxx_bottom,Vel

  INTEGER :: Nx,Nz
  REAL    :: dx,dt,fpeak
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

  INTEGER :: beg_iz,end_iz,iz1,iz2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRINT*,'anisot - Modeling seismic waves in anisotropic media'
  PRINT*,''

  ! Input
  CALL inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,nsnaps,snapmin,dsnap,&
       nb,F,fsf)

  ALLOCATE(C11(Nx,Nz),C13(Nx,Nz),C33(Nx,Nz),C44(Nx,Nz))
  ALLOCATE(rhox(Nx,Nz),rhoz(Nx,Nz))

  CALL inputmodel(Nx,Nz,C11,C13,C33,C44,rhox)
  do j=1,Nz-1
     do i=1,Nx-1
        rhoz(i,j) = 0.25*(rhox(i,j) &
             + rhox(i+1,j) + rhox(i+1,j+1) + rhox(i,j+1))
     end do
  end do
  rhoz(Nx,:) = rhox(Nx,:)
  rhoz(:,Nz) = rhox(:,Nz)   
  do j=1,Nz
     do i=1,Nx-1
        C11(i,j) = 0.5*(C11(i,j) + C11(i+1,j))
        C13(i,j) = 0.5*(C13(i,j) + C13(i+1,j))
        C33(i,j) = 0.5*(C33(i,j) + C33(i+1,j))
     end do
  end do
  do j=1,Nz-1
     do i=1,Nx
        C44(i,j) = 0.5*(C44(i,j) + C44(i,j+1))
     end do
  end do

  C11 = C11*1.e10
  C13 = C13*1.e10
  C33 = C33*1.e10
  C44 = C44*1.e10

 ALLOCATE(Vel(Nx,Nz)) 
 Vel = SQRT(C33/rhox)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(trs(itmax))
  CALL source_init(trs,dt,itmax,fpeak,tdelay)
  PRINT*,'time delay=',tdelay

  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,F)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  beg_iz=1; end_iz=Nz
  iz1=beg_iz ; iz2=end_iz

  ALLOCATE(act_Sxx_left(1:4,beg_iz:end_iz),act_Sxx_right(Nx-3:Nx,beg_iz:end_iz))
  ALLOCATE(old_Sxx_left(1:4,beg_iz:end_iz),old_Sxx_right(Nx-3:Nx,beg_iz:end_iz))

  act_Sxx_left = 0.;  act_Sxx_right = 0.
  old_Sxx_left = 0.;  old_Sxx_right = 0.

  !IF (rank==0) THEN
     IF (fsf.ne.1) THEN  !(free surface: fsf=1)
        ALLOCATE(act_Sxx_top(1:Nx,1:4))
        act_Sxx_top = 0.
     END IF
  !END IF
  !IF (rank==(np-1)) THEN
     ALLOCATE(act_Sxx_bottom(1:Nx,Nz-3:Nz))
     act_Sxx_bottom = 0.
  !END IF
  !IF (rank==0) THEN
     IF (fsf.ne.1) THEN  !(free surface: fsf=1)
        ALLOCATE(old_Sxx_top(1:Nx,1:4))
        old_Sxx_top = 0.
     END IF
  !END IF
  !IF (rank==(np-1)) THEN
     ALLOCATE(old_Sxx_bottom(1:Nx,Nz-3:Nz))
     old_Sxx_bottom = 0.
  !END IF

  ALLOCATE(Vx(Nx,Nz),Vz(Nx,Nz),Sxx(Nx,Nz),Sxz(Nx,Nz),Szz(Nx,Nz))
  Sxx=0.; Sxz=0.; Szz=0.; Vx=0.; Vz=0.
  rhox=1./rhox; rhoz=1./rhoz
  isnap=0

  ! Beggining time loop
  DO it=1,itmax
     t=it*dt
     OPEN(77,FILE="status.txt",STATUS='UNKNOWN',ACTION='WRITE')
     WRITE(77,*) &
          'anisot - Modeling seismic waves in anisotropic media'
     WRITE(77,*) &
          'Iteracao',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)
     CLOSE(77)

     ! Insert source function
     Sxx(lx,lz) = Sxx(lx,lz)+trs(it)
     Szz(lx,lz) = Szz(lx,lz)+trs(it)

      !/* paraxial border (Reynolds)
     CALL bpar_left(Sxx(1:2,iz1:iz2),act_Sxx_left(1:4,iz1:iz2),&
          old_Sxx_left(1:4,iz1:iz2),Vel(1:4,iz1:iz2),iz1,iz2,&
          dx,dx,dt)
     CALL bpar_right(Sxx(Nx-1:Nx,iz1:iz2),act_Sxx_right(Nx-3:Nx,iz1:iz2),&
          old_Sxx_right(Nx-3:Nx,iz1:iz2),&
          Vel(Nx-3:Nx,iz1:iz2),Nx,iz1,iz2,dx,dx,dt)

    ! IF (rank==(np-1)) THEN
    !     CALL bpar_bottom(prev_p(1:Nx,Nz-1:Nz),actl_p(1:Nx,Nz-3:Nz),&
    !         old_p_bottom(1:Nx,Nz-3:Nz),&
    !         Vel(1:Nx,Nz-3:Nz),Nx,Nz,dx,dx,dt)
    ! END IF
    ! IF (rank==0) THEN
    !    IF (fsf.ne.1) THEN
    !       CALL bpar_top(prev_p(1:Nx,1:2),actl_p(1:Nx,1:4),&
    !            old_p_top(1:Nx,1:4),Vel(1:Nx,1:4),&
    !            Nx,Nz,dx,dx,dt)
    !    END IF
    ! END IF

     ! Do one time step
     CALL prop_anisot_2D(Vx,Vz,Sxx,Sxz,Szz,&
          C11,C13,C33,C44,rhox,rhoz,dx,dt,Nx,Nz)

     PRINT*,'it',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)

     ! Output
     CALL save_shotgather_n_snapshots(csg_Vx,csg_Vz,&
          csg_P,csg_S,Vx,Vz,&
          Sxx,Sxz,Szz,nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
          isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,n7,n8,&
          Nx,Nz,it,itmax)

     ! Boundary taper
     CALL bt_apply_multiple(Vx,Vz,Sxx,Sxz,Szz,Nx,Nz,nb,taper,fsf)

     old_Sxx_left(1:4,iz1:iz2) = act_Sxx_left(1:4,iz1:iz2)
     act_Sxx_left(1:4,iz1:iz2) = Sxx(1:4,iz1:iz2)

     old_Sxx_right(Nx-3:Nx,iz1:iz2) = act_Sxx_right(Nx-3:Nx,iz1:iz2)
     act_Sxx_right(Nx-3:Nx,iz1:iz2) = Sxx(Nx-3:Nx,iz1:iz2)

    ! IF (rank==0) THEN
    !    IF (fsf.ne.1) THEN  !(free surface: fsf=1)
    !       old_p_top(1:Nx,1:4)=Sxx(1:Nx,1:4)
    !    END IF
    ! END IF
    ! IF (rank==(np-1)) THEN
    !    old_p_bottom(1:Nx,Nz-3:Nz)=Sxx(1:Nx,Nz-3:Nz)
    ! END IF


  END DO

  PRINT*,'Successful run.'

END PROGRAM anisot
