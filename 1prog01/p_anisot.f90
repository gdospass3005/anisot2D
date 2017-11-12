PROGRAM anisot

  USE mdle_io_utils                                                  
  USE mdle_source
  USE mdle_prop
  USE mdle_taper

  IMPLICIT NONE
  REAL, DIMENSION(:,:), ALLOCATABLE :: Vx,Vz,Sxx,Sxz,Szz

  INTEGER :: Nx,Nz
  REAL    :: dx,dt,fpeak
  INTEGER :: itmax
  INTEGER :: lx, lz
  REAL    :: r
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
  PRINT*,'anisot - Modeling seismic waves in anisotropic media'
  PRINT*,''

  ! Input
  CALL inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,r,&
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(trs(itmax))
  CALL source_dervgauss_init(trs,dt,itmax,fpeak,tdelay)
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

  END DO

  PRINT*,'Successful run.'

END PROGRAM anisot
