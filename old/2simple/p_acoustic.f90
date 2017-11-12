PROGRAM acoustic

  USE mdle_source
  USE mdle_prop
  USE mdle_bpar
  USE mdle_taper
  USE mdle_io_utils

  !Last edit:              SSA, 17/04/2002. gfp.

  IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! VARIABLES DEFINITION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: Nx,Nz           ! grid shape
  REAL    :: dx, dz          ! x, z sample intervals
  INTEGER :: itmax           ! maximum time sample to compute
  REAL, DIMENSION(:,:), ALLOCATABLE  :: Vel
  REAL, DIMENSION(:,:), TARGET, ALLOCATABLE  :: pm,p,pp
  REAL, DIMENSION(:,:), POINTER :: prev_p,next_p,actl_p,swap
  REAL, DIMENSION(:), ALLOCATABLE :: taper
  REAL, DIMENSION(:), ALLOCATABLE :: trs       ! source memory function
  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_P
  REAL :: dt              ! time sample interval
  REAL :: fpeak           ! peak frequency of ricker wavelet
  REAL :: t,tdelay        ! time, time shift of the source function
  INTEGER :: lx, lz       ! location of the source in grid units
  INTEGER :: nb           ! n. of boundary points
  REAL    :: F            ! absorption factor
  INTEGER :: dsnap             ! time samples between snapshots
  INTEGER :: nsnaps,snapmin    ! n. snapshots, 1st. snapshot
  INTEGER :: nphones,npmin_x,npmin_z,dnp_x,dnp_z  ! phone locations
  CHARACTER(LEN=20) :: infile='data.dat'
  CHARACTER(LEN=20) :: Velfile='model.ad'
  CHARACTER(LEN=20) :: reportfile='report.txt'  ! file for report
  INTEGER :: i,it,isnap   ! counters
  INTEGER :: n1,n2
  REAL :: maxVel
  INTEGER :: fsf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PRINT*, 'ACOUSTIC: 2-D Finite difference modeling the acoustic wave &
       &equation'

  OPEN(25,FILE=infile,STATUS='UNKNOWN')

  READ(25,'(t10,i10)') Nx
  READ(25,'(t10,i10)') Nz
  READ(25,'(t10,i10)') itmax
  READ(25,'(t10,f10.4)') dx
  READ(25,'(t10,f10.4)') dt
  READ(25,'(t10,f10.4)') fpeak
  READ(25,'(t10,i10)') lx
  READ(25,'(t10,i10)') lz
  READ(25,'(t10,i10)') nsnaps
  READ(25,'(t10,i10)') snapmin
  READ(25,'(t10,i10)') dsnap
  READ(25,'(t10,i10)') nphones
  READ(25,'(t10,i10)') npmin_x
  READ(25,'(t10,i10)') npmin_z
  READ(25,'(t10,i10)') dnp_x
  READ(25,'(t10,i10)') dnp_z
  READ(25,'(t10,i10)') nb
  READ(25,'(t10,f10.4)') F
  READ(25,'(t10,i10)') fsf

  dz=dx

  ALLOCATE(Vel(Nx,Nz))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Uploading the velocity panel
  OPEN(34,FILE=Velfile,STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
  DO i=1,Nx 
     READ(34, REC=i) Vel(i,:)
  END DO
  CLOSE(34)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !/* compute density*velocity^2 and 1/density */
  maxVel = MAXVAL(Vel)
  Vel = Vel*Vel

  ! Initialize boundary taper
  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,F)

  ! Initialize source vector
  ALLOCATE(trs(itmax))
  call source_init(trs,dt,itmax,fpeak,tdelay)

  n1=31; n2=32
  OPEN(n1,FILE="csg_p.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n2,FILE="snapshots_p.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  ALLOCATE(csg_P(itmax,nphones))

  ALLOCATE(pm(Nx,Nz))
  ALLOCATE(p(Nx,Nz))
  ALLOCATE(pp(Nx,Nz))
  prev_p => pm
  actl_p  => p
  next_p => pp

  prev_p = 0.
  actl_p = 0.
  next_p = 0.

  isnap = 0

  !/* Loop over time steps
  DO it=1,itmax
     !PRINT*, 'it=',it,'/',itmax,'  actl_p(lx,lz+10)=',actl_p(lx,lz+10)

     t = (it-1.)*dt
     OPEN(77,FILE="status.txt",STATUS='UNKNOWN',ACTION='WRITE')
     WRITE(77,*) &
          'acoustic - Modelling seismic waves in acoustic media'
     WRITE(77,*) &
          'Iteration',it,'/',itmax,' actl_p(lx+5,lz+5)=',actl_p(lx+5,lz+5)
     CLOSE(77)

     prev_p(lx,lz) = prev_p(lx,lz) + trs(it)

     !/* do one time step
     CALL prop_wave(next_p,actl_p,prev_p,Vel,maxVel,&
                     1,Nx,1,Nz,dx,dz,dt,t,lx,lz)

     !/* paraxial border (Reynolds)
     CALL bpar_left(next_p,actl_p,prev_p,Vel,&
                     Nx,Nz,dx,dz,dt)
     CALL bpar_right(next_p,actl_p,prev_p,Vel,&
                     Nx,Nz,dx,dz,dt)
     IF (fsf.ne.1) THEN
        CALL bpar_top(next_p,actl_p,prev_p,Vel,&
                     Nx,Nz,dx,dz,dt)
     END IF
     CALL bpar_bottom(next_p,actl_p,prev_p,Vel,&
                     Nx,Nz,dx,dz,dt)

     !/* boundary taper (Cerjan)
     !CALL bt_apply_left(next_p,Nx,Nz,nb,taper)
     !CALL bt_apply_right(next_p,Nx,Nz,nb,taper)
     ! CALL bt_apply_top(next_p,Nx,Nz,nb,taper)
     !CALL bt_apply_bottom(next_p,Nx,Nz,nb,taper)

     !/* save snapshot (actl_p) and update/save csg
     CALL save_shotgather_n_snapshots(csg_P,&
             actl_p,nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
             isnap,nsnaps,snapmin,dsnap,n1,n2,Nx,Nz,it,itmax)

     !/* roll time slice pointers
     swap   => prev_p
     prev_p => actl_p
     actl_p => next_p
     next_p => swap

     !/* end of loop over time steps
  END DO


  !/* saves useful information in the report file
  OPEN(28,FILE=reportfile,STATUS='UNKNOWN')
  WRITE(28,*) 'ACOUSTIC: 2-D Finite difference modelling the acoustic wave &
       &equation'
  WRITE(28,*) ' '
  WRITE(28,*) '---REPORT--- '
  WRITE(28,*) ' '
  WRITE(28,*) 'Nx    =    ',Nx
  WRITE(28,*) 'Nz    =    ',Nz
  WRITE(28,*) 'dx    =    ',dx
  WRITE(28,*) 'dz    =    ',dz
  WRITE(28,*) ' '
  WRITE(28,*) ' '
  WRITE(28,*) 'itmax =    ',itmax
  WRITE(28,*) 'dt    =    ',dt
  WRITE(28,*) ' '
  WRITE(28,*) ' '
  WRITE(28,*) 'fpeak =    ',fpeak
  WRITE(28,*) 'tdelay=    ',tdelay


END PROGRAM acoustic







