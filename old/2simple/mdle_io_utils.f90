MODULE mdle_io_utils

CONTAINS
  SUBROUTINE save_shotgather_n_snapshots(csg_ux,&
       next_ux,nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
       isnap,nsnaps,snapmin,dsnap,n1,n3,Nx,Nz,it,itmax)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: csg_ux
    REAL, DIMENSION(:,:), INTENT(IN)    :: next_ux
    INTEGER, INTENT(IN)    :: nphones
    INTEGER, INTENT(IN)    :: npmin_x,npmin_z,dnp_x,dnp_z
    INTEGER, INTENT(INOUT) :: isnap
    INTEGER, INTENT(IN)    :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(IN)    :: n1,n3
    INTEGER, INTENT(IN)    :: Nx,Nz,it,itmax
    INTEGER :: i,j,ix,iz

    do i=1,nphones
       ix = (i-1)*dnp_x + npmin_x
       iz = (i-1)*dnp_z + npmin_z
       csg_ux(it,i) = next_ux(ix,iz)
    end do

    if (it == itmax) then  
       PRINT*,'Saving CSG'
       do i=1,nphones
          WRITE(n1, REC=i) csg_ux(:,i)
       end do
    end if

    if (it == (isnap * dsnap) + snapmin) then
       isnap=isnap+1
       if (isnap <= nsnaps) then
          PRINT*, 'Saving snapshot ',isnap,'/',nsnaps
          do i=1,Nx
             WRITE(n3, REC=((isnap -1)*Nx + i)) next_ux(i,:)
          end do
       end if
    end if

  END SUBROUTINE save_shotgather_n_snapshots

END MODULE mdle_io_utils
