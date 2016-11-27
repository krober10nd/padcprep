MODULE MESSENGER 
!THIS MODULE SUPPLIES THE MESSAGE PASSING INTERFACE 
USE MPI 
IMPLICIT NONE 

INTEGER :: MYPROC  
INTEGER :: IERR  ! error status of mpi subroutine call
INTEGER :: NPROC !number of processors 
CONTAINS 

!---------------------------------------------------------------------
!  STARTS UP MPI COMMUNCIATION   
!---------------------------------------------------------------------
      SUBROUTINE MSG_INIT
      IMPLICIT NONE
 
      INTEGER, ALLOCATABLE    :: RANKS(:)  ! array of mpi ranks for compute processors
      INTEGER                 :: I
!
!
!.......Initialize MPI
!
      CALL MPI_INIT (IERR)
!
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPROC,IERR)   ! Get number of procs
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPROC,IERR)      ! Get MPI rank
      !PRINT*, 'PE', MYPROC, ': online'
     ! ALLOCATE(RANKS(NPROC+1))

    !  DO I=1,NPROC
    !     RANKS(I) = I-1
    !  ENDDO
      RETURN 
!---------------------------------------------------------------------
      END SUBROUTINE MSG_INIT
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!      S U B R O U T I N E   M S G _ F I N I 
!---------------------------------------------------------------------
!  SHUTS DOWN MPI COMMUNCIATION   
!---------------------------------------------------------------------
      SUBROUTINE MSG_FINI
      
      IMPLICIT NONE

        CALL MPI_FINALIZE(IERR)
          
        IF (MYPROC.EQ.0) THEN 
            PRINT *, "MPI terminated with Status = ",IERR
         ENDIF

     RETURN
!---------------------------------------------------------------------      
      END SUBROUTINE MSG_FINI
!---------------------------------------------------------------------

END MODULE MESSENGER 
