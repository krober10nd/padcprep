PROGRAM Driver

!
USE PRE, ONLY: READGRAPH 
USE PARPREP, ONLY: DECOMPOSE_PAR,PREPGRID
USE MESSENGER 
IMPLICIT NONE
!---------------------------------------------------------------------
!      P R O G R A M   D R I V E R
!---------------------------------------------------------------------
!  KJR, UND, CHL, 2016 
!---------------------------------------------------------------------

CALL MSG_INIT !start up MPI
CALL READGRAPH 
CALL DECOMPOSE_PAR !call parmetis
CALL PREPGRID ! puts the decomposed graph in each directory
CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
