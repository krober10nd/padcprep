PROGRAM Driver

USE PRE, ONLY: READGRAPH,IT,VTXWGTS_LOC,VSIZES_LOC,X_G,Y_G,DP_G,CHUNK
USE PARPREP, ONLY: DECOMPOSE_PAR,EIND
USE MESSENGER 
IMPLICIT NONE

!---------------------------------------------------------------------
!      P R O G R A M   D R I V E R
!---------------------------------------------------------------------
!  KJR, UND, CHL, 2016 
!---------------------------------------------------------------------


CALL MSG_INIT !start up MPI
CALL READGRAPH 
CALL DECOMPOSE_PAR 
CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
