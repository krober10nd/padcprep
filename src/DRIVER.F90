PROGRAM Driver
!
USE PRE, ONLY: X,Y,NNELG,READGRAPH 
USE PARPREP
USE MESSENGER, ONLY: MSG_INIT, MSG_FINI 
IMPLICIT NONE
!---------------------------------------------------------------------
!      P R O G R A M   D R I V E R
!---------------------------------------------------------------------
!  EXECUTES PROGRAM TO SOLVE 2-D WAVE EQUATION WITH DAMPING IN PARALLEL USING
!  GHOST NODES. THE EQUATION IS FORMULATED  USING THE WEAK GALERKIN METHOD WITH LINEAR TRIANGULAR ELEMENTS.   
!  KJR, UND, CHL, 2016 
!---------------------------------------------------------------------

CALL MSG_INIT !start up MPI
CALL READGRAPH !prompt the user for the graph read into all PE's
!DECOMPOSE USING PARMETIS 
CALL DECOMPOSE
!
!THEN TIMESTEP 
CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
