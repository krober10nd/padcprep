PROGRAM Driver
!
USE PRE, ONLY: X,Y,NNELG,READGRAPH 
USE PARPREP
USE MESSENGER, ONLY: MSG_INIT, MSG_FINI 
IMPLICIT NONE
!---------------------------------------------------------------------
!      P R O G R A M   D R I V E R
!---------------------------------------------------------------------
!  A program to model a traveling wave,  U0(x + ct) where c is
!  a parameter denoting the speed of the wave propogation, in parallel while
!  moving vertices of the graph between domains. The domain decomposition is
!  accomplished by using parmetis 4.0 to decompose the underlying triangular mesh.
!  KJR, UND, CHL, 2016 
!---------------------------------------------------------------------

CALL MSG_INIT !start up MPI
IF(MYPROC.EQ.0) THEN
CALL READGRAPH !prompt the user for the graph read into all PE's
CALL DECOMPOSE_SERIAL !build adj and xadj
ENDIF 
!THEN TIMESTEP 
CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
