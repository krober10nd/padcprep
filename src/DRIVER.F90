PROGRAM Driver
!
USE PRE, ONLY: X,Y,NNELG,READGRAPH 
USE PARPREP, ONLY: BUILD,DECOMPOSE_PAR
USE MESSENGER 
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

CALL READGRAPH 
!CALL BUILD !build adj and xadj
! CALL DECOMPOSE_PAR   
!THEN TIMESTEP 

CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
