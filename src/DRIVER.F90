PROGRAM Driver

!
USE PRE, ONLY: READGRAPH,IT,VTXWGTS_LOC,VSIZE_LOC  
USE PARPREP, ONLY: DECOMPOSE_PAR,PREPGRID,X,Y,DP,FREE_MEM
USE MESSENGER 
IMPLICIT NONE

!---------------------------------------------------------------------
!      P R O G R A M   D R I V E R
!---------------------------------------------------------------------
!  KJR, UND, CHL, 2016 
!  MANIPULATE VERTEX WEIGHTS TO REPRESENT COMPUTATIONAL IMBALANCE. 
!  REDECOMPOSE THE DOMAIN TO LOAD BALANCE USING PARMETIS API 
!  THAT USES THE UNIFIED REPARTITIONING ALGORITHM
!---------------------------------------------------------------------


CALL MSG_INIT !start up MPI
CALL READGRAPH 
DO IT = 1,1 !timestep loop 
  CALL DECOMPOSE_PAR(VTXWGTS_LOC,VSIZE_LOC) 
  VTWGTS = !manipulate the weights based on a function of x,y and time 
  CALL FREE_MEM 
ENDDO
CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
