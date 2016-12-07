PROGRAM Driver

USE PRE, ONLY: READGRAPH,IT,VTXWGTS_LOC,VSIZES_LOC,X_G,Y_G,DP_G,CHUNK
USE PARPREP, ONLY: DECOMPOSE_PAR,EIND
USE MESSENGER 
IMPLICIT NONE

REAL(8) :: X1,X2,X3,AVGX
REAL(8) :: XFRONT,C
INTEGER :: NOIT,O,I
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
IF(MYPROC.EQ.0) THEN 
  PRINT *,"Enter no. of time steps"
  read *, NOIT
  PRINT *,"Enter speed of front" 
  read *, C 
  PRINT *, "Enter X location of front" 
  read *, XFRONT
ENDIF
CALL MPI_BCAST(NOIT,1,MPI_INT,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(C,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(XFRONT,1,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
DO IT = 1,NOIT !timestep loop 
  CALL DECOMPOSE_PAR 
   !manipulate the weights based on a function of x,y and time 
  O=1
  DO I = 1,CHUNK
    ! average elemental position
    X1=X_G(EIND(O))
    O=O+1
    X2=X_G(EIND(O))
    O=O+1
    X3=X_G(EIND(O))
    O=O+1
    AVGX=(X1+X2+X3)/3D0
    IF(AVGX.GT.XFRONT) THEN 
       VTXWGTS_LOC(I) = 0
    ELSE 
       VTXWGTS_LOC(I) = 1
    ENDIF
   ENDDO
   IF(MYPROC.EQ.0) THEN 
     PRINT *, "THE FRONT IS AT ",XFRONT 
   ENDIF
   XFRONT=XFRONT-C*IT
ENDDO
CALL MSG_FINI !shut MPI down  

END PROGRAM Driver
