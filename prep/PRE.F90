MODULE PRE
USE MESSENGER
IMPLICIT NONE 
!This module contains most of the global variables, read/write subroutines for
!the fort.14 and prepares the data for the call to PARMETIS by building all
!relevant information locally on each MYPROC.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(len=80)   :: dmy !garbage  
CHARACTER(len=*), PARAMETER ::FILEBASE = "MYPROC_"
CHARACTER(LEN=*), PARAMETER ::FILEEND  = ".14"
CHARACTER(LEN=20)           :: FILENAME !filename used to write the subdomain grids
INTEGER               :: NE, NP,NE_G,NP_G  !number of elements,nodes    
INTEGER,ALLOCATABLE   :: X_G(:),Y_G(:),DP_G(:)
INTEGER,ALLOCATABLE   :: X_LOC(:),Y_LOC(:),DP_LOC(:) 
INTEGER,ALLOCATABLE   :: VTXWGTS_LOC(:),VTXWGTS_G(:) !these are defined for the elements
INTEGER,ALLOCATABLE   :: VSIZES_LOC(:),VSIZES_G(:) 
INTEGER               :: CHUNK,CHUNK_NP,LEFTOVER  !num of local ele, num. of local nodes
INTEGER               :: NOPE,NETA,NBOU,NVEL !boundary information
INTEGER,ALLOCATABLE   :: IEL(:),NNEL(:,:),EIND(:),EPTR(:),ELMDIST(:) !local ele conn.
!
INTEGER               :: IT ! time step counter 

CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT CONNECTIVITY   
!---------------------------------------------------------------------
  SUBROUTINE READGRAPH                                      

    IMPLICIT NONE 

    INTEGER :: I,J,K,O,ITEMP !counters
    REAL(8) :: T1,T2

    OPEN(13, FILE='fort.14',STATUS='OLD')
    !Read title 
    READ(13,*) dmy
    !Read number of elements and nodes
    READ(13,*) NE,NP
     NE_G = NE 
     NP_G = NP
    !determine the chunk size or number of eles on each rank
    IF(MYPROC.NE.(NPROC-1)) THEN !if not the last PE 
        CHUNK = NE/NPROC
    ELSE !if the last PE 
        CHUNK = NE/NPROC
        LEFTOVER = NE - (CHUNK*NPROC)         
        CHUNK = CHUNK + LEFTOVER 
    ENDIF

    ALLOCATE(X_G(NP_G),Y_G(NP_G),DP_G(NP_G)) 
    ALLOCATE(IEL(CHUNK),NNEL(3,CHUNK),EIND(3*CHUNK),EPTR(CHUNK+1))

    !! Read nodal connectivity table
    IF(MYPROC.EQ.0) THEN 
      DO I = 1,NP 
        READ(13,*) J,X_G(I),Y_G(I),DP_G(I)
      ENDDO
    ELSE 
      DO I = 1,NP 
        READ(13,*)
      ENDDO 
    ENDIF
  CALL MPI_BCAST(X_G,NP_G,MPI_REAL,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST(Y_G,NP_G,MPI_REAL,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST(DP_G,NP_G,MPI_REAL,0,MPI_COMM_WORLD,IERR)
! Read element connectivity table in chunks
! Here we read in the fort.14 ele. connectivity table by skipping over the parts
! MYPROC isn't getting. 
IF(MYPROC.EQ.0) THEN 
  PRINT *, "READING IN THE ELE TABLE ..."
ENDIF
CALL CPU_TIME(T1)
   NNEL=0  
   EIND=0
   EPTR=0 
   IF(MYPROC.NE.(NPROC-1)) THEN  !if not the last PE
    K=1 
    O=1 
    DO I = 1,NE
      IF(I.GE.(MYPROC*CHUNK+1).AND.I.LE.((MYPROC*CHUNK)+CHUNK)) THEN 
        READ(13,*) IEL(K),ITEMP,NNEL(1,K),NNEL(2,K),NNEL(3,K)
        EIND(O)=NNEL(1,K)
        O=O+1 
        EIND(O)=NNEL(2,K)
        O=O+1  
        EIND(O)=NNEL(3,K)
        O=O+1 
        K = K + 1 
      ELSE 
        READ(13,*)
      ENDIF
    ENDDO
  ELSE
   K=1
   O=1
    DO I = 1,NE 
     IF(I.GE.MYPROC*(CHUNK-LEFTOVER)+1) THEN 
       READ(13,*) IEL(K),ITEMP,NNEL(1,K),NNEL(2,K),NNEL(3,K)
       EIND(O)=NNEL(1,K) 
       O=O+1 
       EIND(O)=NNEL(2,K) 
       O=O+1 
       EIND(O)=NNEL(3,K) 
       O=O+1 
       K = K + 1
     ELSE 
       READ(13,*)
     ENDIF
    ENDDO
  ENDIF
  CLOSE(13)
CALL CPU_TIME(T2)
IF(MYPROC.EQ.0) THEN 
  PRINT *, "FINISHED READING IN THE ELE TABLE IN ",T2-T1
ENDIF

  K=1
   DO I = 1,SIZE(EIND)+1,3 !eles are triangular  
     EPTR(K)=I
     K=K+1 
   ENDDO

   ALLOCATE(ELMDIST(NPROC+1)) !same as vtxdist since dual graph 
   ELMDIST=0 
   DO I = 0,NPROC !the ele numbers on each rank
     ELMDIST(I+1)=(I*(CHUNK-LEFTOVER))+1
   ENDDO
   ELMDIST(NPROC+1)=NE_G+1 

! read in vertex weights from input file in working dir
  ALLOCATE(VTXWGTS_G(NE_G))
  ALLOCATE(VTXWGTS_LOC(CHUNK)) ! initially vtwgts is size chun
! IDEA IS THAT FILE I/O IS SLOWER THAN NETWORK
IF(MYPROC.EQ.0) THEN 
  OPEN(13,file='VW.txt',STATUS='old') !open file containing vertex weights here 
  DO I = 1,NE_G 
    READ(13,*) VTXWGTS_G(I)
  ENDDO
  CLOSE(13) 
ENDIF 
! broadcast the global vertex weights to all ranks
  CALL MPI_BCAST(VTXWGTS_G,NE_G,MPI_INT,0,MPI_COMM_WORLD,IERR)
!localize the vertex weights on each rank for the naive decomp, this is only done for IT==1 
J=1
DO I = 1,NE_G
  IF(I.GE.ELMDIST(MYPROC+1).and.I.LT.ELMDIST(MYPROC+2)) THEN 
     VTXWGTS_LOC(J)=VTXWGTS_G(I)
     J=J+1
  ELSE  
  ENDIF
ENDDO  
ALLOCATE(VSIZES_LOC(CHUNK)) 
VSIZES_LOC = 1 !this never changes so no need to localize YET 


!#ifdef DEBUG 
!!have each processor write its local grid so we can check for connectivity
!!problems.
!     WRITE(FILENAME,'(A,I3,A)') 'MYPROC_',MYPROC,'.14'
!     OPEN(UNIT=MYPROC+105,FILE=FILENAME)
!!write the title 
!     WRITE(MYPROC+105,85)  dmy
!!write the specs of the grid
!     WRITE(MYPROC+105,100) CHUNK,CHUNK_NP
!     !write nodes
!     DO I = 1,CHUNK_NP 
!       WRITE(MYPROC+105,95) I,X(I),Y(I),DP(I)
!     ENDDO
!     !write ele 
!     DO I = 1,CHUNK 
!       WRITE(MYPROC+105,90) I,3,NNEL_L(1,I),NNEL_L(2,I),NNEL_L(3,I)
!     ENDDO
!     !write boundary information to EOF
!      WRITE(MYPROC+105,80) NOPE! num open bou 
!      WRITE(MYPROC+105,75) NETA ! num of open bou. nodes 
!      WRITE(MYPROC+105,70) NBOU ! num of land bou 
!      WRITE(MYPROC+105,65) NVEL ! num of land bou. nodes    
!!
!100    FORMAT(2I10) !header 
!95     FORMAT(1I10,f10.5,f10.5,f10.5) !nodal connectivity table
!90     FORMAT(5I10) !element connectivity table
!85     FORMAT(A4)  !agrid
!!boundary information 
!80     FORMAT(1I5,"=Number of open boundaries") 
!75     FORMAT(1I5,"=Total number of open boundary nodes.") 
!70     FORMAT(1I5,"=Number of land boundaries.") 
!65     FORMAT(1I5,"=Total number of land boundary nodes.")
!     CLOSE(MYPROC+105)
!#endif DEBUG
  END SUBROUTINE READGRAPH
!---------------------------------------------------------------------------
!                       S U B R O U T I N E    S O R T 
!---------------------------------------------------------------------------
!  Sorts array RA of length N into ascending order using Heapsort algorithm.
!  N is input; RA is replaced on its output by its sorted rearrangement.
!  Ref: Numerical Recipes
!---------------------------------------------------------------------------
!
  SUBROUTINE SORT(N,RA)
      IMPLICIT NONE
      INTEGER N, L, IR, RRA, I, J
      INTEGER RA(N)

      L = N/2 + 1
      IR = N
10    CONTINUE
      IF (L.GT.1)THEN
        L=L-1
        RRA = RA(L)
      ELSE
        RRA=RA(IR)
        RA(IR)=RA(1)
        IR=IR-1
        IF (IR.EQ.1) THEN
          RA(1)=RRA
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
20    IF (J.LE.IR) THEN
        IF (J.LT.IR) THEN
          IF(RA(J).LT.RA(J+1)) J=J+1
        ENDIF
        IF (RRA.LT.RA(J)) THEN
          RA(I)=RA(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
        GO TO 20
      ENDIF
      RA(I)=RRA
      GO TO 10
   END SUBROUTINE SORT 
END MODULE PRE
