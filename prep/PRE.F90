MODULE PRE
USE MESSENGER
IMPLICIT NONE 
!This module contains most of the global variables, read/write subroutines for
!the fort.14 and prepares the data for the call to PARMETIS by building all
!relevant information locally on each MYPROC.  
!I convert to local node numbers from global nnum to create a dense xadj and adj
CHARACTER(len=80)   :: dmy !garbage  
CHARACTER(len=*), PARAMETER ::FILEBASE = "MYPROC_"
CHARACTER(LEN=*), PARAMETER ::FILEEND  = ".14"
CHARACTER(LEN=20) :: FILENAME !filename used to write the subdomain grids
INTEGER*8               :: NE, NP,NE_G,NP_G  !number of elements,nodes    
INTEGER*8               :: CHUNK,CHUNK_NP,LEFTOVER  !num of local ele, num. of local nodes
INTEGER*8             :: NOPE,NETA,NBOU,NVEL !boundary information
INTEGER*8,ALLOCATABLE :: IEL(:),NNEL(:,:),EIND(:),EPTR(:),ELMDIST_L(:),ELMDIST(:) !local ele conn.
INTEGER*8,ALLOCATABLE :: DMY_UQNNUM(:),UQNNUM(:),UQNNUM_G(:)
REAL(8),ALLOCATABLE :: X(:),Y(:),DP(:),X_G(:),Y_G(:),DP_G(:) !local position, local bathy,etc.
!_G denotes a global variable, otherwise it is local i.e., only for MYPROC but
!still using global node numbers.  
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT CONNECTIVITY   
!---------------------------------------------------------------------
  SUBROUTINE READGRAPH                                      
    IMPLICIT NONE 
    INTEGER :: I,J,K,O,ITEMP !counters
    INTEGER,ALLOCATABLE :: ITVECT1(:)!vector of node numbers
    OPEN(13, FILE='fort.14',STATUS='OLD')
!Read title 
    READ(13,*) dmy
!Read number of elements and nodes
    READ(13,*) NE,NP
     NE_G = NE 
     NP_G = NP
!determine the chunk size
    IF(MYPROC.NE.(NPROC-1)) THEN !if not the last PE 
        CHUNK = NE/NPROC
    ELSE !if the last PE 
        CHUNK = NE/NPROC
        LEFTOVER = NE - (CHUNK*NPROC)         
        CHUNK = CHUNK + LEFTOVER 
    ENDIF
!     IF(MYPROC.EQ.(NPROC-1)) THEN 
!       PRINT *, "USING THIS MANY CORES",NPROC
!       PRINT *, "NUMBER OF ELEMENTS:",  NE 
!       PRINT *, "NUMBER OF NODES:",     NP
!       PRINT *, "LOCAL ELE. CHUNK SIZE",CHUNK,"ON MYPROC", MYPROC
!       PRINT *, "LEFTOVER OF",LEFTOVER,"ON MYPROC",MYPROC
!    ENDIF  
!allocate some global arrays
    ALLOCATE(UQNNUM_G(NP),X_G(NP),Y_G(NP),DP_G(NP)) 
    ALLOCATE(IEL(CHUNK),NNEL(3,CHUNK),ITVECT1(3*CHUNK),EIND(3*CHUNK),EPTR(CHUNK+1))
! Read nodal connectivity table
     DO I = 1,NP 
       READ(13,*) UQNNUM_G(I),X_G(I),Y_G(I),DP_G(I)
     ENDDO 
! Read element connectivity table in chunks
! Here we read in the fort.14 ele. connectivity table by skipping over the parts
! MYPROC isn't getting. 
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
   K=1
   DO I = 1,SIZE(EIND)+1,3 
     EPTR(K)=I
     K=K+1 
   ENDDO

   ALLOCATE(ELMDIST(NPROC+1)) 
   ELMDIST=0 
   DO I = 0,NPROC
     ELMDIST(I+1)=(I*(CHUNK-LEFTOVER))+1
   ENDDO
   !PRINT *, ELMDIST_L," ON MYPROC ", MYPROC
   !CALL MPI_ALLREDUCE(ELMDIST_L,ELMDIST,NPROC+1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,IERR)
   ELMDIST(NPROC+1)=NE_G+1 
   !PRINT *, ELMDIST 
! Determine the unique node numbers on MYPROC
! I chose to allocate 6*CHUNK for DMY_UQNNUM but that is arbitrary
! We just need a large enough array to make sure we don't exceed
! the bounds.
   ALLOCATE(DMY_UQNNUM(6*CHUNK))
   DMY_UQNNUM(:)=0
! Steps necessary to determine the local node numbers associated with chunk of
! elements you just read into MYPROC.
!step 1. find the node numbers associated with the elements that are on myproc.
!step 2. sort them in ascending order. 
!step 3. determine only the unique node numbers from that sorted array.
    ITVECT1=0
    ITVECT1=RESHAPE(NNEL,(/CHUNK*3/))
    CALL SORT(CHUNK*3,ITVECT1)
!only keep unique node numbers 
    J=1
    DO I = 1,SIZE(ITVECT1)-1
      IF(ITVECT1(I).EQ.ITVECT1(I+1).AND.ITVECT1(I).NE.0) THEN
          !skip
      ELSE !store
        DMY_UQNNUM(J)=ITVECT1(I)
        J=J+1
      ENDIF
    ENDDO
! Handle the last position of the sorted array since we stop one before the end.
     DMY_UQNNUM(J)=ITVECT1(SIZE(ITVECT1))
     CHUNK_NP=J !This is the number of unique nodes on MYPROC.
! chunk_np is the number of nodes on each partition      
    ALLOCATE(UQNNUM(CHUNK_NP))
    ALLOCATE(X(CHUNK_NP),Y(CHUNK_NP),DP(CHUNK_NP))  
!uqnnum is the the unique node numbers on each partition 
    UQNNUM=0
    K=1 !only keep non zero entries 
    DO I = 1,CHUNK_NP 
      UQNNUM(K)=DMY_UQNNUM(I) 
      K=K+1
    ENDDO 
   
   !get the local nodes position and bathy on each PE
    DO I = 1,CHUNK_NP
     X(I)=X_G(UQNNUM(I))
     Y(I)=Y_G(UQNNUM(I))
     DP(I)=DP_G(UQNNUM(I))    
   ENDDO

! these are assumed to be zero for now
     NOPE = 0 
     NETA = 0 
     NBOU = 0 
     NVEL = 0  
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
    RETURN 
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
