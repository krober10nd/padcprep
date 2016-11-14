MODULE PRE
USE MESSENGER
IMPLICIT NONE 
!This module contains all the global variables, read subroutines for
!the graph, and the metis libaries.  
CHARACTER(len=80)   :: dmy !garbage  
CHARACTER(len=*), PARAMETER ::FILEBASE = "MYPROC_"
CHARACTER(LEN=*), PARAMETER ::FILEEND  = ".14"
CHARACTER(LEN=20) :: FILENAME
INTEGER             :: NE, NP  !number of elements,nodes    
INTEGER             :: CHUNK,CHUNK_NP,LEFTOVER 
INTEGER,ALLOCATABLE :: NNEL(:,:) !element connectivity 
INTEGER,ALLOCATABLE :: ITVECT1(:)!vector of node numbers
INTEGER,ALLOCATABLE :: DMY_UQNNUM(:),UQNNUM(:) !unique node numbers contained in elements
REAL(8),ALLOCATABLE :: X(:),Y(:),X_G(:),Y_G(:) !position     
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT CONNECTIVITY   
!---------------------------------------------------------------------
  SUBROUTINE READGRAPH                                      
    IMPLICIT NONE 
    INTEGER :: I,J,K,ITEMP !counters
!
    OPEN(13, FILE='fort.14',STATUS='OLD')
    !--Read title 
    READ(13,*) dmy
    !--Read number of elements and nodes
    READ(13,*) NE,NP
   
IF(MYPROC.NE.(NPROC-1)) THEN !if not the last processor 
    CHUNK = NE/NPROC
ELSE 
    CHUNK = NE/NPROC
    LEFTOVER = NE - CHUNK*NPROC         
    CHUNK = CHUNK + LEFTOVER 
ENDIF

IF(MYPROC.EQ.1) THEN 
  PRINT *, "USING THIS MANY CORES",NPROC
  PRINT *, "NUMBER OF ELEMENTS:",NE 
  PRINT *, "NUMBER OF NODES:",NP
  PRINT *, "CHUNK SIZE",CHUNK
  PRINT *, "LEFTOVER",leftover
ENDIF     
! I choose to allocate 6*CHUNK for DMY_UQNNUM but that is arbitrary
 ALLOCATE(NNEL(3,CHUNK),ITVECT1(3*CHUNK),DMY_UQNNUM(6*CHUNK))
 ! skip nodes the first time, we'll come back
    DO I = 1,NP 
        READ(13,*)
    ENDDO
    !--Read Element Connectivity Table in chunks
   NNEL=0   
   IF(MYPROC.NE.(NPROC-1)) THEN 
    K=1
    DO I = 1,NE
      IF(I.GE.(MYPROC*CHUNK+1).AND.I.LE.((MYPROC*CHUNK)+CHUNK)) THEN 
        READ(13,*) J,ITEMP,NNEL(1,K),NNEL(2,K),NNEL(3,K)
        K = K + 1 
      ELSE 
        READ(13,*)
      ENDIF
     ENDDO
   ELSE
   K=1
     DO I = 1,NE 
      IF(I.GE.MYPROC*(CHUNK-LEFTOVER)+1) THEN 
         READ(13,*) J,ITEMP,NNEL(1,K),NNEL(2,K),NNEL(3,K)
         K = K + 1
      ELSE 
         READ(13,*)
      ENDIF
     ENDDO
   ENDIF
   CLOSE(13)
!step 1. find the node numbers associated with the elements that are on myproc
!step 2. sort them in ascending order 
!step 3. determine only the unique node numbers
   !reshape NNEL to vector 
    ITVECT1=0
    ITVECT1=RESHAPE(NNEL,(/CHUNK*3/))
    CALL SORT(CHUNK*3,ITVECT1)
    DMY_UQNNUM(:)=0
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
! have to handle the last entry 
    IF(ITVECT1(SIZE(ITVECT1)).NE.ITVECT1(SIZE(ITVECT1)-1).AND.ITVECT1(SIZE(ITVECT1)).NE.0) THEN 
     DMY_UQNNUM(J)=ITVECT1(SIZE(ITVECT1))
     CHUNK_NP=J
    ELSE 
     CHUNK_NP=CHUNK_NP-1
    ENDIF
    ALLOCATE(UQNNUM(CHUNK_NP))
    UQNNUM=0
    K=1 !only keep non zero entries 
    DO I = 1,CHUNK_NP 
        UQNNUM(K)=DMY_UQNNUM(I) 
        K=K+1
    ENDDO 
!now read in the node table associated with those elements 
!rewind the file to the top 
   ALLOCATE(X(CHUNK_NP),Y(CHUNK_NP))
   ALLOCATE(X_G(NP),Y_G(NP))
!
   OPEN(13,FILE='fort.14',STATUS='OLD')
   READ(13,*) 
   READ(13,*) 
     DO I = 1,NP 
       READ(13,*) ITEMP,X_G(I),Y_G(I)
     ENDDO 
   CLOSE(13)
   
   
   DO I = 1,CHUNK_NP 
     X(I)=X_G(UQNNUM(I))
     Y(I)=Y_G(UQNNUM(I))
     !PRINT *,I 
   ENDDO 
!have each processor write its local grid so we can check it's right
     WRITE(FILENAME,'(A,I3,A)') 'MYPROC_',MYPROC,'.14'
     OPEN(UNIT=MYPROC+105,FILE=FILENAME)
     WRITE(MYPROC+105,100) CHUNK,CHUNK_NP
     DO I = 1,CHUNK_NP 
       WRITE(MYPROC+105,95) I,X(I),Y(I)
     ENDDO 
     DO I = 1,CHUNK 
       WRITE(MYPROC+105,90) I,3,NNEL(1,I),NNEL(2,I),NNEL(3,I)
     ENDDO
100    FORMAT(2I10)
95     FORMAT(1I10,f10.5,f10.5)
90     FORMAT(5I10)
     CLOSE(MYPROC+105)
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
