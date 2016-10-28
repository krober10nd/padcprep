MODULE PRE
USE MESSENGER
IMPLICIT NONE 
!This module contains all the global variables, read subroutines for
!the graph, and the metis libaries.  
CHARACTER(len=80)   :: AGRAPH,dmy !graph name, garbage  
INTEGER             :: NE, NP  !number of elements,nodes    
INTEGER             :: CHUNK 
INTEGER,ALLOCATABLE :: NNELG(:,:) !element connectivity 
INTEGER,ALLOCATABLE :: NNUM(:) !node numbers
INTEGER,ALLOCATABLE :: ITVECT1(:)!vector of elements numbers
INTEGER,ALLOCATABLE :: DMY_UQNNUM(:),UQNNUM(:) !unique node numbers contained in elements
REAL(8),ALLOCATABLE :: X(:),Y(:) !position     
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT CONNECTIVITY   
!---------------------------------------------------------------------
  SUBROUTINE READGRAPH                                      
    IMPLICIT NONE 
    INTEGER             :: I,J,ITEMP !counters
    OPEN(13, FILE='fort.14',STATUS='OLD')
    !--Read title 
    READ(13,*) dmy
    !--Read number of elements and nodes
    READ(13,*) NE,NP
    CHUNK = NE/NPROC
   !
    IF(MYPROC.EQ.0) THEN 
      PRINT *, "NUMBER OF ELEMENTS:",NE 
      PRINT *, "NUMBER OF NODES:",NP
      PRINT *, "CHUNK SIZE",CHUNK
    ENDIF     
    
IF(MYPROC.NE.(NPROC-1)) !if not the last processor 
    ALLOCATE(NNELG(3,CHUNK),ITVECT1(3*CHUNK),UQNNUM(6*CHUNK))
    UQNNUM(:)=0
    ! skip nodes the first time, we'll come back
    DO I = 1,NP 
        READ(13,*)
    ENDDO
    !
    !--Read Element Connectivity Table in chunks
    !
    DO I = 1,NE
     IF(I.GE.(MYPROC*CHUNK+1).AND.I.LE.(MYPROC*CHUNK)+CHUNK) THEN 
       READ(13,*) J,ITEMP,NNELG(1,I-(MYPROC*CHUNK)),NNELG(2,I-(MYPROC*CHUNK)),NNELG(3,I-(MYPROC*CHUNK))
     ELSE 
       READ(13,*)
     ENDIF
    ENDDO

!step 1. find the node numbers associated with the elements that are on myproc
!step 2. sort them in ascending order 
!step 3. determine only the unique node numbers
    IF(MYPROC.EQ.0) THEN 
      !reshape NNELG to vector 
      ITVECT1=RESHAPE(NNELG,(/CHUNK*3/))
      CALL SORT(CHUNK*3,ITVECT1)
    ENDIF
    !only keep unique node numbers 
    J=1
    DO I = 1,SIZE(ITVECT1)-1
        IF(I.EQ.1) THEN 
         DMY_UQNNUM(1)=ITVECT1(1)
        ENDIF 
        IF(ITVECT1(I).EQ.ITVECT1(I+1)) THEN
          !skip
        ELSE !store
          DMY_UQNNUM(J)=ITVECT1(I)
          J=J+1
        ENDIF
    ENDDO
    CLOSE(13)
! chop off the extra zeros in dmy_uqnnum 
! cycle if zero
    J=1
    DO I = 1,SIZE(DMY_UQNNUM) 
       IF(DMY_UQNNUM(I).EQ.0) CYCLE
        J=J+1
    ENDDO
    ALLOCATE(uqnnum(J)) !j is now # of non-zero enteries
    K=1 !only keep non zero enteries 
    DO I = 1,SIZE(DMY_UQNNUM) 
       IF(DMY_UQNNUM(I).EQ.0) CYCLE
        UQNNUM(K)=DMY_UQNNUM(I) 
         K=K+1
    ENDDO 


!TODO: handle the case when NE/NPROC has non-zero modulus.
!rewind the file to the top 
    OPEN(13,FILE='fort.14',STATUS='OLD')
    !--Read title 
    READ(13,*)
    !--Read number of elements and nodes
    READ(13,*) 

    ALLOCATE(X(NP),Y(NP),NNUM(NP))
    !
    !--Read nodal connectivity and bathymetry
    !
    I = 1 
    J = 1  
    DO I = 1,NP !read in only the sections of the graph you need 
       READ(13,*) J,X(I),Y(I)
       NNUM(I) = J
    ENDDO
ELSE !must be last processor, then the chunk size chunk + leftover






ENDIF !






80    FORMAT(A95)
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
