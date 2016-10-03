MODULE PRE
USE MESSENGER
IMPLICIT NONE 
!This module contains all the global variables, read subroutines for
!the graph, and the metis libaries.  
CHARACTER(len=80)   :: AGRAPH,dmy !graph name, garbage  
INTEGER             :: NE, NP  !number of elements,nodes    
INTEGER,ALLOCATABLE :: NNELG(:,:) !global element connectivity 
INTEGER,ALLOCATABLE :: NNUM(:) !node numbers
REAL(8),ALLOCATABLE :: X(:),Y(:) !position     
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT, AND BOUNDARY CONNECTIVITY   
!---------------------------------------------------------------------
  SUBROUTINE READGRAPH                                      
    IMPLICIT NONE 
    INTEGER             :: I,J,ITEMP !counters
IF(MYPROC.EQ.0) THEN 
     WRITE(*,*)'ENTER THE NAME OF THE GRAPH'
     READ(*,'(A)')AGRAPH 
    OPEN(13, FILE=AGRAPH,STATUS='OLD')
    !--Read title 
    READ(13,*) dmy
    !--Read number of elements and nodes
    READ(13,*) NE,NP
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "NUMBER OF ELEMENTS:",NE 
    PRINT *, "NUMBER OF NODES:",NP
  ENDIF 
    !
    ALLOCATE(X(NP),Y(NP),NNUM(NP))
    ALLOCATE(NNELG(3,NE))
    !
    !--Read nodal connectivity and bathymetry
    !
    I = 1 
    J = 1  
    DO I = 1,NP
       READ(13,*) J,X(I),Y(I)
       NNUM(I) = J
       IF (J.NE.I) THEN
          print *, I,J
          STOP 'Node Numbering not in Sequential Order'
       ENDIF
    ENDDO
    !
    !--Read Element Connectivity Table
    !
    DO I = 1,NE
       READ(13,*) J,ITEMP,NNELG(1,I),NNELG(2,I),NNELG(3,I)
       IF (J.NE.I) THEN
         print *, I,J
         STOP 'Element Numbering not Sequential'
       ENDIF
    ENDDO
  CLOSE(13)
80    FORMAT(A95)
ENDIF
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
