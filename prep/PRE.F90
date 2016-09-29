MODULE PRE
USE MESSENGER
IMPLICIT NONE 
!This module contains all the global variables, read subroutines for
!the graph, and the metis libaries.  
CHARACTER(len=80)   :: AGRAPH,dmy !graph name, garbage  
INTEGER             :: NE, NP !number of elements,nodes    
INTEGER,ALLOCATABLE :: NNELG(:,:) !global element connectivity 
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
    ALLOCATE(X(NP),Y(NP))
    ALLOCATE(NNELG(3,NE))
    !
    !--Read nodal connectivity and bathymetry
    !
    I = 1 
    J = 1  
    DO I = 1,NP
       READ(13,*) J,X(I),Y(I)
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
END MODULE PRE
