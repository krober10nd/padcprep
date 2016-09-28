MODULE PRE
IMPLICIT NONE 
!This module contains all the global variables and read subroutines for
!the graph. 
CHARACTER(len=80)   :: AGRAPH !graph name 
INTEGER             :: NELG, NNODG !number of elements and nodes     
INTEGER             :: I,J,ITEMP !counters
INTEGER,ALLOCATABLE :: NNEG(:,:) !global element connectivity 
REAL(8),ALLOCATABLE :: DP(:),X(:),Y(:) !bathymetry, lon and lat
      
CONTAINS 
  SUBROUTINE READGRAPH                                      
  !--Read the graph
  OPEN(unit=14, file="fort.14")
  READ(14,80) AGRAPH
  !--Read number of elements and nodes
  READ(14,*) NELG,NNODG
  PRINT *, "NUMBER OF ELEMENTS:",NELG 
  PRINT *, "NUMBER OF NODES:",NNODG
  ALLOCATE(DP(NNODG),X(NNODG),Y(NNODG))
  ALLOCATE(NNEG(3,NELG))
  !
  !--Read nodal connectivity and bathymetry (DP)
  ! 
  DO I = 1,NNODG
     READ(14,*) J,X(I),Y(I),DP(I)
     IF (J.NE.I) THEN
        print *, I,J
        STOP 'Node Numbering not in Sequential Order'
     ENDIF
  ENDDO
  !
  !--Read Element Connectivity Table
  !
  DO I = 1,NELG
     READ(14,*) J,ITEMP,NNEG(1,I),NNEG(2,I),NNEG(3,I)
     IF (J.NE.I) THEN
       print *, I,J
       STOP 'Element Numbering not Sequential'
     ENDIF
  ENDDO
  CLOSE(14)
  80    FORMAT(A95)
  RETURN 
  END SUBROUTINE READGRAPH
END MODULE PRE
