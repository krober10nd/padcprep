      PROGRAM DynamicExchange
      IMPLICIT NONE 
!******************************************************
! This is an experiment for testing how to exchange information
! dynamically between two domains on a triangular finite element 
! mesh.
! kjr,started 2016 09, 28
!*****************************************************
      CHARACTER(len=80) :: AGRID !grid name 
      INTEGER           :: NELG, NNODG !number of elements and nodes     
! read the mesh 
      OPEN(unit=14, file="fort.14")
      READ(14,80) AGRID

! Read elements and nodes
      READ(14,*) NELG,NNODG
      
      PRINT *, "THE GRID IS:",AGRID
      PRINT *, "NUMBER OF ELEMENTS:",NELG 
      PRINT *, "NUMBER OF NODES:",NNODG
      
      
      CLOSE(14)
80    FORMAT(A95)

      END PROGRAM DynamicExchange  
