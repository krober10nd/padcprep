PROGRAM DynamicExchange

USE PRE, ONLY: X,Y,DP,NNEG,READGRAPH 
IMPLICIT NONE 
!******************************************************
! This is an experiment for testing how to exchange information
! dynamically between two domains on a triangular finite element 
! mesh.
! kjr,started 2016 09, 28
!*****************************************************
CALL READGRAPH  
END PROGRAM DynamicExchange  
