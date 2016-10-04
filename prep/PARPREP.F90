MODULE PARPREP
USE PRE
USE MESSENGER  
IMPLICIT NONE 

   INTEGER,ALLOCATABLE    :: CO_NODES(:,:),XADJ(:),RXADJ(:),ADJNCY(:),VWGTS(:),EWGTS(:)  
CONTAINS 

!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  S E R I A L 
!---------------------------------------------------------------------
!  ROUTINE TO BUILD XADJ AND ADJ 
!  DEVELOP XADJ AND ADJ FOR ENTIRE GRID 
!---------------------------------------------------------------------
SUBROUTINE DECOMPOSE_SERIAL 
   IMPLICIT NONE 
!
   INTEGER                :: INODE,JNODE,K,J,IEL,ITOT  !COUNTERS
   INTEGER                :: MNED,NCOUNT  
   INTEGER                :: MNEDLOC     ! maximum number of edges
   INTEGER                :: NEDGETOT
!
   INTEGER,ALLOCATABLE    :: NEDLOC(:)   ! number of edges on each node 
   INTEGER,ALLOCATABLE    :: NEDGES(:)   ! also number of edges on each node,
   INTEGER,ALLOCATABLE    :: ITVECT1(:),ITVECT2(:)
   LOGICAL                :: FOUND,SYMMETRIC 
!
   INTEGER                :: CHUNK,LEFTOVER       
   INTEGER                :: SCHUNK(0:NPROC-1),DISPS(0:NPROC-1)
!
        ALLOCATE( NEDGES(NP), NEDLOC(NP) )
        ALLOCATE ( ITVECT1(NP),ITVECT2(NP) )
        ALLOCATE ( XADJ(NP+1), VWGTS(NP)) 
!
!-------------------------------------------------------------
!  COMPUTES THE TOTAL NUMBER OF EDGES        -->    MNED
!  AND THE MAX NUMBER OF EDGES FOR ANY NODE  -->    MNEDLOC
!-------------------------------------------------------------
!make sure all counter variables are assigned 
        MNED    =   0  
        NCOUNT =    0 
        IEL    =    1 
        INODE  =    1
        J      =    1
! 
        MNED = 0
        DO INODE = 1,NP
           NEDLOC(INODE) = 0
        ENDDO
!
        DO J=1,3
           DO IEL=1, NE
              INODE = NNELG(J,IEL)
              NCOUNT = NEDLOC(INODE) + 2
              MNED = MNED + 2
              NEDLOC(INODE) = NCOUNT
           ENDDO
        ENDDO
!
        MNEDLOC = 0
        DO INODE=1, NP
           IF (NEDLOC(INODE).GE. MNEDLOC) THEN 
              MNEDLOC = NEDLOC(INODE)
           ENDIF  
       ENDDO

IF(MYPROC.EQ.0) THEN
        print *, "total number of edges = ", MNED
        print *, "maximum co-nodes for any node = ", MNEDLOC
ENDIF 
        ALLOCATE ( ADJNCY(MNED), EWGTS(MNED) )
        ALLOCATE ( CO_NODES(MNEDLOC,NP) )
! 
!-------------------------------------------------------------
!--COMPUTE CO_NODES LISTS AND NUMBER OF EDGES CONTAINING A NODE
!-------------------------------------------------------------
!
        DO INODE = 1,NP
           NEDGES(INODE) = 0
        ENDDO
!
        DO IEL=1, NE 
           INODE = NNELG(1,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNELG(2,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNELG(3,IEL)
           NCOUNT = NEDGES(INODE) + 2
           NEDGES(INODE) = NCOUNT
        ENDDO
!
        DO IEL=1, NE  
           INODE = NNELG(2,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNELG(3,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNELG(1,IEL)
           NCOUNT = NEDGES(INODE) + 2
           NEDGES(INODE) = NCOUNT
        ENDDO
!
        DO IEL=1, NE 
           INODE = NNELG(3,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNELG(1,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNELG(2,IEL)
           NCOUNT = NEDGES(INODE) + 2
           NEDGES(INODE) = NCOUNT
        ENDDO
!
!-------------------------------------------------------------
!  REMOVE REDUNDANCY IN NODE LISTS
!-------------------------------------------------------------
!
         NEDGETOT = 0           !  This will be twice number of edges
         DO INODE = 1,NP   
           DO J=1, NEDGES(INODE)
              ITVECT1(J) = CO_NODES(J,INODE)
           ENDDO
!
           IF (NEDGES(INODE).GT.1) THEN
             NCOUNT = NEDGES(INODE)
              CALL SORT(NCOUNT,ITVECT1)
              JNODE = ITVECT1(1)
              CO_NODES(1,INODE) = JNODE
              NCOUNT = 1
              DO J=2, NEDGES(INODE)
                 IF (ITVECT1(J).NE.JNODE) THEN
                   NCOUNT = NCOUNT + 1
                   JNODE = ITVECT1(J)
                   CO_NODES(NCOUNT,INODE) = JNODE
                 ENDIF
              ENDDO
            ELSE
             IF(MYPROC.EQ.0) THEN 
              PRINT *, "node = ",INODE," is isolated"
              stop 'KEITH'
             ENDIF    
            ENDIF
            NEDGES(INODE) = NCOUNT
            NEDGETOT = NEDGETOT + NCOUNT
            IF (NEDGES(INODE) == 0) THEN
             IF(MYPROC.EQ.0) THEN
              PRINT *, "inode = ", INODE, " belongs to no edges"
              STOP 'KEITH'
             ENDIF
            ENDIF
         ENDDO
         NEDGETOT = NEDGETOT/2
        IF(MYPROC.EQ.0) THEN
         PRINT *, "edge count = ",NEDGETOT
        ENDIF
!
!  CHECK THAT ADJACENCY MATRIX IS SYMMETRIC 
!
       SYMMETRIC = .TRUE.
       DO INODE = 1, NP
       DO J = 1, NEDGES(INODE)
          JNODE = CO_NODES(J,INODE)
          FOUND = .FALSE.
          DO K= 1, NEDGES(JNODE)
             IF (CO_NODES(K,JNODE) == INODE) THEN
               FOUND = .TRUE.
               EXIT
             ENDIF
          ENDDO
          IF (.not. FOUND) THEN
            SYMMETRIC = .FALSE.
         IF(MYPROC.EQ.0) THEN 
           print *, "node ",inode," adjacent to ",jnode," but not visa-versa"
         ENDIF 
         ENDIF
       ENDDO
       ENDDO
       IF (.NOT. SYMMETRIC) THEN
         IF(MYPROC.EQ.0) THEN  
         stop 'bad adjacency matrix: not symmetric!'
         ENDIF
       ENDIF
!
!  COMPUTE WEIGHTS OF THE GRAPH VERTICES
!
        DO INODE = 1,NP   
           VWGTS(INODE) = NEDGES(INODE)
        ENDDO
!
!--COMPUTE ADJACENCY LIST OF GRAPH AND ITS EDGE WEIGHTS
!
        XADJ(1) = 1
        ITOT = 0
        DO INODE = 1,NP
        DO J = 1, NEDGES(INODE)
           ITOT = ITOT + 1
           JNODE = CO_NODES(J,INODE)
           ADJNCY(ITOT) = JNODE
           EWGTS(ITOT)  = (VWGTS(JNODE)+VWGTS(INODE))
        ENDDO
        XADJ(INODE+1) = ITOT+1
        ENDDO
!
! DUMP GRAPH TO A FILE FOR DEBUGGING
IF(MYPROC.EQ.0) THEN 
        OPEN(FILE='metis_graph.txt',UNIT=99)
        WRITE(99,100) NP, NEDGETOT, 11, 1
        DO INODE=1, NP
           WRITE(99,200) VWGTS(INODE),(CO_NODES(J,INODE), EWGTS(XADJ(INODE)+J-1),J=1,NEDGES(INODE))
        ENDDO
        CLOSE(99)
ENDIF       
! 
! SCATTERV XADJ TO ALL  PE'S
!
        CHUNK=NP/NPROC 
        LEFTOVER=MOD(NP,CHUNK) 
        IF(LEFTOVER.EQ.0) THEN
            DO J = 1,NPROC !NP/NPROC has zero modulus 
               SCHUNK(J-1)=CHUNK !equal chunks to scatter
            ENDDO
          ELSE
            DO J = 1,NPROC-1 !NP/NPROC has nonzero modulus
               SCHUNK(J-1)=CHUNK 
            ENDDO
             !make the last one a little longer 
               SCHUNK(NPROC-1)=CHUNK+LEFTOVER
        ENDIF 
            DO J = 1,NPROC 
              DISPS(J-1)=0 !displacements for scatterv, none
            ENDDO
        ALLOCATE(RXADJ(SCHUNK(MYPROC))) !receive buffer size 

        CALL MPI_SCATTERV(XADJ,SCHUNK,DISPS,MPI_INT,RXADJ,SCHUNK(MYPROC),MPI_INT,0,MPI_COMM_WORLD,IERR)
 
100    FORMAT(4I10)
200    FORMAT(100I10)
      
      RETURN
      END SUBROUTINE DECOMPOSE_SERIAL 


!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  P A R 
!---------------------------------------------------------------------
!  ROUTINE TO CALL PARMETIS 4.0 AFTER THE CALL TO DECOMPOSE_SERIAL
!---------------------------------------------------------------------
!
      SUBROUTINE DECOMPOSE_PAR
         IMPLICIT NONE 
     
      PRINT*,"MYPROC ",MYPROC," XADJ RECEIVED",SIZE(RXADJ(:))



      RETURN
      END SUBROUTINE DECOMPOSE_PAR 
END MODULE PARPREP  
