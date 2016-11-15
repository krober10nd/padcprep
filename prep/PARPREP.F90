MODULE PARPREP
USE MESSENGER  
IMPLICIT NONE 

   INTEGER,ALLOCATABLE    :: CO_NODES(:,:),XADJ(:),LOCXADJ(:),ADJNCY(:),RADJNCY(:),VWGTS(:),EWGTS(:)  
CONTAINS 

!---------------------------------------------------------------------
!      S U B R O U T I N E   B U I L D 
!---------------------------------------------------------------------
!  ROUTINE TO BUILD XADJ AND ADJ 
! **** MOVE THIS PRE? AND MAKE THIS MODULE ONLY FOR CALLING PARMETIS??
!---------------------------------------------------------------------
SUBROUTINE BUILD 
USE PRE, only : NNEL_L,UQNNUM_L,X,Y,CHUNK,CHUNK_NP,LEFTOVER,g2l,SORT
   IMPLICIT NONE 
   INTEGER                :: NE,NP,INODE,JNODE,K,J,IEL,ITOT
   INTEGER                :: MNED,NCOUNT,MNEDLOC,NEDGETOT
   INTEGER,ALLOCATABLE    :: NEDLOC(:)   ! number of edges on each node 
   INTEGER,ALLOCATABLE    :: NEDGES(:)   ! also number of edges on each node,
   INTEGER,ALLOCATABLE    :: ITVECT1(:),ITVECT2(:)!temp vectors used for sorting
   LOGICAL                :: FOUND,SYMMETRIC 
   INTEGER                :: VTXDIST(NPROC+1)
!
!******NP and NE are reassigned to the number of local nodes and elements.
!  adj is built using local node numbers i.e., 1 : NP so this much be switched
!  back to global node numbers at the end before passing to parmetis

   NP = CHUNK_NP
   NE = CHUNK 
   ALLOCATE ( NEDGES(NP), NEDLOC(NP) )
   ALLOCATE ( XADJ(NP+1), VWGTS(NP)  ) 
!-------------------------------------------------------------
!  COMPUTES THE TOTAL NUMBER OF EDGES        -->    MNED
!  AND THE MAX NUMBER OF EDGES FOR ANY NODE  -->    MNEDLOC
!-------------------------------------------------------------
! build NEDLOC... the number of nodes connected locally to all other nodes       
        NEDLOC = 0
        MNED = 0
        DO J=1,3
           DO IEL=1,NE !loop over ele's part of sections
              INODE  = NNEL_L(J,IEL) !this retrieves the node numbers
              NCOUNT = NEDLOC(INODE) + 2 !since undirected, we represent edge twice
              MNED   = MNED + 2
              NEDLOC(INODE) = NCOUNT
           ENDDO
        ENDDO
!
       MNEDLOC = 0
        DO INODE=1,NP
           IF (NEDLOC(INODE).GE.MNEDLOC) THEN 
              MNEDLOC = NEDLOC(INODE)
           ENDIF  
       ENDDO

!        print *, "total number of edges = ", MNED,"ON MYPROC",MYPROC
!        print *, "maximum co-nodes for any node = ", MNEDLOC,"ON MYPROC",MYPROC
     
        ALLOCATE ( ADJNCY(MNED), EWGTS(MNED) )
        ALLOCATE ( CO_NODES(MNEDLOC,NP) )
! 
!-------------------------------------------------------------
!--COMPUTE CO_NODES LISTS AND NUMBER OF EDGES CONTAINING A NODE
!-------------------------------------------------------------
!
           NEDGES = 0
        DO IEL=1, NE 
           INODE = NNEL_L(1,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEL_L(2,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEL_L(3,IEL)
           NCOUNT = NEDGES(INODE) + 2
           NEDGES(INODE) = NCOUNT
        ENDDO
!
        DO IEL=1, NE
           INODE = NNEL_L(2,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEL_L(3,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEL_L(1,IEL)
           NCOUNT = NEDGES(INODE) + 2
           NEDGES(INODE) = NCOUNT
        ENDDO
!
        DO IEL=1, NE
           INODE = NNEL_L(3,IEL)
           CO_NODES(NEDGES(INODE)+1,INODE) = NNEL_L(1,IEL)
           CO_NODES(NEDGES(INODE)+2,INODE) = NNEL_L(2,IEL)
           NCOUNT = NEDGES(INODE) + 2
           NEDGES(INODE) = NCOUNT
        ENDDO
!
!-------------------------------------------------------------
!  REMOVE REDUNDANCY IN NODE LISTS
!-------------------------------------------------------------
! Allocate these sorting vectors based on the max. number of edges that
! intersect a node on the graph. 
! This avoids the issue when NP is on the order of the degree of the vertices.
   ALLOCATE ( ITVECT1(MNED),ITVECT2(MNED))
         NEDGETOT = 0           
         DO INODE = 1,NP   
           DO J=1, NEDGES(INODE)
              ITVECT1(J) = CO_NODES(J,INODE) !this is a temp. vec that stores
                                             !the nnum that are connected to INODE
          ENDDO
!
           IF (NEDGES(INODE).GT.1) THEN !if connected to more than one node
             NCOUNT = NEDGES(INODE) !number of nodes connectd to
              CALL SORT(NCOUNT,ITVECT1) !sort it in ascending order
              JNODE = ITVECT1(1) !first node 
              CO_NODES(1,INODE) = JNODE !sort CO-NODES so all the columns have
                                        !asending nnum order
              NCOUNT = 1
              DO J=2, NEDGES(INODE) !only retain unique nnum
                 IF (ITVECT1(J).NE.JNODE) THEN
                   NCOUNT = NCOUNT + 1
                   JNODE = ITVECT1(J)
                   CO_NODES(NCOUNT,INODE) = JNODE
                 ENDIF
              ENDDO
            ELSE 
             PRINT *, "node = ",INODE," is isolated on MYPROC = ", MYPROC
             stop
           ENDIF
          NEDGES(INODE) = NCOUNT
          NEDGETOT = NEDGETOT + NCOUNT
           IF (NEDGES(INODE) == 0) THEN
             PRINT *, "inode = ", INODE, " belongs to no edges on MYPROC = ", MYPROC
            ENDIF
         ENDDO
         NEDGETOT = NEDGETOT/2

!  CHECK THAT ADJACENCY MATRIX IS SYMMETRIC 

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
           print *, "node ",inode," adjacent to ",jnode," but not visa-versa on MYPROC", MYPROC
         ENDIF
       ENDDO
       ENDDO
       IF (.NOT. SYMMETRIC) THEN
          print *, "bad adjacency matrix: not symmetric on MYPROC", MYPROC
          stop
       ENDIF


! COMPUTE WEIGHTS OF THE GRAPH VERTICES
        DO INODE = 1,NP   
           VWGTS(INODE) = NEDGES(INODE)
        ENDDO
! COMPUTE ADJACENCY LIST OF GRAPH AND ITS EDGE WEIGHTS
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
          PRINT *, "FINISHED BUILDING LOCAL XADJ AND ADJ ON MYPROC", MYPROC
      RETURN
      END SUBROUTINE BUILD


!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  P A R 
!---------------------------------------------------------------------
!  ROUTINE TO CALL PARMETIS 4.0 AFTER THE CALL TO DECOMPOSE_SERIAL
!---------------------------------------------------------------------
!
      SUBROUTINE DECOMPOSE_PAR
         IMPLICIT NONE 
     

      RETURN
      END SUBROUTINE DECOMPOSE_PAR 
END MODULE PARPREP  
