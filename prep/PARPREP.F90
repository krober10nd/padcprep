MODULE PARPREP
USE PRE
USE MESSENGER  
IMPLICIT NONE 

CONTAINS 

!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E 
!---------------------------------------------------------------------
!  ROUTINE TO DECOMPOSE A GRID USING PARMETIS 4.0  
!---------------------------------------------------------------------
SUBROUTINE DECOMPOSE 
   IMPLICIT NONE 
! 
   INTEGER                :: WGTFLAG, NUMFLAG, NPARTS
   INTEGER,ALLOCATABLE    :: CO_NODES(:,:),XADJG(:),ADJNCYG(:)  
   INTEGER,ALLOCATABLE    :: ITVECT1(:),ITVECT2(:)
!
   INTEGER                :: INODE,JNODE,J,IEL  !COUNTERS
   INTEGER                :: MNED,NCOUNT  
   INTEGER                :: MNEDLOC     ! maximum number of edges
   INTEGER                :: NEDGETOT
!
   INTEGER,ALLOCATABLE    :: NEDLOC(:)   ! number of edges on each node 
   INTEGER,ALLOCATABLE    :: NEDGES(:)   ! also number of edges on each node,
!                                        ! used in different places
   INTEGER,ALLOCATABLE    :: OPTIONS     ! options to pass  
   INTEGER,ALLOCATABLE    :: VTXDIST(:)  ! how the nodes are distributed across all PEs.
   INTEGER,ALLOCATABLE    :: XL(:),YL(:) ! nodal positions local to this PE 
!
   INTEGER                :: BGIND,ENDIND! here I split the nodes equally among PEs
   INTEGER                :: CHUNK       ! amount of nodes passed to each PE
   INTEGER                :: LEFTOVER    ! some vectors may not be of equal len
   INTEGER,ALLOCATABLE    :: NPLdmy(:),NPL(:)  ! vectors of local nodes on this PE
   INTEGER,ALLOCATABLE    :: REMAINDER(:)!used to deal with cases when NP/NPROC
!                                        !has modulus
!
! break apart vector nodes into chunks
   CHUNK=NP/NPROC ! integer division rounds down
   LEFTOVER=MOD(NP,NPROC) !modulus
ALLOCATE( NPLdmy(CHUNK) )
! nodes that are sent to this PE
   BGIND=(MYPROC*CHUNK)+1 
   ENDIND=(MYPROC+1)*CHUNK
   NPLdmy(:)=NP(BGIND:ENDIND)
! if NP/NPROC has a nonzero modulus, then send these nodes to the last PE (NPROC) 
IF(MYPROC.EQ.NPROC) THEN 
   IF(LEFTOVER.GE.0) THEN
      ALLOCATE( REMAINDER(LEFTOVER) )
      ALLOCATE( NPL(LEFTOVER+CHUNK) ) 
      REMAINDER(:)=NP( (CHUNK*MYPROC)+1:(CHUNK*MYPROC)+LEFTOVER )
      NPL = [NPLdmy,REMAINDER] 
   ENDIF
  ELSE !i.e., not the last PE (NPROC) 
   ALLOCATE(NPL(CHUNK))
   NPL = NPLdmy
ENDIF
 
   ALLOCATE( NEDGES(CHUNK), NEDLOC(CHUNK) )
   ALLOCATE( ITVECT1(CHUNK),ITVECT2(CHUNK) )

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
        DO J=1, 3
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


        print *, "total number of edges = ", MNED
        print *, "maximum co-nodes for any node = ", MNEDLOC
!
        ALLOCATE ( ADJNCYG(MNED) )
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
!        NEDGETOT = 0           !  This will be twice number of edges
!        DO INODE = 1,NP   
!           DO J=1, NEDGES(INODE)
!              ITVECT1(J) = CO_NODES(J,INODE)
!           ENDDO
!
!           IF (NEDGES(INODE).GT.1) THEN
!             NCOUNT = NEDGES(INODE)
!             CALL SORT(NCOUNT,ITVECT1)
!             JNODE = ITVECT1(1)
!             CO_NODES(1,INODE) = JNODE
!             NCOUNT = 1
!             DO J=2, NEDGES(INODE)
!                IF (ITVECT1(J).NE.JNODE) THEN
!                  NCOUNT = NCOUNT + 1
!                  JNODE = ITVECT1(J)
!                  CO_NODES(NCOUNT,INODE) = JNODE
!                ENDIF
!             ENDDO
!           ELSE
!             print *, "node = ",INODE," is isolated"
!             stop 'vic'
!           ENDIF
!           NEDGES(INODE) = NCOUNT
!           NEDGETOT = NEDGETOT + NCOUNT
!           if (nedges(inode) == 0) then
!             print *, "inode = ", inode, " belongs to no edges"
!             stop 'vic'
!           endif
!        ENDDO
!        NEDGETOT = NEDGETOT/2
!        print *, "edge count = ",nedgetot
!
!C  check that adjacency matrix is symmetric
!C
!      SYMMETRIC = .true.
!      DO INODE = 1, MNP
!      DO J = 1, NEDGES(INODE)
!         JNODE = CO_NODES(J,INODE)
!         FOUND = .false.
!         DO K= 1, NEDGES(JNODE
!            IF (CO_NODES(K,JNODE) == INODE) THEN
!              FOUND = .true.
!              EXIT
!            ENDIF
!         ENDDO
!         IF (.not. FOUND) THEN
!           SYMMETRIC = .false.
!      print *, "node ",inode," adjacent to ",jnode," but not visa-versa"
!         ENDIF
!      ENDDO
!      ENDDO
!      IF (.not. SYMMETRIC) THEN
!         stop 'bad adjacency matrix: not symmetric!'
!      ENDIF
ENDIF
      RETURN
      END SUBROUTINE DECOMPOSE  
END MODULE PARPREP  
