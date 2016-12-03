MODULE PARPREP
! interface for parmetis adaptive repartition and redistributes nodal
! information
USE PRE,ONLY: VTXWGTS_G,VTXWGTS_LOC,VSIZES_LOC,IT,ELMDIST,EPTR,EIND,CHUNK,CHUNK_NP,NP_G,NE_G,X_G,Y_G,DP_G,X,Y,DP
USE MESSENGER
USE ISO_C_BINDING 
IMPLICIT NONE

         INTEGER,POINTER,SAVE         :: XADJ2(:),ADJNCY2(:) ! dual graph structure 
         INTEGER,SAVE                 :: WGTFLAG=2,NUMFLAG=1,NCON=1,EDGECUT,NPARTS ! args for parmetis
         INTEGER,SAVE                 :: NCOMMONNODES=2 !args for parmetis 
         INTEGER,ALLOCATABLE,SAVE     :: OPTS(3),ELMWGT(:),ADJWGT_LOC(:) ! args for parmetis
         INTEGER,ALLOCATABLE,SAVE     :: G2L(:) ! graph to local node table
         REAL(8),SAVE                 :: UBVEC,ITR=1000D0 !args for parmetis 
         REAL(8),ALLOCATABLE,SAVE     :: TPWGTS(:,:) !args for parmetis
         TYPE(C_PTR),ALLOCATABLE,SAVE :: PTRXADJ(:),PTRADJNCY(:) !args for parmetis 
         INTEGER,ALLOCATABLE,SAVE     :: PARTE_LOC(:),PARTE_G(:),PARTE_GDMY(:)!partition labels for elements/dual graph nodes
         INTEGER,ALLOCATABLE,SAVE     :: PARTN_G(:),PARTN_GDMY(:)!partition labels for nodes/nodal graph
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  P A R 
!---------------------------------------------------------------------
!  ROUTINE TO CALL PARMETIS 4.0
!---------------------------------------------------------------------

 SUBROUTINE DECOMPOSE_PAR
        IMPLICIT NONE 
 
         INTEGER                            :: I,J,O
         REAL(8)                            :: T1,T2 !for timing 
         ! api calls from the static lib
         EXTERNAL ParMETIS_V3_Mesh2Dual
         EXTERNAL ParMETIS_V3_AdaptiveRepart 
         EDGECUT=0
         NPARTS = NPROC         
     IF(IT.EQ.1) THEN
         ALLOCATE(PARTE_G(CHUNK),PARTN_G(NP_G),PARTN_GDMY(NP_G))
         PARTE_G=-1 !this are based on FEM nodes 
         PARTE_GDMY=-1 !this is based on FEM nodes
         PARTE=-1 !this is for the dual graph nodes 
         ! build the dual graph 
         ALLOCATE(PTRXADJ(1))
         ALLOCATE(PTRADJNCY(1))
         CALL ParMETIS_V3_Mesh2Dual(ELMDIST,EPTR,EIND,NUMFLAG,NCOMMONNODES,PTRXADJ,PTRADJNCY,MPI_COMM_WORLD)
         CALL C_F_POINTER(PTRXADJ(1),XADJ2,[CHUNK+1])
         CALL C_F_POINTER(PTRADJNCY(1),ADJNCY2,[XADJ2(CHUNK+1)-1])
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "BUILT DUAL GRAPH" 
         ENDIF
         ALLOCATE(ADJWGT_LOC(XADJ2(CHUNK+1)-1)) !edge weights
         ADJWGT =1 ! anyway, this is ignored because of the wgtflag==2
         UBVEC=1.05D0
         ALLOCATE(TPWGTS(NCON,NPARTS))  
         TPWGTS=1D0/DBLE(NPARTS) !partitions have weights that sum to one
         OPTS(1)=0 !default options used for first IT  
         ! first call to parmetis 
         CALL CPU_TIME(T1)
        CALL ParMETIS_V3_PartMeshKway(ELMDIST,EPTR,EIND,VTXWGTS_LOC,WGTFLAG,NUMFLAG,NCON,NCOMMONNODES,NPARTS,TPWGTS,UBVEC,OPTS,EDGECUT,PART,MPI_COMM_WORLD)
         CALL CPU_TIME(T2)
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "KWay ParMeTiS CUT THIS MANY EDGES ", EDGECUT," in ", T2-T1
         ENDIF

        CALL ASSIGNRANKLABELS 
        ! localize the graph
        CALL LOCALIZE 
 
ELSE !IT.NE.1 
      OPTS(1)=1
      OPTS(3)=2 ! we pass the previous k-way partitioning to adaptive repartition
      CALL CPU_TIME(T1)
      CALL ParMETIS_V3_AdaptiveRepart(ELMDIST,XADJ2,ADJNCY2,VTXWGTS_LOC,VSIZES_LOC,ADJWGT,WGTFLAG,NUMFLAG,NCON,NPARTS,TPWGTS,UBVEC,ITR,OPTS,EDGECUT,PART_G,MPI_COMM_WORLD) 
      CALL CPU_TIME(T2)
      IF(MYPROC.EQ.0) THEN 
        PRINT *, "Adaptive repart ParMeTiS CUT THIS MANY EDGES ", EDGECUT, " at IT = ",IT, " in ", T2-T1
      ENDIF
      
      CALL ASSIGNRANKLABELS 
      CALL LOCALIZE
      
ENDIF 
RETURN
END SUBROUTINE DECOMPOSE_PAR 

!---------------------------------------------------------------------
!      S U B R O U T I N E   F R E E    M E M
!---------------------------------------------------------------------
!  deallocate arrays that store local data so they can be rebuilt
!---------------------------------------------------------------------
SUBROUTINE FREE_MEM
IMPLICIT NONE
IF(IT.NE.1) THEN ! all local data is deallocated  
  DEALLOCATE(G2L,DP_LOC,X_LOC,Y_LOC,VTXWGTS_LOC,VSIZES_LOC,PARTE_LOC,PARTN_LOC,EIND,EPTR) ! this frees up the memory associated with the
ELSE 
  DEALLOCATE(VTXWGTS_LOC,VSIZES_LOC,PART,EIND,EPTR)  
ENDIF
RETURN
END SUBROUTINE

!---------------------------------------------------------------------
!      S U B R O U T I N E   A S S I G N    R A N K   L A B E L S
!---------------------------------------------------------------------
!  develops nodal and elemental partition labels for local and global
!---------------------------------------------------------------------
SUBROUTINE ASSIGNRANKLABELS
IMPLICIT NONE 
INTEGER :: I,J,K,O
! convert dual graph vertices to rank labels for the FEM nodes 
! eind contains all the nodes in each element. 
O=1
DO I = 1,CHUNK
  PARTN_GDMY(EIND(O))=PARTE_LOC(I) 
  PARTN_LOC(O)=PARTE_LOCY
  O=O+1 
  PARTN_GDMY(EIND(O))=PARTE_LOC(I)
  PARTN_LOC(O)=PART(I)
  O=O+1 
  PARTN_GDMY(EIND(O))=PARTE_LOC(I)
  PARTN_LOC(O)=PART(I)
  O=O+1 
ENDDO
! use mpi_max to reduce to a global vector of partition labels by taking the max 
CALL MPI_ALLREDUCE(PARTN_GDMY,PARTN_G,NP_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
!cat together all the partE_LOC arrays to form a vector of dual graph node/ele
!partition labels then bcast it to all ranks. 
ALLOCATE(PARTE_ROOT(NE_G))
CALL MPI_GATHER(PARTE_LOC,SIZE(PARTE_LOC),MPI_INT,PARTE_ROOT,SIZE(PARTE_LOC),MPI_INT,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST() 
IF(EDGECUT.GT.0) THEN 
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "writing nodal graph partition to file: partmesh.txt"
    OPEN(990,FILE='partmesh.txt')
    DO I=1, NP_G
      WRITE(990,*) PART_G(I)
    ENDDO
    CLOSE(990)
  ENDIF
ENDIF
RETURN
END SUBROUTINE

!---------------------------------------------------------------------
!      S U B R O U T I N E   L O C A L I Z E
!---------------------------------------------------------------------
!  in order to use partmetis_adaptiveRepart, one must redistribute the graph
!  according to the partition labels and then call mesh2dual to rebuild the
!  dual graph again-- corresponding with the redistributed graph. 
!---------------------------------------------------------------------
SUBROUTINE LOCALIZE
IMPLICIT NONE 
INTEGER :: O,I,J,K 
INTEGER :: NO_LOC_NODE,NO_LOC_ELE
EXTERNAL ParMETIS_V3_Mesh2Dual
!local variables/arrays must be deallocated at the end of the timestep
CALL FREE_MEM
!determine the number of dual graph nodes/eles on each rank
K = 1
DO I = 1,SIZE(PARTE_LOC)
  IF(MYPROC+1.EQ.PART(I)) THEN
    K = K+1
  ENDIF
ENDDO
NO_LOC_ELE = K 
!reallocate eind and eptr 
ALLOCATE(EIND(NO_LOC_ELE),EPTR(3*NO_LOC_ELE))
!rebuild EIND and EPTR corresponding with redistributed graph 
EIND=-1
EPTR=-1 
K=1 
O=1 
DO I = 1,SIZE(PARTE_G)
  IF(MYPROC+1.EQ.PARTE_G(I)) THEN
    EIND(O)=NNEL(1,I)
    O=O+1 
    EIND(O)=NNEL(2,I)
    O=O+1 
    EIND(O)=NNEL(3,I)
    O=O+1 
  ENDIF
ENDDO
K=1
DO I = 1,SIZE(EIND)+1,3 !eles are triangular  
  EPTR(K)=I
  K=K+1 
ENDDO

!rebuild elmdist
ALLOCATE(ELMDIST(NPROC+1))  
! THIS NEEDS TO BE DONE 


!create dual graph 
CALL ParMETIS_V3_Mesh2Dual(ELMDIST,EPTR,EIND,NUMFLAG,NCOMMONNODES,PTRXADJ,PTRADJNCY,MPI_COMM_WORLD)
CALL C_F_POINTER(PTRXADJ(1),XADJ2,[NO_LOC_ELE+1])
CALL C_F_POINTER(PTRADJNCY(1),ADJNCY2,[XADJ2(NO_LOC_ELE+1)-1])
!determine how many nodes are on each rank now 
!STOPPED HERE DEC 3 3 AM
O=0
DO I = 1,NP_G 
  IF(PART_G(I).EQ.MYPROC+1) THEN 
    O = O + 1 !this is the number of FEM nodes associated with each rank  
  ENDIF
ENDDO
LOC_NODE=O
!localize nodal data, vertex weights and create a global to local
!element map 
 
! and reallocated here with the new size
ALLOCATE(G2L(NP_G),X(LOC_NODE),Y(LOC_NODE),DP(LOC_NODE),VTXWGTS_LOC(LOC_NODE),VSIZES_LOC(LOC_NODE))
! the global X,Y and DP are already in memory.
! identify the nodes that are associated with each rank
    O=1 
    DO I = 1,NP_G
      IF(PART_G(I).EQ.MYPROC+1) THEN !DATA BELONGS TO MYPROC 
        X(O)=X_G(I)
        Y(O)=Y_G(I)
        DP(O)=DP_G(I)
        VTXWGTS_LOC(O)=VTXWGTS_G(PART_G(I)) 
        VSIZES_LOC(O) = 1 !this is constant for now
        G2L(I)=O !global to local node map 
        O=O+1
      ENDIF
    ENDDO 
  RETURN 
END SUBROUTINE

!---------------------------------------------------------------------
!      S U B R O U T I N E   P R E P   G R I D 
!---------------------------------------------------------------------
!  CREATE AN INDUVIDUAL DIRECTORY FOR EACH RANK AND PLACE THE PORTION OF THE
!  GRID IN EACH DIRECTORY
!---------------------------------------------------------------------
   SUBROUTINE PREPGRID
   IMPLICIT NONE 

   INTEGER :: I,PE
   CHARACTER(LEN=10) :: DIRNAME,PENUM

   ! EACH RANK MAKES IT OWN DIR CALLED PE RANK#
         PENUM  = 'PE0000'

  RETURN
   END SUBROUTINE 
END MODULE PARPREP  
