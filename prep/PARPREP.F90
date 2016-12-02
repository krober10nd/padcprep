MODULE PARPREP
! interface for parmetis adaptive repartition and redistributes nodal
! information
USE PRE,ONLY: VTXWGTS_G,IT,ELMDIST,EPTR,EIND,CHUNK,CHUNK_NP,NP_G,NE_G,X_G,Y_G,DP_G,VTXWGTS_LOC,VSIZES_LOC 
USE MESSENGER
USE ISO_C_BINDING 
IMPLICIT NONE

         INTEGER,POINTER,SAVE:: XADJ2(:),ADJNCY2(:) ! dual graph structure 
         INTEGER,SAVE        :: WGTFLAG,NUMFLAG,NCON,EDGECUT,NPARTS ! args for parmetis
         INTEGER,SAVE        :: NCOMMONNODES !args for parmetis 
         INTEGER,ALLOCATABLE,SAVE :: OPTS(:),EPART(:),ELMWGT(:),ADJWGT(:),VTXDIST(:) ! args for parmetis
         INTEGER,ALLOCATABLE,SAVE :: G2L(:) !localized nodal attribute data 
         REAL(8),ALLOCATABLE,SAVE :: X(:),Y(:),DP(:)
         REAL(8),SAVE             :: UBVEC,ITR !args for parmetis 
         REAL(8),ALLOCATABLE,SAVE :: TPWGTS(:,:) !args for parmetis
         TYPE(C_PTR),ALLOCATABLE,SAVE :: PTRXADJ(:),PTRADJNCY(:) !args for parmetis 
         INTEGER,ALLOCATABLE,SAVE  :: PART(:),PART_G(:),PART_GLoc(:) !partition labels for dualg graph nodes and FEM nodes
         !the size of these arrays never changes since we do not redistribute
         !the graph amongst the ranks
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

     IF(IT.EQ.1) THEN
         ALLOCATE(PART(CHUNK),PART_G(NP_G),PART_GLoc(NP_G))
         PART_G=-1 !this are based on FEM nodes 
         PART_GLoc=-1 !this is based on FEM nodes
         PART=-1 !this is for the dual graph nodes 
         ALLOCATE(OPTS(3)) 
         !!! opts to pass to parmetis
         WGTFLAG=2 !use dual graph vertex weights only 
         NUMFLAG=1 !fortran style 
         
         ! build the dual graph 
         ALLOCATE(PTRXADJ(1))
         ALLOCATE(PTRADJNCY(1))
        CALL ParMETIS_V3_Mesh2Dual(ELMDIST,EPTR,EIND,NUMFLAG,NCOMMONNODES,PTRXADJ,PTRADJNCY,MPI_COMM_WORLD)
        CALL C_F_POINTER(PTRXADJ(1),XADJ2,[CHUNK+1])
        CALL C_F_POINTER(PTRADJNCY(1),ADJNCY2,[XADJ2(CHUNK+1)-1])
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "SUCCESS...BUILT DUAL GRAPH" 
         ENDIF
         ALLOCATE(ADJWGT(XADJ2(CHUNK+1)-1)) !edge weights
         ADJWGT =1 ! anyway, this is ignored because of the wgtflag==2
         ALLOCATE(VTXDIST(NPROC+1)) ! distribution of verticies on each rank
         VTXDIST = ELMDIST ! note that elmdist is identically vtxdist for the dual graph
         NPARTS = NPROC ! this is coupled for parmetis         
         NCON=1 !number of vertex constraints
         UBVEC=1.05D0
         ALLOCATE(TPWGTS(NCON,NPARTS))  
         TPWGTS=1D0/DBLE(NPARTS) !partitions have weights that sum to one
         NCOMMONNODES=2 !triangular mesh 
         OPTS(1)=0 !default options used for first IT  
         ITR = 1000 
         ! first call to parmetis 
         CALL CPU_TIME(T1)
         CALL ParMETIS_V3_PartMeshKway(ELMDIST,EPTR,EIND,VTXWGTS_LOC,WGTFLAG,NUMFLAG,NCON,NCOMMONNODES,NPARTS,TPWGTS,UBVEC,OPTS,EDGECUT,PART,MPI_COMM_WORLD)
         !CALL ParMETIS_V3_AdaptiveRepart(VTXDIST,XADJ2,ADJNCY2,VTXWGTS_LOC,VSIZES_LOC,ADJWGT,WGTFLAG,NUMFLAG,NCON,NPARTS,TPWGTS,UBVEC,ITR,OPTS,EDGECUT,PART,MPI_COMM_WORLD) 
         CALL CPU_TIME(T2)
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "ParMeTiS CUT THIS MANY EDGES ", EDGECUT 
           PRINT *, "In ", T2-T1
         ENDIF
! convert dual graph vertices to rank labels for the FEM nodes 
! eind contains all the nodes in each element. 
! so there are chunk*3 nodes since we have triangular ele
! each row in part corresponds with one element
       O=1
         DO I = 1,CHUNK
          PART_GLoc(EIND(O))=PART(I) 
          O=O+1 
          PART_GLoc(EIND(O))=PART(I)
          O=O+1 
          PART_GLoc(EIND(O))=PART(I)
          O=O+1 
         ENDDO
       ! use mpi_max to reduce to a global vector of partition labels by taking the max 
       CALL MPI_ALLREDUCE(PART_GLoc,PART_G,NP_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
      IF(EDGECUT.GT.0) THEN 
        IF(MYPROC.EQ.0) THEN 
           print *, "writing mesh partition to file: partmesh.txt"
           OPEN(990,FILE='partmesh.txt')
           DO I=1, NP_G
             WRITE(990,*) PART_G(I)
           ENDDO
          CLOSE(990)
        ENDIF
      ENDIF
      ! localize nodal attribute data based on partition labels in PART_G
     CALL LOCALIZE 
 
ELSE !IT.NE.1 
OPTS(1)=1
OPTS(3)=2 ! we must pass the k-way partitioning to adaptive
                               ! repartition
! inherited vertex weights are the only thing that has changed from IT-1 to IT 
CALL CPU_TIME(T1)
CALL ParMETIS_V3_AdaptiveRepart(VTXDIST,XADJ2,ADJNCY2,VTXWGTS_LOC,VSIZES_LOC,ADJWGT,WGTFLAG,NUMFLAG,NCON,NPARTS,TPWGTS,UBVEC,ITR,OPTS,EDGECUT,PART,MPI_COMM_WORLD) 
CALL CPU_TIME(T2)
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "ParMeTiS CUT THIS MANY EDGES ", EDGECUT, " at IT = ",IT 
           PRINT *, "In ", T2-T1
         ENDIF
! convert dual graph vertices to rank labels for the FEM nodes 
       O=1
      DO I = 1,CHUNK
          PART_GLoc(EIND(O))=PART(I) 
          O=O+1 
          PART_GLoc(EIND(O))=PART(I)
          O=O+1 
          PART_GLoc(EIND(O))=PART(I)
          O=O+1 
      ENDDO
!
CALL MPI_ALLREDUCE(PART_GLoc,PART_G,NP_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)

      IF(EDGECUT.GT.0) THEN 
        IF(MYPROC.EQ.0) THEN 
           print *, "writing mesh partition to file: partmesh.txt"
           OPEN(990,FILE='partmesh.txt')
           DO I=1, NP_G
             WRITE(990,*) PART_G(I)
           ENDDO
          CLOSE(990)
        ENDIF
      ENDIF

! free mem localize nodal attribute data 
CALL LOCALIZE

ENDIF 
END SUBROUTINE DECOMPOSE_PAR 

!---------------------------------------------------------------------
!      S U B R O U T I N E   F R E E    M E M
!---------------------------------------------------------------------
!  deallocate arrays that store local data   
!---------------------------------------------------------------------
SUBROUTINE FREE_MEM
USE MESSENGER 
IMPLICIT NONE
IF(IT.NE.1) THEN ! all local data is deallocated  
  DEALLOCATE(G2L,DP,X,Y,VTXWGTS_LOC,VSIZES_LOC) ! this frees up the memory associated with the
ELSE !OTHERWISE JUST DEALLOCATE VTXWGTS_LOC AND VSIZES_LOC SINCE THEY WERE
     !NAIVELY ALLOCATED BASED ON THE WAY THE GRAPH WAS DISTRIBUTED 
  DEALLOCATE(VTXWGTS_LOC,VSIZES_LOC)  
ENDIF
RETURN
END SUBROUTINE

!---------------------------------------------------------------------
!      S U B R O U T I N E   L O C A L I Z E
!---------------------------------------------------------------------
!  localizes nodal attributes   
!---------------------------------------------------------------------
SUBROUTINE LOCALIZE
USE MESSENGER 
IMPLICIT NONE 
INTEGER :: O,I,J 
INTEGER :: LOC_NODE
!these must be deallocated at the end of the timestep
CALL FREE_MEM
!determine how many nodes are on each rank now 
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
   USE MESSENGER 
   IMPLICIT NONE 

   INTEGER :: I,PE
   CHARACTER(LEN=10) :: DIRNAME,PENUM

   ! EACH RANK MAKES IT OWN DIR CALLED PE RANK#
         PENUM  = 'PE0000'

  RETURN
   END SUBROUTINE 
END MODULE PARPREP  
