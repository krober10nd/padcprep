MODULE PARPREP
! interface for parmetis adaptive repartition and kway
USE PRE,ONLY: VTXWGTS_G,VTXWGTS_LOC,VSIZES_LOC,IT,ELMDIST,EPTR,EIND,CHUNK,CHUNK_NP &
NP_G,NE_G,X_G,Y_G,DP_G,X_LOC,Y_LOC,DP_LOC,L2G,G2L
USE MESSENGER
USE ISO_C_BINDING 

IMPLICIT NONE

INTEGER,POINTER,SAVE         :: XADJ2(:),ADJNCY2(:) ! dual graph structure 
TYPE(C_PTR),ALLOCATABLE,SAVE :: PTRXADJ(:),PTRADJNCY(:) !args for parmetis 
INTEGER,SAVE                 :: WGTFLAG=2,NUMFLAG=1,NCON=1,EDGECUT=0,NPARTS ! args for parmetis
INTEGER,SAVE                 :: NCOMMONNODES=2 !args for parmetis 
INTEGER,ALLOCATABLE,SAVE     :: OPTS(3),ELMWGT(:),ADJWGT_LOC(:) ! args for parmetis
REAL(8),SAVE                 :: UBVEC,ITR=1000D0 !args for parmetis 
REAL(8),ALLOCATABLE,SAVE     :: TPWGTS(:,:) !args for parmetis
INTEGER,ALLOCATABLE,SAVE     :: PARTE_LOC(:),PARTE_LOC_OLD(:),PARTE_G(:),PARTE_GDMY(:)!partition labels for elements/dual graph nodes
INTEGER,ALLOCATABLE,SAVE     :: PARTN_G(:),PARTN_GDMY(:)!partition labels for nodes/nodal graph
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  P A R 
!---------------------------------------------------------------------
!  routine to call parmetis 4.0
!---------------------------------------------------------------------
SUBROUTINE DECOMPOSE_PAR

IMPLICIT NONE 

INTEGER                      :: I,J,O
REAL(8)                      :: T1,T2 !for timing 
! api calls from the static lib
EXTERNAL ParMETIS_V3_Mesh2Dual
EXTERNAL ParMETIS_V3_AdaptiveRepart 
EXTERNAL ParMETIS_V3_PartMeshKWay
EDGECUT=0
!THIS DOESN'T HAVE TO EQUAL 
NPARTS = NPROC         
IF(IT.EQ.1) THEN
  ALLOCATE(PARTE_G(NE_G),PARTE_LOC(CHUNK),PARTE_LOC_OLD(CHUNK),PARTN_G(NP_G),PARTN_GDMY(NP_G))
  PARTE_G=-1 !this are based on FEM nodes 
  PARTE_GDMY=-1 !this is based on FEM nodes
  PARTE_LOC=-1 !this is for the dual graph nodes 
  PARTE_LOC_OLD=-1
  PARTN_G=-1
  PARTN_GDMY=-1
  ! build the dual graph corresponding with naive decomp 
  ALLOCATE(PTRXADJ(1))
  ALLOCATE(PTRADJNCY(1))
  CALL ParMETIS_V3_Mesh2Dual(ELMDIST,EPTR,EIND,NUMFLAG,NCOMMONNODES,PTRXADJ,PTRADJNCY,MPI_COMM_WORLD)
  CALL C_F_POINTER(PTRXADJ(1),XADJ2,[CHUNK+1])
  CALL C_F_POINTER(PTRADJNCY(1),ADJNCY2,[XADJ2(CHUNK+1)-1])
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "ParMeTIS built the dual graph." 
  ENDIF
  ALLOCATE(ADJWGT_LOC(XADJ2(CHUNK+1)-1)) !edge weights
  ADJWGT =1 ! anyway, this is ignored because of the wgtflag==2
  UBVEC=1.05D0
  ALLOCATE(TPWGTS(NCON,NPARTS))  
  TPWGTS=1D0/DBLE(NPARTS) !partitions have weights that sum to one
  OPTS(1)=0 !default options used for first IT  
  ! first call to parmetis 
  CALL CPU_TIME(T1)
  CALL ParMETIS_V3_PartMeshKway(ELMDIST,EPTR,EIND,VTXWGTS_LOC,WGTFLAG,NUMFLAG,NCON,NCOMMONNODES,NPARTS,TPWGTS,UBVEC,OPTS,EDGECUT,PARTE_LOC,MPI_COMM_WORLD)
  CALL CPU_TIME(T2)
  !save partition labels in case nproc!=nparts 
  PARTE_LOC_OLD=PARTE_LOC
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "ParMeTiS KWay cut this many edges ", EDGECUT," in ", T2-T1
  ENDIF
  CALL FREE_MEM
  CALL ASSIGNRANKLABELS 
  CALL LOCALIZE 
    
ELSE !IT.NE.1 
  OPTS(1)=1
  OPTS(3)=2 ! we pass the previous k-way partitioning to adaptive repartition
  CALL CPU_TIME(T1)
  CALL ParMETIS_V3_AdaptiveRepart(ELMDIST,XADJ2,ADJNCY2,VTXWGTS_LOC,VSIZES_LOC,ADJWGT,WGTFLAG,NUMFLAG,NCON,NPARTS,TPWGTS,UBVEC,ITR,OPTS,EDGECUT,PARTE_LOC_OLD,MPI_COMM_WORLD) 
  CALL CPU_TIME(T2)
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "ParMeTiS Adaptive repart cut this many edges ", EDGECUT, " at IT = ",IT, " in ", T2-T1
  ENDIF
  CALL FREE_MEM
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

DEALLOCATE(IEL_LOC,DP_LOC,X_LOC,Y_LOC,VTXWGTS_LOC,VSIZES_LOC,PARTE_LOC,PARTN_LOC,EIND,EPTR,ELMDIST)
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

! map PARTE_LOC from local to global element numbers using the l2g map that was
! created before parmetis was called 
! Note that on the first iteration L2G and G2L are the same so no change in ele
! numbers occurs
K=1
DO I = (MYPROC*CHUNK)+1,(MYPROC*CHUNK)+CHUNK
  L2G(I)=J
  PARTE_LOC(J)=PARTE_LOC(K)
  K=K+1
ENDDO
!cat together all the partE_LOC arrays to form a vector of dual graph node/ele
!partition labels then bcast it to all ranks. 
ALLOCATE(PARTE_ROOT(NE_G))
CALL MPI_GATHER(PARTE_LOC,SIZE(PARTE_LOC),MPI_INT,PARTE_ROOT,SIZE(PARTE_LOC),MPI_INT,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST() 
! convert dual graph vertices to rank labels for the FEM nodes 
! Note that eind contains all the nodes in each element. 
O=1
DO I = 1,SIZE(PARTE_LOC)
  PARTN_GDMY(EIND(O))=PARTE_LOC(I) 
  PARTN_LOC(O)=PARTE_LOC
  O=O+1 
  PARTN_GDMY(EIND(O))=PARTE_LOC(I)
  PARTN_LOC(O)=PART(I)
  O=O+1 
  PARTN_GDMY(EIND(O))=PARTE_LOC(I)
  PARTN_LOC(O)=PART(I)
  O=O+1 
ENDDO
! use mpi_max to reduce to a global vector of partition labels 
CALL MPI_ALLREDUCE(PARTN_GDMY,PARTN_G,NP_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
IF(EDGECUT.GT.0) THEN 
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "writing nodal graph partition to file: partmesh.txt"
    OPEN(990,FILE='partmesh.txt')
    DO I=1, NP_G
      WRITE(990,*) PARTN_G(I)
    ENDDO
    CLOSE(990)
  ENDIF
ENDIF
!update global to local and local to global maps
!determine the new number of dual graph nodes/eles on each rank
K = 1
DO I = 1,CHUNK
  IF(MYPROC+1.EQ.PARTE_LOC(I)) THEN
    IEL_LOC(K)=I !ele number on MYPROC
    K = K+1
  ENDIF
ENDDO
CHUNK = K
G2L=0
L2G= 
DO I = 1,CHUNK 
  G2L(IEL_LOC(I))=I 
ENDDO
CALL MPI_ALLREDUCE(G2L,G2L,NE_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
DO I = 1,CHUNK
  L2G(I)=IEL_LOC(I)
ENDDO
CALL MPI_ALLREDUCE(L2G,L2G,NE_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
RETURN
END SUBROUTINE

!---------------------------------------------------------------------
!      S U B R O U T I N E   L O C A L I Z E
!---------------------------------------------------------------------
!  in order to use partmetis_adaptiveRepart, one must redistribute the graph
!  according to the partition labels and then call mesh2dual to rebuild the
!  dual graph again-- corresponding with the redistributed graph. 
!  however, parmetis also requires that the elmdist be contiguous, so the call to
!  assignranklabels maps the local element numbers to a contiguous always
!  increasing element numbers and then to use the data it maps it back to the
!  global node numbes before writing the partmesh.txt file. 
!---------------------------------------------------------------------
SUBROUTINE LOCALIZE

IMPLICIT NONE 

INTEGER :: O,I,J,K 
EXTERNAL ParMETIS_V3_Mesh2Dual
ALLOCATE(EIND(CHUNK),EPTR(3*CHUNK),ELMDIST(NPROC+1))
!rebuild EIND and EPTR corresponding with redistributed graph 
EIND=-1
EPTR=-1 
ELMDIST=-1
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
ELMDIST(MYPROC+1)=(MYPROC*CHUNK)+1
ELMDIST(NPROC+1)=NE_G+1 
PRINT *, ELMDIST 
!create dual graph 
CALL ParMETIS_V3_Mesh2Dual(ELMDIST,EPTR,EIND,NUMFLAG,NCOMMONNODES,PTRXADJ,PTRADJNCY,MPI_COMM_WORLD)
CALL C_F_POINTER(PTRXADJ(1),XADJ2,[CHUNK+1])
CALL C_F_POINTER(PTRADJNCY(1),ADJNCY2,[XADJ2(CHUNK+1)-1])
!determine how many nodes are on each rank now 
K=0
DO I = 1,NP_G 
  IF(PARTN_G(I).EQ.MYPROC+1) THEN 
    K = K + 1 !this is the number of FEM nodes associated with each rank  
  ENDIF
ENDDO
CHUNK_NP=K
!localize nodal data, vertex weights and create a!element map 
! and reallocated here with the new size
ALLOCATE(X_LOC(CHUNK_NP),Y_LOC(CHUNK_NP),DP_LOC(CHUNK_NP),VTXWGTS_LOC(CHUNK_NP),VSIZES_LOC(CHUNK_NP))
! the global X,Y and DP are already in memory.
! identify the nodes that are associated with each rank
O=1 
DO I = 1,NP_G
  IF(PARTN_G(I).EQ.MYPROC+1) THEN !DATA BELONGS TO MYPROC 
    X_LOC(O)=X_G(I)
    Y_LOC(O)=Y_G(I)
    DP_LOC(O)=DP_G(I)
    VTXWGTS_LOC(O)=VTXWGTS_G(PARTN_G(I)) 
    VSIZES_LOC(O) = 1 !this is constant for now
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
