MODULE PARPREP
! interface for parmetis adaptive repartition and kway
USE PRE,ONLY:NNEL_LOC,NNEL_G,VTXWGTS_G,VTXWGTS_LOC,VSIZES_LOC,IT,ELMDIST,EPTR,EIND,CHUNK,CHUNK_NP,NP_G,NE_G
USE MESSENGER
USE ISO_C_BINDING 

IMPLICIT NONE

INTEGER,POINTER              :: XADJ2(:),ADJNCY2(:) ! dual graph structure 
TYPE(C_PTR),ALLOCATABLE      :: PTRXADJ(:),PTRADJNCY(:) !args for parmetis 
INTEGER                      :: OPTS(3),OPTS_RP(4),WGTFLAG=2,NUMFLAG=1,NCON=1,EDGECUT=0,NPARTS ! args for parmetis
INTEGER                      :: NCOMMONNODES=2 !args for parmetis 
INTEGER,ALLOCATABLE          :: ELMWGT(:),ADJWGT_LOC(:) ! args for parmetis
REAL(8)                      :: UBVEC=1.05D0,ITR=1000D0 !args for parmetis 
REAL(8),ALLOCATABLE          :: TPWGTS(:,:) !args for parmetis
INTEGER,ALLOCATABLE          :: PARTE_LOC(:),PARTE_LOC_OLD(:),PARTE_G(:)!partition labels for elements/dual graph nodes
INTEGER,ALLOCATABLE          :: PARTN_LOC(:),PARTN_G(:)!partition labels for nodes/nodal graph
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

!THIS DOESN'T HAVE TO EQUAL 
NPARTS = NPROC         
IF(IT.EQ.1) THEN
  ALLOCATE(PARTE_G(NE_G),PARTE_LOC(CHUNK),PARTE_LOC_OLD(CHUNK),PARTN_G(NP_G))
  PARTE_G=-1 !this are based on FEM nodes 
  PARTE_LOC=-1 !this is for the dual graph nodes 
  PARTE_LOC_OLD=-1 !this is so we can run uncoupled NPARTS ~= NPROC
  PARTN_G=-1 !global FEM node partition labels 
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
  ADJWGT_LOC =1 ! anyway, this is ignored because of the wgtflag==2
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
  CALL ASSIGNRANKLABELS 
    
ELSE !IT.NE.1 
  OPTS_RP(1)=1
  OPTS_RP(4)=2 !parmetis uncoupled from subdomains
  CALL CPU_TIME(T1)
  CALL ParMETIS_V3_AdaptiveRepart(ELMDIST,XADJ2,ADJNCY2,VTXWGTS_LOC,VSIZES_LOC,ADJWGT_LOC,WGTFLAG,NUMFLAG,NCON,NPARTS,TPWGTS,UBVEC,ITR,OPTS_RP,EDGECUT,PARTE_LOC_OLD,MPI_COMM_WORLD) 
  CALL CPU_TIME(T2)
  PARTE_LOC=PARTE_LOC_OLD
  IF(MYPROC.EQ.0) THEN 
    PRINT *, "ParMeTiS Adaptive repart cut this many edges ", EDGECUT, " at IT = ",IT, " in ", T2-T1
  ENDIF
  CALL ASSIGNRANKLABELS 
ENDIF 
RETURN
END SUBROUTINE DECOMPOSE_PAR 

!---------------------------------------------------------------------
!      S U B R O U T I N E   A S S I G N    R A N K   L A B E L S
!---------------------------------------------------------------------
!  develops nodal and elemental partition labels for local and global
!---------------------------------------------------------------------
SUBROUTINE ASSIGNRANKLABELS
IMPLICIT NONE 

INTEGER :: I,J,K,O,RECVCOUNTS(NPROC),DISPL(NPROC)
!gather up the element/dual graph node rank labels 
!from the last to parmetis
RECVCOUNTS=0
DISPL=0
RECVCOUNTS(MYPROC+1)=CHUNK 
CALL MPI_ALLREDUCE(MPI_IN_PLACE,RECVCOUNTS,NPROC,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
CALL MPI_GATHERV(PARTE_LOC,CHUNK,MPI_INT,PARTE_G,RECVCOUNTS,DISPL,MPI_INT,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(PARTE_G,NE_G,MPI_INT,0,MPI_COMM_WORLD,IERR) 

! convert dual graph vertices/eles rank labels to rank labels for the FEM nodes 
O=1
DO I = 1,CHUNK
  PARTN_G(EIND(O))=PARTE_LOC(I) !this goes in global nodal position  
  O=O+1 
  PARTN_G(EIND(O))=PARTE_LOC(I)
  O=O+1 
  PARTN_G(EIND(O))=PARTE_LOC(I)
  O=O+1 
ENDDO
! use mpi_max to reduce to a global vector of partition labels for FEM nodes.  
CALL MPI_ALLREDUCE(MPI_IN_PLACE,PARTN_G,NP_G,MPI_INT,MPI_MAX,MPI_COMM_WORLD,IERR)
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
!reset PARTN_G and PARTE_G 
PARTN_G=-1
PARTE_G=-1
RETURN
END SUBROUTINE

END MODULE PARPREP  
