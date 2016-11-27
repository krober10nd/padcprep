MODULE PARPREP
USE PRE,ONLY: ELMDIST,EPTR,EIND,CHUNK,NP_G
IMPLICIT NONE
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  P A R 
!---------------------------------------------------------------------
!  ROUTINE TO CALL PARMETIS 4.0
!---------------------------------------------------------------------

      SUBROUTINE DECOMPOSE_PAR
      USE MESSENGER 
      USE ISO_C_BINDING 
         IMPLICIT NONE 
         INTEGER             :: WGTFLAG,NUMFLAG,NCON,EDGECUT,NPARTS 
         INTEGER             :: I,O
         INTEGER             :: NCOMMONNODES 
         INTEGER,ALLOCATABLE :: OPTS(:),EPART(:),PART_G(:),PART_GLoc(:),ELMWGT(:) 
         REAL(8)             :: UBVEC,T1,T2
         REAL(8),ALLOCATABLE :: TPWGTS(:,:)
! for building the graph in parallel 
         TYPE(C_PTR),ALLOCATABLE :: PTRXADJ(:),PTRADJNCY(:)
         INTEGER,POINTER :: XADJ2(:),ADJNCY2(:)
! api calls from the static lib
         EXTERNAL ParMETIS_V3_PartMeshKway
         EXTERNAL ParMETIS_V3_Mesh2Dual
         EDGECUT=0 
         ALLOCATE(EPART(CHUNK)) !partition labels for locally stored elements 
         ALLOCATE(PART_G(NP_G)) ! partition labels global node numbers
         ALLOCATE(PART_GLoc(NP_G)) !used in the all reduce op
         EPART=-1
         PART_G=-1
         PART_GLoc=-1
         ALLOCATE(ELMWGT(CHUNK)) 
         ALLOCATE(OPTS(3)) 
         !!! opts to pass to parmetis
         ELMWGT=1 ! this is ignored since wgtflag 
         WGTFLAG=0 !use edge weights? 
         NUMFLAG=1 !fortran style 
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "ENTER NUMBER OF PARTITIONS " 
           READ (*,*) NPARTS
         ENDIF
         !NPARTS=480
         CALL MPI_BCAST(NPARTS,1,MPI_INT,0,MPI_COMM_WORLD,IERR)
         NCON=1 !number of vertex constraints
         UBVEC=1.05D0
         ALLOCATE(TPWGTS(NCON,NPARTS))  
         TPWGTS=1D0/DBLE(NPARTS) !partitions have weights that sum to one
         NCOMMONNODES=2 !triangular mesh 
         OPTS(1)=0 !default options used 
         CALL CPU_TIME(T1)
         CALL ParMETIS_V3_PartMeshKway(ELMDIST,EPTR,EIND,ELMWGT,WGTFLAG,NUMFLAG,NCON,NCOMMONNODES,NPARTS,TPWGTS,UBVEC,OPTS,EDGECUT,EPART,MPI_COMM_WORLD)
         CALL CPU_TIME(T2)
         IF(MYPROC.EQ.0) THEN 
           PRINT *, "ParMeTiS CUT THIS MANY EDGES ", EDGECUT 
           PRINT *, "In ", T2-T1
         ENDIF

      ! label vertices associated with elements 
       O=1
      DO I = 1,CHUNK
          PART_GLoc(EIND(O))=EPART(I) 
          O=O+1 
          PART_GLoc(EIND(O))=EPART(I)
          O=O+1 
          PART_GLoc(EIND(O))=EPART(I)
          O=O+1 
      ENDDO
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
     ! build the graph 
     ALLOCATE(PTRXADJ(1))
     ALLOCATE(PTRADJNCY(1))
     CALL ParMETIS_V3_Mesh2Dual(ELMDIST,EPTR,EIND,NUMFLAG,NCOMMONNODES,PTRXADJ,PTRADJNCY,MPI_COMM_WORLD)
     CALL C_F_POINTER(PTRXADJ(1),XADJ2,[CHUNK+1])
     CALL C_F_POINTER(PTRADJNCY(1),ADJNCY2,[XADJ2(CHUNK+1)-1])
     IF(MYPROC.EQ.0) THEN 
       PRINT *, "BUILT DUAL GRAPH" 
     ENDIF
  RETURN 
   END SUBROUTINE DECOMPOSE_PAR 
END MODULE PARPREP  
