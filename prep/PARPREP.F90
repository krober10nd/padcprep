MODULE PARPREP
IMPLICIT NONE 
CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   D E C O M P O S E  P A R 
!---------------------------------------------------------------------
!  ROUTINE TO CALL PARMETIS 4.0
!---------------------------------------------------------------------

      SUBROUTINE DECOMPOSE_PAR
      USE PRE,ONLY: ELMDIST,EPTR,EIND,CHUNK,CHUNK_NP
      USE MESSENGER 
         IMPLICIT NONE 
         INTEGER*8             :: WGTFLAG,NUMFLAG,NCON,EDGECUT,NPARTS 
         INTEGER               :: I
         INTEGER*8             :: NCOMMONNODES 
         INTEGER*8,ALLOCATABLE :: OPTS(:),PART(:),ELMWGT(:) 
         REAL(8)             :: UBVEC
         REAL(8),ALLOCATABLE :: TPWGTS(:,:)
         EXTERNAL ParMETIS_V3_PartMeshKway
         EDGECUT=0 
         ALLOCATE(PART(CHUNK_NP))
         ALLOCATE(ELMWGT(CHUNK)) 
         ALLOCATE(OPTS(3)) 
         !PRINT *, "THE ELEMENT DISTRIBUTION IS ", ELMDIST 
         !PRINT *, "THE EPTR IS ", EPTR 
         !PRINT *, "THE EIND IS ", EIND 
         !!! opts to pass to parmetis
         ELMWGT=8 ! this is ingored since 
         WGTFLAG=0 !use edge weights? 
         NUMFLAG=1 !fortran style 
         NPARTS=NPROC  !decompose into this many parts 
         NCON=1 !number of vertex constraints
         UBVEC=1.05D0 !allowable deviation in partition size
         ALLOCATE(TPWGTS(NCON,NPARTS))  
         TPWGTS=1/DBLE(NPARTS) !partitions have weights that sum to one
         NCOMMONNODES=2 !triangular mesh 
         OPTS(1)=0 !default options used 
         !print *,"ELMDIST ", size(elmdist) 
         !print *,"EPTR ", size(eptr) 
         !print *,"EIND ", size(eind) 
         !!!
CALL ParMETIS_V3_PartMeshKway(ELMDIST,EPTR,EIND,ELMWGT,WGTFLAG,NUMFLAG,NCON,NCOMMONNODES,NPARTS,TPWGTS,UBVEC,OPTS,EDGECUT,PART,MPI_COMM_WORLD)
           IF(MYPROC.EQ.0) THEN 
             PRINT *, "ParMETIS CUT THIS MANY EDGES ", EDGECUT 
           ENDIF
          IF(EDGECUT.GT.0) THEN 
           IF(MYPROC.EQ.0) THEN 
             print *, "writing mesh partition to file: partmesh.txt"
             OPEN(990,FILE='partmesh.txt')
             DO I=1, SIZE(PART)
               WRITE(990,*) PART(I)
             ENDDO
           CLOSE(990)
           ENDIF
         ENDIF
         RETURN
      END SUBROUTINE DECOMPOSE_PAR 
END MODULE PARPREP  
