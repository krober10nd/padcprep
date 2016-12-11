MODULE PRE
USE MESSENGER

IMPLICIT NONE 

!This module contains most of the global variables, read/write subroutines for
!the fort.14 and prepares the data for the call to PARMETIS by building all
!relevant information locally on each MYPROC.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(len=80)           :: dmy !garbage  
CHARACTER(len=*), PARAMETER ::FILEBASE = "MYPROC_"
CHARACTER(LEN=*), PARAMETER ::FILEEND  = ".14"
CHARACTER(LEN=20)           :: FILENAME !filename used to write the subdomain grids
INTEGER                     :: NE, NP,NE_G,NP_G  !number of elements,nodes    
REAL(8),ALLOCATABLE         :: X_G(:),Y_G(:),DP_G(:)
REAL(8),ALLOCATABLE         :: X_LOC(:),Y_LOC(:),DP_LOC(:) 
INTEGER,ALLOCATABLE         :: DMYCOUNT(:)
INTEGER,ALLOCATABLE         :: VTXWGTS_LOC(:),VTXWGTS_LOC2(:),VTXWGTS_G(:) !these are defined for the elements
INTEGER,ALLOCATABLE         :: VSIZES_LOC(:),VSIZES_LOC2(:),VSIZES_G(:) 
INTEGER                     :: CHUNK,CHUNK_NP,LEFTOVER  !num of local ele, num. of local nodes
INTEGER                     :: NOPE,NETA,NBOU,NVEL !boundary information
INTEGER,ALLOCATABLE         :: IEL_LOC(:),NNEL_LOC(:,:),NNEL_G(:,:),EIND(:),EPTR(:),ELMDIST(:) !local ele conn.
INTEGER                     :: IT ! time step counter 

CONTAINS 
!---------------------------------------------------------------------
!      S U B R O U T I N E   R E M O V E   D R Y 
!---------------------------------------------------------------------
!  REMOVES NEVER DRY ELEMENTS AND NODES 
!---------------------------------------------------------------------
SUBROUTINE RM_DRY                  

IMPLICIT NONE 

INTEGER :: I,J,K,O,ITEMP,chunk_temp !counters
INTEGER,ALLOCATABLE :: NNEL_LOC_TRIM(:,:)
REAL(8) :: MINDEPTH,NM1_DP,NM2_DP,NM3_DP
REAL(8) :: AVGDEPTH_LOC(CHUNK) 

MinDepth=0.10D0
AvgDepth_LOC=0D0 
O=1
K=1
DO I = 1,CHUNK 
  NM1_DP=DP_G(EIND(O))
  O=O+1 
  NM2_DP=DP_G(EIND(O)) 
  O=O+1
  NM3_DP=DP_G(EIND(O))
  O=O+1
  AvgDepth_LOC(K)=(NM1_DP+NM2_DP+NM3_DP)/3D0
  K=K+1
ENDDO
K=0
DO I =1,CHUNK !if it's ge than mindepth, keep it
  IF(AVGDEPTH_LOC(I).GE.MinDEPTH) THEN 
    K=k+1
  ENDIF
ENDDO
CHUNK_TEMP=K !RESIZED
ALLOCATE(NNEL_LOC_TRIM(3,CHUNK_TEMP)) 
NNEL_LOC_TRIM= -1
CALL MPI_ALLREDUCE(CHUNK_TEMP,NE_G,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,IERR)
IF(MYPROC.EQ.0) THEN 
  PRINT *, "The resized no. of eles in domain is ", NE_G
ENDIF

K=1
O=1
J=1
DO I = (MYPROC*CHUNK)+1,(MYPROC*CHUNK)+CHUNK !this is global ele number 
  IF(AVGDEPTH_LOC(O).GE.MinDepth) THEN ! keep it
    NNEL_LOC_TRIM(1,K)=NNEL_LOC(1,O)
    NNEL_LOC_TRIM(2,K)=NNEL_LOC(2,O)
    NNEL_LOC_TRIM(3,K)=NNEL_LOC(3,O)
    K=K+1
  ENDIF 
  O=O+1
ENDDO

CHUNK = CHUNK_TEMP
DEALLOCATE(EIND,NNEL_LOC) 
ALLOCATE(NNEL_LOC(3,CHUNK),EIND(3*CHUNK)) 
EIND=-1 
!rebuild eind and return back 
O=1 
DO I = 1,CHUNK
  EIND(O)=NNEL_LOC_TRIM(1,I)
  O=O+1 
  EIND(O)=NNEL_LOC_TRIM(2,I)
  O=O+1  
  EIND(O)=NNEL_LOC_TRIM(3,I)
  O=O+1 
ENDDO
NNEL_LOC=NNEL_LOC_TRIM
RETURN 
END SUBROUTINE RM_DRY

!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT CONNECTIVITY   
!---------------------------------------------------------------------
SUBROUTINE READGRAPH                                      

IMPLICIT NONE 

INTEGER :: I,J,K,O,ITEMP !counters
REAL(8) :: T1,T2

OPEN(13, FILE='fort.14',STATUS='OLD')
!Read title 
READ(13,*) dmy
!Read number of elements and nodes
READ(13,*) NE,NP
NE_G = NE 
NP_G = NP

!naive decomposition is contiguous 
!determine the chunk size or number of eles on each rank
IF(MYPROC.NE.(NPROC-1)) THEN !if not the last PE 
    CHUNK = NE/NPROC
ELSE !if the last PE 
    CHUNK = NE/NPROC
    LEFTOVER = NE - (CHUNK*NPROC)         
    CHUNK = CHUNK + LEFTOVER 
ENDIF
ALLOCATE(X_G(NP_G),Y_G(NP_G),DP_G(NP_G)) 
ALLOCATE(IEL_LOC(CHUNK),NNEL_LOC(3,CHUNK),NNEL_G(3,NE_G),EIND(3*CHUNK),EPTR(CHUNK+1))
!! Read nodal connectivity table
IF(MYPROC.EQ.0) THEN 
  DO I = 1,NP_G 
    READ(13,*) J,X_G(I),Y_G(I),DP_G(I)
  ENDDO
ELSE 
  DO I = 1,NP 
    READ(13,*)
  ENDDO 
ENDIF
CALL MPI_BCAST(X_G,NP_G,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(Y_G,NP_G,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(DP_G,NP_G,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
! Read element connectivity table in chunks
! Here we read in the fort.14 ele. connectivity table by skipping over the parts
! MYPROC isn't getting. 
IF(MYPROC.EQ.0) THEN 
  PRINT *, "READING IN THE ELE TABLE ..."
ENDIF
CALL CPU_TIME(T1)
NNEL_LOC=0
NNEL_G=0  
EIND=0
EPTR=0 
IF(MYPROC.NE.(NPROC-1)) THEN  !if not the last PE
  K=1 
  O=1 
  DO I = 1,NE
    IF(I.GE.(MYPROC*CHUNK+1).AND.I.LE.((MYPROC*CHUNK)+CHUNK)) THEN 
      READ(13,*) IEL_LOC(K),ITEMP,NNEL_LOC(1,K),NNEL_LOC(2,K),NNEL_LOC(3,K)
      EIND(O)=NNEL_LOC(1,K)
      O=O+1 
      EIND(O)=NNEL_LOC(2,K)
      O=O+1  
      EIND(O)=NNEL_LOC(3,K)
      O=O+1 
      K = K + 1 
    ELSE 
      READ(13,*)
    ENDIF
  ENDDO
  ELSE
  K=1
  O=1
   DO I = 1,NE 
     IF(I.GE.MYPROC*(CHUNK-LEFTOVER)+1) THEN 
       READ(13,*) IEL_LOC(K),ITEMP,NNEL_LOC(1,K),NNEL_LOC(2,K),NNEL_LOC(3,K)
       EIND(O)=NNEL_LOC(1,K) 
       O=O+1 
       EIND(O)=NNEL_LOC(2,K) 
       O=O+1 
       EIND(O)=NNEL_LOC(3,K) 
       O=O+1 
       K = K + 1
     ELSE 
       READ(13,*)
     ENDIF
   ENDDO
ENDIF
CLOSE(13)
CALL CPU_TIME(T2)
IF(MYPROC.EQ.0) THEN 
  PRINT *, "FINISHED READING IN THE ELE TABLE IN ",T2-T1
ENDIF

#ifdef TRIM_DRY
IF(MYPROC.EQ.0) THEN
  PRINT *, "Trimming dry elements..."
ENDIF
CALL RM_DRY 
#endif 

K=1
DO I = 1,SIZE(EIND)+1,3 !eles are triangular  
  EPTR(K)=I
  K=K+1 
ENDDO

ALLOCATE(ELMDIST(NPROC+1)) !same as vtxdist since dual graph 
ELMDIST=0 
#ifndef TRIM_DRY 
DO I = 0,NPROC !the ele numbers on each rank
  ELMDIST(I+1)=(I*(CHUNK-LEFTOVER))+1
ENDDO
#else           
DO I = 0,NPROC ! leftover isn't valid anymore. the ele numbers on each rank
  ELMDIST(I+1)=(I*CHUNK)+1
ENDDO
#endif
ELMDIST(NPROC+1)=NE_G+1 

ALLOCATE(VTXWGTS_G(NE_G))
ALLOCATE(VTXWGTS_LOC(CHUNK),VSIZES_LOC(CHUNK)) ! initially vtwgts is size chun
!IF(MYPROC.EQ.0) THEN 
!  OPEN(13,file='VW.txt',STATUS='old') !open file containing vertex weights here 
!  DO I = 1,NE_G 
!    READ(13,*) VTXWGTS_G(I)
!  ENDDO
!  CLOSE(13) 
!ENDIF 
! broadcast the global vertex weights to all ranks
!CALL MPI_BCAST(VTXWGTS_G,NE_G,MPI_INT,0,MPI_COMM_WORLD,IERR)
VTXWGTS_G=1
J=1
DO I = 1,NE_G
  IF(I.GE.ELMDIST(MYPROC+1).and.I.LT.ELMDIST(MYPROC+2)) THEN 
     VTXWGTS_LOC(J)=VTXWGTS_G(I)
     J=J+1
  ENDIF
ENDDO  
VSIZES_LOC = 1 !this doesn't change (yet) so no need to localize YET 
! localize the nodal attributes (X_LOC,Y_LOC,DP_LOC)
! determine the number of unique nodes on each rank 
ALLOCATE(DMYCOUNT(NP_G))
DMYCOUNT=0 
DO I = 1,3*CHUNK 
  DMYCOUNT(EIND(I))=1 
ENDDO
CHUNK_NP=SUM(DMYCOUNT) !number of nodes on each rank 
ALLOCATE(X_LOC(CHUNK_NP),Y_LOC(CHUNK_NP),DP_LOC(CHUNK_NP))
K=1
DO I = 1,NP_G 
  IF(DMYCOUNT(I).NE.0) THEN 
    X_LOC(K)=X_G(I) 
    Y_LOC(K)=Y_G(I) 
   DP_LOC(K)=DP_G(I) 
   K=K+1 
  ENDIF
ENDDO

!#ifdef DEBUG 
!!have each processor write its local grid so we can check for connectivity
!!problems.
!     WRITE(FILENAME,'(A,I3,A)') 'MYPROC_',MYPROC,'.14'
!     OPEN(UNIT=MYPROC+105,FILE=FILENAME)
!!write the title 
!     WRITE(MYPROC+105,85)  dmy
!!write the specs of the grid
!     WRITE(MYPROC+105,100) CHUNK,CHUNK_NP
!     !write nodes
!     DO I = 1,CHUNK_NP 
!       WRITE(MYPROC+105,95) I,X(I),Y(I),DP(I)
!     ENDDO
!     !write ele 
!     DO I = 1,CHUNK 
!       WRITE(MYPROC+105,90) I,3,NNEL_L(1,I),NNEL_L(2,I),NNEL_L(3,I)
!     ENDDO
!     !write boundary information to EOF
!      WRITE(MYPROC+105,80) NOPE! num open bou 
!      WRITE(MYPROC+105,75) NETA ! num of open bou. nodes 
!      WRITE(MYPROC+105,70) NBOU ! num of land bou 
!      WRITE(MYPROC+105,65) NVEL ! num of land bou. nodes    
!!
!100    FORMAT(2I10) !header 
!95     FORMAT(1I10,f10.5,f10.5,f10.5) !nodal connectivity table
!90     FORMAT(5I10) !element connectivity table
!85     FORMAT(A4)  !agrid
!!boundary information 
!80     FORMAT(1I5,"=Number of open boundaries") 
!75     FORMAT(1I5,"=Total number of open boundary nodes.") 
!70     FORMAT(1I5,"=Number of land boundaries.") 
!65     FORMAT(1I5,"=Total number of land boundary nodes.")
!     CLOSE(MYPROC+105)
!#endif DEBUG
RETURN 
END SUBROUTINE READGRAPH
END MODULE PRE
