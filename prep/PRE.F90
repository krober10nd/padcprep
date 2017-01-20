MODULE PRE
USE MESSENGER
USE KDTREE2_MODULE  
 
IMPLICIT NONE 

!This module contains most of the global variables, read/write subroutines for
!the fort.14 and prepares the data for the call to PARMETIS by building all
!relevant information locally on each MYPROC.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHARACTER(len=80)           :: dmy !garbage  
!CHARACTER(LEN=20)           :: FILENAME !filename used to write the subdomain grids
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
LOGICAL,SAVE                :: isfirst=.true.
! for kd tree 
TYPE(KDTREE2), POINTER :: TREE
TYPE(KDTREE2_RESULT), ALLOCATABLE :: KDRESULTS(:)


CONTAINS 
!--------------------------------------------------------------------
!      S U B R O U T I N E   B U I L D  KD  T  R  E  E
!--------------------------------------------------------------------
!  BUILDS KD TREE TO STORE NODES FOR QUICK LOOK-UP
!--------------------------------------------------------------------
SUBROUTINE BUILDKDTREE

IMPLICIT NONE 

REAL(8),ALLOCATABLE :: bcxy(:,:) 
allocate(bcxy(2,NP_G))

bcxy(1,:)=X_G
bcxy(2,:)=Y_G 

! build kdtree 
tree => kdtree2_create(bcxy,rearrange=.true.,sort=.true.)
ALLOCATE(KDRESULTS(100)) !get 100 neighbors

return 
end subroutine 

!---------------------------------------------------------------------
!      S U B R O U T I N E   R E M O V E   D R Y 
!---------------------------------------------------------------------
!  REMOVES DRY ELEMENTS ABOVE MIN DEPTH. Additionally, the element must feature all their neighbors
!  also dry above some min-depth, THEN RE-WRITES
!  THE GRAPH TO A FILE CALLED 'FORT.14.trim'
!---------------------------------------------------------------------
SUBROUTINE RM_DRY                  

IMPLICIT NONE 

INTEGER :: I,J,JJ,K,O !counters
INTEGER :: print_no,pe,numNeighs
INTEGER,ALLOCATABLE :: eleNeigh(:) 
INTEGER,ALLOCATABLE :: NNEL_LOC_TRIM(:,:),EIND_TRIM(:)
REAL(8) :: MINDEPTH,NM1_DP,NM2_DP,NM3_DP,NM1_X,NM2_X,NM3_X,NM1_Y,NM2_Y,NM3_Y
REAL(8) :: AVGDEPTH_LOC(chunk),AVGX_LOC(chunk),AVGY_LOC(chunk) 
LOGICAL :: finish=.false.
INTEGER                     :: KOUNT

if(myproc.eq.0) then 
  PRINT *, "TRIMMING GRAPH..."
  print *, "Enter min. depth to consider initially dry"
  read *, MinDepth
endif

call mpi_bcast(MinDepth,1,mpi_double,0,mpi_comm_world,ierr)
AvgDepth_LOC=0D0 

if(myproc.eq.0) then 
   print *, "building kd tree..." 
endif
call buildkdtree
if(myproc.eq.0) then 
  print *, "built kd tree!"
endif

!calculate avgerage elemental depth below geoid
O=1
DO I = 1,CHUNK
  NM1_DP=DP_G(EIND(O))
  NM1_X=X_G(EIND(O))
  NM1_Y=Y_G(EIND(O))
  O=O+1 
  NM2_DP=DP_G(EIND(O)) 
  NM2_X=X_G(EIND(O)) 
  NM2_Y=Y_G(EIND(O)) 
  O=O+1
  NM3_DP=DP_G(EIND(O))
  NM3_X=X_G(EIND(O))
  NM3_Y=Y_G(EIND(O))
  O=O+1
  AvgDepth_LOC(I)=(NM1_DP+NM2_DP+NM3_DP)/3D0
  AvgX_LOC(I)=(NM1_X+NM2_X+NM3_X)/3D0
  AvgY_LOC(I)=(NM1_Y+NM2_Y+NM3_Y)/3D0
ENDDO

kount=1
print_no=0
DO I =1,chunk !if it's gt than mindepth, keep it
  IF(AVGDEPTH_LOC(I).gt.MinDEPTH) THEN
       ! search the nearest 100 points 
        call kdtree2_n_nearest(tp=tree,qv=(/AvgX_LOC(I),AvgY_LOC(I)/), nn=100,results=KDRESULTS)
!test the k-d tree to make sure its working
#ifdef DEBUG
        if(myproc.eq.0.and.kount.eq.1) then 
          print *, "for element ",I," located at ",AvgX_LOC(I),AvgY_LOC(I)   
          open(unit=13,file='nneigh.debug') 
          write(13,5)"nearest neighbors"  
          do K = 1,100
            write(13,66) X_G(KDRESULTS(K)%idx),Y_G(KDRESULTS(K)%idx)
          enddo
          close(13)
          kount=kount+1
        endif
#endif
        print_no=print_no+1
  ENDIF
ENDDO

IF(print_no.ne.0) then 
  ALLOCATE(NNEL_LOC_TRIM(3,PRINT_NO))
endif

K=1 
jj=1
DO I = 1,chunk
  IF(AVGDEPTH_LOC(I).GT.MinDepth) THEN
       NNEL_LOC_TRIM(1,jj)=NNEL_LOC(1,K)
       NNEL_LOC_TRIM(2,jj)=NNEL_LOC(2,K)
       NNEL_LOC_TRIM(3,jj)=NNEL_LOC(3,K)
       jj = jj + 1
  ENDIF
  K = K + 1 
ENDDO

CALL mpi_allreduce(jj-1,ne_g,1,mpi_int,mpi_sum,mpi_comm_world,ierr)
if(myproc.eq.0) then 
  print *, "resized the problem to ",ne_g," elements."
endif

if(myproc.eq.0) then 
  open(unit=13,file='fort.14.trim') 
  write(13,5)"trimmed graph"  
  write(13,6) ne_g,np_g 
  do i = 1,np_g 
    write(13,55) i,x_g(i),y_g(i),dp_g(i)
  enddo
close(13)
endif 
! wait until rank 0 finishes writing nodal connectivity
call mpi_barrier(mpi_comm_world,ierr) 

pe = -1 
jj=1
do while(finish.eqv..false.) 
  pe = pe + 1 
  if(myproc.eq.pe) then 
     !print *, MYPROC ," has opened the file. Printing ", print_/no 
     open(unit=13,file='fort.14.trim',access='append')
     do i = 1,print_no
       write(13,77) jj,3,nnel_loc_trim(1,i),nnel_loc_trim(2,i),nnel_loc_trim(3,i) !write nnel 
       jj=jj+1 
     enddo
     close(13) 
   endif

  call mpi_bcast(jj,1,mpi_int,pe,mpi_comm_world,ierr)
  ! wait until all ranks finish writing ele connectivity
  call mpi_barrier(mpi_comm_world,ierr) 
  if(pe.eq.nproc-1) then 
    finish=.true. 
  endif 
enddo 

! this will be rebuilt in read_graph 
deallocate(IEL_LOC,NNEL_LOC,NNEL_G,EIND,EPTR)

5  format(a25)
6  format(2I15)
66 format(2f25.16)
55 format(I10,3f25.16)
77 format(5I9) 


RETURN 
END SUBROUTINE RM_DRY

!---------------------------------------------------------------------
!      S U B R O U T I N E   R E A D _ G R A P H 
!---------------------------------------------------------------------
!  READS IN NODAL, ELEMENT CONNECTIVITY   
!---------------------------------------------------------------------
SUBROUTINE READGRAPH                                      

IMPLICIT NONE 

INTEGER :: I,J,K,O,ITEMP
REAL(8) :: T1,T2
CHARACTER(LEN=256) :: FILENAME
 
filename='fort.14'
#ifdef TRIM_DRY 
9876 continue  
if(isfirst.eqv..false.) then
  filename='fort.14.trim'
endif
#endif

OPEN(13, FILE=FILENAME,STATUS='OLD')
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

!important allocate statement
ALLOCATE(IEL_LOC(CHUNK),NNEL_LOC(3,CHUNK),NNEL_G(3,NE_G),EIND(3*CHUNK),EPTR(CHUNK+1))

if(isfirst.eqv..true.) then 
  allocate(X_G(NP_G),Y_G(NP_G),DP_G(NP_G)) 
  ! Read nodal connectivity table
  IF(MYPROC.EQ.0) THEN 
    DO I = 1,NP_G 
      READ(13,*) J,X_G(I),Y_G(I),DP_G(I)
    ENDDO
  ELSE 
    DO I = 1,NP_g 
      READ(13,*)
    ENDDO 
  ENDIF
  CALL MPI_BCAST(X_G,NP_G,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST(Y_G,NP_G,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST(DP_G,NP_G,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
else ! don't re-read nodal table in since it's already in memory 
  do i=1,np_g
    read(13,*)
  enddo
endif

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
if(isfirst.eqv..true.) then 
  CALL RM_DRY ! this rewrites the graph to a new file and restarts to top of PRE
  isfirst=.FALSE.
  GOTO 9876 
endif
#endif

K=1
DO I = 1,SIZE(EIND)+1,3 !eles are triangular  
  EPTR(K)=I
  K=K+1 
ENDDO

ALLOCATE(ELMDIST(NPROC+1)) !same as vtxdist since dual graph 
ELMDIST=0 
DO I = 0,NPROC !the ele numbers on each rank
  ELMDIST(I+1)=(I*(CHUNK-LEFTOVER))+1
ENDDO
ELMDIST(NPROC+1)=NE_G+1
! initial localization of dual graph vertex weights 
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
