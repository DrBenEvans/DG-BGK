      SUBROUTINE OUTPUT(NNODE,NPOIN_PP,NPOIN,NELEM,NELEM_PP,&
     &              maxNELEM_PP,COORD,INTMA , DISNF_PP, VNPNT,IPCOM_PP,&
     &                     ELGRP,NEGRP,IVD,&
     &                     MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                     MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                     VSPACE_FIRST,VSPACE_LAST) 
! 
		IMPLICIT NONE 
        INCLUDE 'mpif.h' 
!
! CREATE STRUCTURE FOR INPUT VARIABLES
!
        TYPE InputVariables
          INTEGER :: TORDER !Quadrature order in Theta for V-SPACE
          INTEGER :: NTIME !Number of timesteps
          INTEGER :: IMMAT !Lumped mass matrices (no/yes)
          INTEGER :: RS !Using Restart data (no/yes)
          INTEGER :: INF !Is there an inflow (no/yes)
          INTEGER :: NVSPACEPART ! number of VSPACE partitions
          REAL :: CSAFM !Safety factor applied to timestep (Courant)
          REAL :: rv !Radial extent of the V-SPACE
          REAL :: T1 !Initial condition temp
          REAL :: P1 !Initial condition pressure
          REAL :: U0 !Initial condition X-vel
          REAL :: V0 !Initial condition Y-vel
          REAL :: CINF(4)
          REAL :: W
          REAL :: ALPHA !WALL MOLECULAR REFLECTION PARAMETER
          REAL :: R !GAS CONSTANT
          REAL :: d !MOLECULAR DIAMETER
          REAL :: M !MOLAR MASS
          CHARACTER :: LobattoFile*80
          CHARACTER :: PSpaceFile*80
          CHARACTER :: OutFile*80
          CHARACTER :: PartitionFile*80
          CHARACTER :: RestartInFile*80
          CHARACTER :: ResidualFile*80
          CHARACTER :: ResultsFile1*80
          CHARACTER :: ResultsFile2*80
          CHARACTER :: RestartOutFile*80
          CHARACTER :: GIDMeshFile*80
        END TYPE
! 
	  INTEGER VNPNT,NPOIN,NELEM,NNODE,IP,IN,NRANK,NELEM_PP 
      INTEGER maxNELEM_PP  
      INTEGER IPT,J,IE,MPI_IERR,MPI_STATUS(MPI_STATUS_SIZE) 
      INTEGER IE_PP,TAG 
      INTEGER NPOIN_PP,IPCOM_PP(NPOIN_PP),IV,ELGRP(NELEM,2) 
      REAL COORD(2,NPOIN) 
      INTEGER NEGRP(MPI_SIZE_P-1)
! 
      INTEGER INTMA(NNODE,NELEM) 
      ! mpi-stuff for position space partitioning
      INTEGER MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,GROUP_P
      ! mpi-stuff for velocity space partitioning
      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,GROUP_V,IVP
      INTEGER VNPNT_PART,VSPACE_FIRST,VSPACE_LAST,IG,TSIZE
      REAL DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,maxNELEM_PP)
      REAL DISNF(NNODE,VSPACE_FIRST:VSPACE_LAST,NELEM)


!
      TYPE(InputVariables) :: IVD
! *** DECLARE THE CHARACTER VARIABLE FILENAME 
      CHARACTER filename*80 
! 
! *** OUTPUT FOR RESTART. 
! 
! *** SET UP CHANNEL TO WRITE TO RESTART FILE 
! 
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN 
        WRITE(*,*) 
        WRITE(*,*) 
        filename = IVD%RestartOutFile
        PRINT*,'OPENING NEW RESTART FILE = ',IVD%RestartOutFile 
        OPEN  (20, file=filename, status='UNKNOWN')  
    
        WRITE(20,185)  
        WRITE(20,*)  NPOIN,'    ',VNPNT  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
    
      VNPNT_PART = VSPACE_LAST-VSPACE_FIRST+1
      WRITE(*,"(A3,4I5)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &           VSPACE_FIRST,VSPACE_LAST

      DO IG=1,MPI_SIZE_P-1 !
        TSIZE = NEGRP(IG)*3*VNPNT_PART
        IF(MPI_RANK_P.EQ.IG)THEN
          WRITE(*,"(A3,2I5,I10,A12)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &            TSIZE,"BS"
          CALL MPI_SEND(DISNF_PP(:,VSPACE_FIRST:VSPACE_LAST,:),&
     &            TSIZE,MPI_REAL,0,IG,&
     &            MPI_COMM_P,MPI_IERR)
          WRITE(*,"(A3,2I8,A12)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &            "AS"

        ENDIF
        IF(MPI_RANK_P.EQ.0)THEN
          WRITE(*,"(A3,2I5,I10,A12)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &            TSIZE,"BR"
          CALL MPI_RECV(DISNF_PP(:,VSPACE_FIRST:VSPACE_LAST,:),&
     &            TSIZE,MPI_REAL,IG,IG,&
     &            MPI_COMM_P,MPI_IERR)
          WRITE(*,"(A3,2I8,A12)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &            "AR"
          ! Copying DISNF_PP into DISNF    
          DO IE=1,NELEM ! scanning through DISNF
          NRANK = ELGRP(IE,1) ! Checking if we have received the 
                              ! necessary data in this iteration of the
                              ! DO IG,MPI_SIZE_P group
            IF(NRANK.EQ.IG)THEN
              IE_PP=ELGRP(IE,2)
              DISNF(:,:,IE) = DISNF_PP(:,:,IE_PP)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      WRITE(*,"(A3,2I8,A12)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &            "OUTPUT2"

      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
     
      ! only the master ranks in the MPI_COMM_P communicators
      IF(MPI_RANK_P.EQ.0)THEN !jdfcdknsdasdjkha
        TSIZE =  NELEM*3*VNPNT_PART
        DO IVP=0,MPI_SIZE_V-1 !dvkjdkjdjkshasdjkh
          IF((IVP.EQ.MPI_RANK_V).AND.(IVP.NE.0)) THEN
              CALL MPI_SEND(DISNF,TSIZE,MPI_REAL,0,1,&
     &             MPI_COMM_V,MPI_IERR)
          ENDIF          
          IF(0.EQ.MPI_RANK_V) THEN !dsaaanccmna
            IF(IVP.NE.0) CALL MPI_RECV(DISNF,TSIZE,MPI_REAL,IVP,1,&
     &             MPI_COMM_V,MPI_IERR)

            DO IV=VSPACE_FIRST,VSPACE_LAST ! 1 to VNPNT_PART for rank 0
              DO IE=1,NELEM
                WRITE(20,*) (DISNF(IN,IV,IE),IN=1,NNODE)
              ENDDO
            ENDDO
          ENDIF ! IF(0.EQ.MPI_RANK_V) THEN !dsaaanccmna
        ENDDO ! DO IVP=0,MPI_SIZE_V-1 !dvkjdkjdjkshasdjkh

      WRITE(*,"(A3,2I8,A12)")"MPI",MPI_RANK_P,MPI_RANK_V,&
     &            "OUTPUT3"
 
        IF(MPI_RANK_V.EQ.0)THEN !lkjfsdccddna
          CLOSE(20) 
! ***     OUTPUT FOR GID MESHFILE 
! ***     SET UP CHANNEL 
          filename = IVD%GIDMeshFile 
          PRINT*,'OPENING GIDMESH FILE = ',IVD%GIDMeshFile 
          OPEN  (19, file=filename, status='unknown') 
          WRITE(*,*) 
! 
          WRITE(19,210)  
          WRITE(19,300) 
          DO IP=1,NPOIN 
            WRITE(19,400)IP,(COORD(J,IP),J=1,2) 
          ENDDO
          WRITE(19,220) 
! 
          WRITE(19,230) 
          DO IE=1,NELEM 
            WRITE(19,200)IE,(INTMA(J,IE),J=1,NNODE) 
          ENDDO
          WRITE(19,240) 
! 
          CLOSE(19) 
          WRITE(*,*) 
          WRITE(*,*) 
        ENDIF ! IF(MPI_RANK_V.EQ.0)THEN !lkjfsdccddna
      ENDIF ! IF(MPI_RANK_P.EQ.0)THEN !jdfcdknsdasdjkha
! 
! *** format statements 
!	 
  200 FORMAT(10I8) 
  300 FORMAT('   COORDINATES') 
  400 FORMAT(I10,2(2X,F10.5)) 
  180 FORMAT('WHAT WOULD YOU LIKE THE RESTART FILE TO BE CALLED?& 
     & (____________.RES) :    ', $) 
  185 FORMAT('NPOIN     VNPNT') 
  190 FORMAT('WHAT WOULD YOU LIKE THE GID MESHFILE TO BE CALLED?& 
     & (____________.RES) :    ', $) 
  210 FORMAT('mesh dimension = 2 elemtype triangle nnode = 3') 
  220 FORMAT('end coordinates') 
  230 FORMAT('elements') 
  240 FORMAT('end elements') 
! 
      RETURN 
      END 
