      SUBROUTINE OUTPUT(NNODE ,NPOIN_PP ,NPOIN,NELEM ,NELEM_PP,COORD,& 
     &               INTMA , DISNF_PP, VNPNT,RANK,IPCOM_PP,ELGRP,IVD) 
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
	  INTEGER VNPNT,NPOIN,NELEM,NNODE,RANK,IP,IN,NRANK,NELEM_PP 
      INTEGER IPT,J,IE,MPI_IERR,MPI_STATUS(MPI_STATUS_SIZE) 
      INTEGER IE_PP,TAG 
      INTEGER NPOIN_PP,IPCOM_PP(NPOIN_PP),IV,ELGRP(NELEM,2) 
      REAL DISNF_PP(NNODE,VNPNT,NELEM_PP),DISNFPARCEL(NNODE) 
      REAL COORD(2,NPOIN) 
! 
      INTEGER INTMA(NNODE,NELEM) 
!
      TYPE(InputVariables) :: IVD
! *** DECLARE THE CHARACTER VARIABLE FILENAME 
      CHARACTER filename*80 
! 
! *** OUTPUT FOR RESTART. 
! 
! *** SET UP CHANNEL TO WRITE TO RESTART FILE 
! 
	  IF(RANK.EQ.0)THEN 
	  WRITE(*,*) 
      WRITE(*,*) 
      filename = IVD%RestartOutFile
      PRINT*,'OPENING NEW RESTART FILE = ',IVD%RestartOutFile 
      OPEN  (20, file=filename, status='UNKNOWN')  
! 
	  WRITE(20,185)  
      WRITE(20,*)  NPOIN,'    ',VNPNT 
      ENDIF !LINKS WITH IF STATEMENT ON LINE 20 
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
! 
! *** COLLECT THE LOCAL DISNF_PP's AND WRITE TO OUTPUT FILE 
! *** LOOP OVER THE V-SPACE NODES AND GLOBAL ELEMENTS 
!   
	  DO 1000 IV=1,VNPNT 
	  DO 1001 IE=1,NELEM 
      	IF(RANK.EQ.0)THEN 
        	NRANK=ELGRP(IE,1) 
            IE_PP=ELGRP(IE,2) 
        ENDIF 
        CALL MPI_BCAST(NRANK,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR) 
        CALL MPI_BCAST(IE_PP,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR) 
        IF(RANK.EQ.NRANK)THEN 
          DO IN=1,NNODE 
            DISNFPARCEL(IN)=DISNF_PP(IN,IV,IE_PP) 
          ENDDO 
          TAG=1 
          CALL MPI_SEND(DISNFPARCEL,NNODE,MPI_REAL,0,TAG,MPI_COMM_WORLD,& 
     &			MPI_IERR) 
        ENDIF 
        IF(RANK.EQ.0)THEN 
          TAG=1 
          CALL MPI_RECV(DISNFPARCEL,NNODE,MPI_REAL,NRANK,TAG,& 
     &				MPI_COMM_WORLD,MPI_STATUS,MPI_IERR) 
     	  WRITE(20,*) (DISNFPARCEL(IN),IN=1,NNODE) 
        ENDIF 
 1001 CONTINUE 
 1000 CONTINUE	   
	  IF(RANK.EQ.0)THEN 
	  CLOSE(20) 
! 
! *** OUTPUT FOR GID MESHFILE 
! 
! *** SET UP CHANNEL 
! 
      filename = IVD%GIDMeshFile 
      PRINT*,'OPENING GIDMESH FILE = ',IVD%GIDMeshFile 
      OPEN  (19, file=filename, status='unknown') 
      WRITE(*,*) 
! 
	  WRITE(19,210)  
      WRITE(19,300) 
      DO 6000 IP=1,NPOIN 
      WRITE(19,400)IP,(COORD(J,IP),J=1,2) 
 6000 CONTINUE 
 	  WRITE(19,220) 
! 
	  WRITE(19,230) 
      DO 6001 IE=1,NELEM 
      WRITE(19,200)IE,(INTMA(J,IE),J=1,NNODE) 
 6001 CONTINUE 
 	  WRITE(19,240) 
! 
	  CLOSE(19) 
      WRITE(*,*) 
      WRITE(*,*) 
      ENDIF !LINKS WITH IF STATEMENT ON LINE 53 
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
