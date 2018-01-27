        SUBROUTINE PTMESH(NELEM,INTMA,NPOIN,MPI_SIZE_P,NSIDE,ISIDE,& 
     &                    NELEM_PP,NPOIN_PP,NBOUN_PP,NSIDE_PP,NEGRP,& 
     &                    LCOMM_PP,GCOMM,MXCOM_PP,ELGRP,& 
     &                    NGEOM,NBNOI,NBOUN,NBNOR,NNODE,NX,NY,EL,GEOME,MMAT,&
     &                    CMMAT,COORD,BSIDO,RSIDO,maxNSIDE_PP,maxNBOUN_PP,& 
     &                    maxNELEM_PP,maxNPOIN_PP,NX_PP,NY_PP,EL_PP,& 
     &                    GEOME_PP,MMAT_PP,CMMAT_PP,INTMA_PP,COORD_PP,BSIDO_PP,& 
     &                    RSIDO_PP,IMMAT,MPI_RANK_P,IPCOM_PP,IBCOM_PP,ISCOM_PP,NPGRP,& 
     &                    ISIDE_PP,NCOMM_PP,SDCOM_PP,IVD,MPI_COMM_P) 
      IMPLICIT NONE 
      INCLUDE 'mpif.h' 
!
! CREATE STRUCTURE FOR INPUT VARIABLES
!
        TYPE InputVariables
          INTEGER :: TORDER !Quadrature order in Theta for V-SPACE
          INTEGER :: NTIME !Number of timesteps
          INTEGER :: FORCEOUT !output forces at each timestep?
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
! *** THIS SUBROUTINE TAKES THE P-SPACE MESH AND NODAL PARTITION GENERATED IN 
! *** 'METIS' AND CALCULATES THE GLOBAL COMM ARRAY (GCOMM) AND THE LOCAL COMM 
! *** ARRAYS (LCOMM_PP) 
!
      INTEGER, ALLOCATABLE :: VFLAG(:) 
      INTEGER IE,IP,IP1,IP2,IP3,G1,G2,G3,IG,IEG,IS,IB,NCOMM_PP 
      INTEGER NELEM,NPOIN,NGRPS,NELEM_PP,NELEM_PP_CP,GRPT,FLAG 
      INTEGER NPOIN_PP,NSIDE_PP,NGEOM,NBNOI,NBOUN,NBNOR 
      INTEGER MPI_SIZE_P,TAG,SIZE1,TAG1,TAG2,MPI_RANK_P 
      INTEGER NNODE,maxNELEM_PP,maxNBOUN_PP,maxNPOIN_PP 
      INTEGER maxNSIDE_PP,NBOUN_PP,IMMAT 
      INTEGER NSIDE,GL,GR,NT,MXCOM_PP,GE,IEL,IER,I 
      INTEGER ELGRP(NELEM,2),NODPT(NPOIN),INTMA(3,NELEM) 
      INTEGER NEGRP(MPI_SIZE_P),ISIDE(8,NSIDE),BSIDO(NBNOI,NBOUN) 
      INTEGER GCOMM(MPI_SIZE_P,MPI_SIZE_P) 
      INTEGER LCOMM_PP(maxNSIDE_PP),SDCOM_PP(3,maxNSIDE_PP) 
      INTEGER INTMA_PP(NNODE,maxNELEM_PP),ISIDE_PP(8,maxNSIDE_PP) 
      INTEGER BSIDO_PP(NBNOI,maxNBOUN_PP) 
      INTEGER IPCOM_PP(maxNPOIN_PP),IBCOM_PP(maxNBOUN_PP) 
      INTEGER ISCOM_PP(maxNSIDE_PP),NPGRP(MPI_SIZE_P) 
      INTEGER MPI_IERR,MPI_STATUS(MPI_STATUS_SIZE) 
      INTEGER MPI_COMM_P
! 
      REAL NX(NSIDE),NY(NSIDE),EL(NSIDE),GEOME(NGEOM,NELEM) 
      REAL MMAT(NNODE,NELEM),CMMAT(3,NNODE,NELEM) 
      REAL COORD(2,NPOIN),RSIDO(NBNOR,NBOUN) 
      REAL NX_PP(maxNSIDE_PP),NY_PP(maxNSIDE_PP),EL_PP(maxNSIDE_PP) 
      REAL COORD_PP(2,maxNPOIN_PP) 
      REAL RSIDO_PP(NBNOR,maxNBOUN_PP) 
      REAL MMAT_PP(NNODE,maxNELEM_PP) 
      REAL CMMAT_PP(3,3,maxNELEM_PP) 
      REAL GEOME_PP(NGEOM,maxNELEM_PP) 
      INTEGER jason,jason1
!
      TYPE(InputVariables) :: IVD
! 
      CHARACTER METISPARTITIONNAME*80 
!
      ALLOCATE(VFLAG(MPI_SIZE_P))
! 
      NGRPS=MPI_SIZE_P
      IF(MPI_RANK_P.EQ.0)THEN  ! fkjfhkafda 
! 
! *** READ IN THE METIS NODAL PARTITION FILE (CHANNEL 12) 
! 
      METISPARTITIONNAME = IVD%PartitionFile
      PRINT*,'READING METIS PARTITION FILE = ',METISPARTITIONNAME 
      OPEN (12, file=METISPARTITIONNAME, status='old') 
      WRITE(*,*) 
! 
! *** INITIALISE ELGRP 
! 
      CALL IFILLM(ELGRP,NELEM,2,0) 
! 
! *** READ IN THE NODAL PARTITION FROM THE METIS FILE 
! 
      DO 102 IP=1,NPOIN 
      READ(12,*) NODPT(IP) 
 102  CONTINUE 
! 
! *** CLOSE CHANNEL 12 
! 
      CLOSE(12) 
! 
! *** INCREASE THE MPI_RANK_P BY 1 FOR EACH NODE SO THAT NO NODES/ELEMENT ARE 
! *** ALLOCATED TO PROCESSOR MPI_RANK_P 0 
! 
      DO 120 IP=1,NPOIN 
      NODPT(IP)=NODPT(IP)+1 
 120  CONTINUE 
! 
! *** INITIALISE NEGRP 
! 
      CALL IFILLV(NEGRP,NGRPS,0) 
! 
! *** LOOP OVER THE ELEMENTS TO DETERMINE THE ELEMENT PARTITION 
! 
      DO 103 IE=1,NELEM 
        IP1=INTMA(1,IE) 
        IP2=INTMA(2,IE) 
        IP3=INTMA(3,IE) 
        G1=NODPT(IP1) 
        G2=NODPT(IP2) 
        G3=NODPT(IP3) 
        GE=MAX(G1,G2,G3) 
        ELGRP(IE,1)=GE 
        NEGRP(GE)=NEGRP(GE)+1 
 103  CONTINUE 
      ENDIF     ! END IF MPI_RANK_P.EQ.0  fkjfhkafda
! 
! *** LOOP OVER ELEMENTS IN EACH GROUP TO ASSIGN A GROUP ELEMENT NUMBER 
! 
      CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR) 
      DO 104 IG=1,NGRPS 
      IF(MPI_RANK_P.EQ.0)THEN 
      NELEM_PP_CP=NEGRP(IG) 
      IF(NELEM_PP.GT.maxNELEM_PP)THEN 
        WRITE(*,205) IG 
        WRITE(*,206) 
        STOP 
      ENDIF  
! 
! *** SEND NELEM_PP TO PROC MPI_RANK_P IG 
!     
      IF(IG.NE.1)THEN 
        CALL MPI_SEND(NELEM_PP_CP,1,MPI_INTEGER,IG-1,100,& 
     &       MPI_COMM_P,MPI_IERR) 
      ELSE
        NELEM_PP = NELEM_PP_CP
      ENDIF
      ENDIF   !LINKS WITH IF STATEMENT ON LINE 94 
! 
! *** RECEIVE NELEM_PP ON PROC MPI_RANK_P IG 
! 
      IF(((MPI_RANK_P+1).EQ.IG).AND.(IG.NE.1))THEN 
        CALL MPI_RECV(NELEM_PP,1,MPI_INTEGER,0,100,& 
     &       MPI_COMM_P,MPI_STATUS,MPI_IERR) 
      ENDIF 
!	    	             
      IF(MPI_RANK_P.EQ.0)THEN  
        DO 105 IEG=1,NELEM_PP_CP
! ***     SCROLL THROUGH 1ST COLUMN OF ELGRP TO FIND THE ELEMENTS IN EACH GROUP 
          DO 106 IE=1,NELEM 
            GRPT=ELGRP(IE,1) 
            FLAG=ELGRP(IE,2) 
            IF((GRPT.EQ.IG).AND.(FLAG.EQ.0))THEN 
              ELGRP(IE,2)=IEG 
              GOTO 113 ! "EXIT" (the innermost loop, 106)
            ENDIF 
 106      CONTINUE 
 113      CONTINUE 
 105    CONTINUE   
      ENDIF      !LINKS WITH IF STATMENT ON LINE 118 

      CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR) 
 104  CONTINUE   ! END DO 104 IG=1,NGRPS 
       CALL MPI_BARRIER(MPI_COMM_P,MPI_IERR) 
! 
! *** BROADCAST ELGRP TO ALL PROCESSORS 
! 
      SIZE1=2*NELEM 
      CALL MPI_BCAST(ELGRP,SIZE1,MPI_INTEGER,0,MPI_COMM_P& 
     &     ,MPI_IERR) 
! 
! *** INITIALISE GCOMM 
! 
      IF(MPI_RANK_P.EQ.0)THEN 
      CALL IFILLM(GCOMM,NGRPS,NGRPS,0) 
! 
! *** LOOP OVER THE MESH EDGES AND FILL IN GCOMM 
! 
      DO 107 I=1,NSIDE 
      IEL=ISIDE(3,I) 
      GL=ELGRP(IEL,1) 
      IER=ISIDE(4,I) 
      IF(IER.EQ.0) GOTO 111 
      GR=ELGRP(IER,1) 
      IF(GL.NE.GR)THEN 
        GCOMM(GL,GR)=GCOMM(GL,GR)+1 
        GCOMM(GR,GL)=GCOMM(GR,GL)+1 
      ENDIF 
 111  CONTINUE 
 107  CONTINUE 
      ENDIF   !LINKS WITH IF STATEMENT ON LINE 140
      
      CALL MPI_BCAST(GCOMM,NGRPS*NGRPS,MPI_INTEGER,0,MPI_COMM_P,&
     &        MPI_IERR)
      CALL MPI_BCAST(NEGRP,NGRPS,MPI_INTEGER,0,MPI_COMM_P,&
     &        MPI_IERR)
!
      CALL GETLOC(NSIDE,NGEOM,NELEM,maxNELEM_PP,NPOIN,NBNOI,& 
     &      NBOUN,NBNOR,NNODE,MXCOM_PP,NGRPS,ELGRP,NEGRP,GCOMM,NX,NY,& 
     &      EL,GEOME,MMAT,CMMAT,INTMA,COORD,BSIDO,RSIDO,ISIDE,& 
     &      NELEM_PP,maxNBOUN_PP,NBOUN_PP,maxNPOIN_PP,NPOIN_PP,& 
     &      maxNSIDE_PP,NSIDE_PP,ISIDE_PP,& 
     &      NX_PP,NY_PP,EL_PP,GEOME_PP,MMAT_PP,CMMAT_PP,& 
     &      INTMA_PP,COORD_PP,BSIDO_PP,RSIDO_PP,IMMAT,& 
     &      MPI_RANK_P,LCOMM_PP,IPCOM_PP,IBCOM_PP,ISCOM_PP,NPGRP,&
     &      NCOMM_PP,SDCOM_PP,MPI_COMM_P)
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
      CALL MPI_BCAST(NPGRP,NGRPS,MPI_INTEGER,0,MPI_COMM_P,&
     &        MPI_IERR)
!
! 
! *** FORMAT STATEMENTS 
! 
 201  FORMAT('What is the name of the METIS partition file?  ',$) 
 202  FORMAT('ERROR: NEED TO INCREASE PARAMETER maxmxcom_pp !') 
 203  FORMAT('     in DG_LIN_CONVEC subroutine   ') 
 204  FORMAT('   PROGRAM STOPPED!     ') 
 205  FORMAT('maxNELEM_PP PARAMETER NOT LARGE ENOUGH IN PTMESH FOR& 
     &      PROCESSOR MPI_RANK_P',I5) 
 206  FORMAT('PROGRAM STOPPED') 
!
      RETURN
      END 
