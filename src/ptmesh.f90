MODULE PTMESH_MODULE

CONTAINS

  SUBROUTINE PTMESH(NELEM, INTMA, NPOIN, MPI_SIZE_P, NSIDE, ISIDE,&
       &                    NELEM_PP, NPOIN_PP, NBOUN_PP, NSIDE_PP, NEGRP,&
       &                    LCOMM_PP, GCOMM, ELGRP,&
       &                    NGEOM, NBNOI, NBOUN, NBNOR, NNODE, NX, NY, EL, GEOME, MMAT,&
       &                    CMMAT, COORD, BSIDO, RSIDO, maxNSIDE_PP, maxNBOUN_PP,&
       &                    maxNELEM_PP, maxNPOIN_PP, NX_PP, NY_PP, EL_PP,&
       &                    GEOME_PP, MMAT_PP, CMMAT_PP, INTMA_PP, COORD_PP, BSIDO_PP,&
       &                    RSIDO_PP, IMMAT, MPI_RANK_P, IPCOM_PP, IBCOM_PP, ISCOM_PP, NPGRP,&
       &                    ISIDE_PP, NCOMM_PP, SDCOM_PP, IVD, MPI_COMM_P)

    USE INPUT_VARIABLE_MODULE, ONLY: INPUTVARIABLES
    USE GETLOC_MODULE, ONLY: GETLOC
    IMPLICIT NONE
    INCLUDE 'mpif.h'

! *** THIS SUBROUTINE TAKES THE P-SPACE MESH AND NODAL PARTITION GENERATED IN
! *** 'METIS' AND CALCULATES THE GLOBAL COMM ARRAY (GCOMM) AND THE LOCAL COMM
! *** ARRAYS (LCOMM_PP)
!
    INTEGER, ALLOCATABLE :: VFLAG(:)
    INTEGER IE, IP, IP1, IP2, IP3, G1, G2, G3, IG, IEG, NCOMM_PP
    INTEGER NELEM, NPOIN, NGRPS, NELEM_PP, NELEM_PP_CP, GRPT, FLAG
    INTEGER NPOIN_PP, NSIDE_PP, NGEOM, NBNOI, NBOUN, NBNOR
    INTEGER MPI_SIZE_P, SIZE1, MPI_RANK_P
    INTEGER NNODE, maxNELEM_PP, maxNBOUN_PP, maxNPOIN_PP
    INTEGER maxNSIDE_PP, NBOUN_PP, IMMAT
    INTEGER NSIDE, GL, GR, GE, IEL, IER, I
    INTEGER ELGRP(NELEM, 2), NODPT(NPOIN), INTMA(3, NELEM)
    INTEGER NEGRP(MPI_SIZE_P), ISIDE(8, NSIDE), BSIDO(NBNOI, NBOUN)
    INTEGER GCOMM(MPI_SIZE_P, MPI_SIZE_P)
    INTEGER LCOMM_PP(maxNSIDE_PP), SDCOM_PP(3, maxNSIDE_PP)
    INTEGER INTMA_PP(NNODE, maxNELEM_PP), ISIDE_PP(8, maxNSIDE_PP)
    INTEGER BSIDO_PP(NBNOI, maxNBOUN_PP)
    INTEGER IPCOM_PP(maxNPOIN_PP), IBCOM_PP(maxNBOUN_PP)
    INTEGER ISCOM_PP(maxNSIDE_PP), NPGRP(MPI_SIZE_P)
    INTEGER MPI_IERR, MPI_STATUS(MPI_STATUS_SIZE)
    INTEGER MPI_COMM_P
!
    REAL NX(NSIDE), NY(NSIDE), EL(NSIDE), GEOME(NGEOM, NELEM)
    REAL MMAT(NNODE, NELEM), CMMAT(3, NNODE, NELEM)
    REAL COORD(2, NPOIN), RSIDO(NBNOR, NBOUN)
    REAL NX_PP(maxNSIDE_PP), NY_PP(maxNSIDE_PP), EL_PP(maxNSIDE_PP)
    REAL COORD_PP(2, maxNPOIN_PP)
    REAL RSIDO_PP(NBNOR, maxNBOUN_PP)
    REAL MMAT_PP(NNODE, maxNELEM_PP)
    REAL CMMAT_PP(3, 3, maxNELEM_PP)
    REAL GEOME_PP(NGEOM, maxNELEM_PP)
!
    TYPE(INPUTVARIABLES) :: IVD
!
    CHARACTER METISPARTITIONNAME*80
!
    ALLOCATE (VFLAG(MPI_SIZE_P))
!
    NGRPS = MPI_SIZE_P
    IF (MPI_RANK_P .EQ. 0) THEN  ! fkjfhkafda
!
! *** READ IN THE METIS NODAL PARTITION FILE (CHANNEL 12)
!
      METISPARTITIONNAME = IVD%PartitionFile
      PRINT *, 'READING METIS PARTITION FILE = ', METISPARTITIONNAME
      OPEN (12, file=METISPARTITIONNAME, status='old')
      WRITE (*, *)
!
! *** INITIALISE ELGRP
!
      CALL IFILLM(ELGRP, NELEM, 2, 0)
!
! *** READ IN THE NODAL PARTITION FROM THE METIS FILE
!
      DO IP = 1, NPOIN
        READ (12, *) NODPT(IP)
      ENDDO
!
! *** CLOSE CHANNEL 12
!
      CLOSE (12)
!
! *** INCREASE THE MPI_RANK_P BY 1 FOR EACH NODE SO THAT NO NODES/ELEMENT ARE
! *** ALLOCATED TO PROCESSOR MPI_RANK_P 0
!
      DO IP = 1, NPOIN
        NODPT(IP) = NODPT(IP) + 1
      ENDDO
!
! *** INITIALISE NEGRP
!
!      CALL IFILLV(NEGRP,NGRPS,0)
      NEGRP = 0
!
! *** LOOP OVER THE ELEMENTS TO DETERMINE THE ELEMENT PARTITION
!

      DO IE = 1, NELEM
        IP1 = INTMA(1, IE)
        IP2 = INTMA(2, IE)
        IP3 = INTMA(3, IE)
        G1 = NODPT(IP1)
        G2 = NODPT(IP2)
        G3 = NODPT(IP3)
        GE = MAX(G1, G2, G3)
        ELGRP(IE, 1) = GE
        NEGRP(GE) = NEGRP(GE) + 1
      ENDDO
    ENDIF     ! END IF MPI_RANK_P.EQ.0  fkjfhkafda
!
! *** LOOP OVER ELEMENTS IN EACH GROUP TO ASSIGN A GROUP ELEMENT NUMBER
!
    CALL MPI_BARRIER(MPI_COMM_P, MPI_IERR)
    DO IG = 1, NGRPS
    IF (MPI_RANK_P .EQ. 0) THEN
      NELEM_PP_CP = NEGRP(IG)
      IF (NELEM_PP_CP .GT. maxNELEM_PP) THEN
        WRITE (*, 205) IG, NELEM_PP, maxNELEM_PP
        WRITE (*, 206)
        STOP
      ENDIF
!
! ***   SEND NELEM_PP TO PROC MPI_RANK_P IG
!
      IF (IG .NE. 1) THEN
        CALL MPI_SEND(NELEM_PP_CP, 1, MPI_INTEGER, IG - 1, 100,&
   &         MPI_COMM_P, MPI_IERR)
      ELSE
        NELEM_PP = NELEM_PP_CP
      ENDIF
    ENDIF   !LINKS WITH IF STATEMENT ON LINE 94
!
! *** RECEIVE NELEM_PP ON PROC MPI_RANK_P IG
!
    IF (((MPI_RANK_P + 1) .EQ. IG) .AND. (IG .NE. 1)) THEN
      CALL MPI_RECV(NELEM_PP, 1, MPI_INTEGER, 0, 100,&
   &       MPI_COMM_P, MPI_STATUS, MPI_IERR)
    ENDIF
!
    IF (MPI_RANK_P .EQ. 0) THEN
      DO IEG = 1, NELEM_PP_CP
! ***     SCROLL THROUGH 1ST COLUMN OF ELGRP TO FIND THE ELEMENTS IN EACH GROUP
        DO IE = 1, NELEM
          GRPT = ELGRP(IE, 1)
          FLAG = ELGRP(IE, 2)
          IF ((GRPT .EQ. IG) .AND. (FLAG .EQ. 0)) THEN
            ELGRP(IE, 2) = IEG
            GOTO 113 ! "EXIT" (the innermost loop, 106)
          ENDIF
        ENDDO
113     CONTINUE
      ENDDO
    ENDIF      !LINKS WITH IF STATMENT ON LINE 118

    CALL MPI_BARRIER(MPI_COMM_P, MPI_IERR)
    ENDDO ! END DO IG=1,NGRPS
    CALL MPI_BARRIER(MPI_COMM_P, MPI_IERR)
!
! *** BROADCAST ELGRP TO ALL PROCESSORS
!
    SIZE1 = 2*NELEM
    CALL MPI_BCAST(ELGRP, SIZE1, MPI_INTEGER, 0, MPI_COMM_P&
   &     , MPI_IERR)
!
! *** INITIALISE GCOMM
!
    IF (MPI_RANK_P .EQ. 0) THEN
      CALL IFILLM(GCOMM, NGRPS, NGRPS, 0)
!
! *** LOOP OVER THE MESH EDGES AND FILL IN GCOMM
!
      DO I = 1, NSIDE
        IEL = ISIDE(3, I)
        GL = ELGRP(IEL, 1)
        IER = ISIDE(4, I)
        IF (IER .EQ. 0) GOTO 111
        GR = ELGRP(IER, 1)
        IF (GL .NE. GR) THEN
          GCOMM(GL, GR) = GCOMM(GL, GR) + 1
          GCOMM(GR, GL) = GCOMM(GR, GL) + 1
        ENDIF
111     CONTINUE
      ENDDO
    ENDIF   !LINKS WITH IF STATEMENT ON LINE 140

    CALL MPI_BCAST(GCOMM, NGRPS*NGRPS, MPI_INTEGER, 0, MPI_COMM_P,&
   &        MPI_IERR)
    CALL MPI_BCAST(NEGRP, NGRPS, MPI_INTEGER, 0, MPI_COMM_P,&
   &        MPI_IERR)
!
    CALL GETLOC(NSIDE=NSIDE, NGEOM=NGEOM, NELEM=NELEM, &
                maxNELEM_PP=maxNELEM_PP, NPOIN=NPOIN, NBNOI=NBNOI, &
                NBOUN=NBOUN, NBNOR=NBNOR, NNODE=NNODE, &
                NGRPS=NGRPS, ELGRP=ELGRP, NEGRP=NEGRP, &
                NX=NX, NY=NY, EL=EL, &
                GEOME=GEOME, MMAT=MMAT, CMMAT=CMMAT, &
                INTMA=INTMA, COORD=COORD, BSIDO=BSIDO, &
                RSIDO=RSIDO, ISIDE=ISIDE, NELEM_PP=NELEM_PP, &
                maxNBOUN_PP=maxNBOUN_PP, NBOUN_PP=NBOUN_PP, maxNPOIN_PP=maxNPOIN_PP, &
                NPOIN_PP=NPOIN_PP, maxNSIDE_PP=maxNSIDE_PP, NSIDE_PP=NSIDE_PP, &
                ISIDE_PP=ISIDE_PP, NX_PP=NX_PP, NY_PP=NY_PP, &
                EL_PP=EL_PP, GEOME_PP=GEOME_PP, MMAT_PP=MMAT_PP, &
                CMMAT_PP=CMMAT_PP, INTMA_PP=INTMA_PP, COORD_PP=COORD_PP, &
                BSIDO_PP=BSIDO_PP, RSIDO_PP=RSIDO_PP, IMMAT=IMMAT, &
                MPI_RANK_P=MPI_RANK_P, LCOMM_PP=LCOMM_PP, IPCOM_PP=IPCOM_PP, &
                IBCOM_PP=IBCOM_PP, ISCOM_PP=ISCOM_PP, NPGRP=NPGRP, &
                NCOMM_PP=NCOMM_PP, SDCOM_PP=SDCOM_PP, MPI_COMM_P=MPI_COMM_P)

    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
    CALL MPI_BCAST(NPGRP, NGRPS, MPI_INTEGER, 0, MPI_COMM_P,&
   &        MPI_IERR)
!
!
! *** FORMAT STATEMENTS
!
205 FORMAT('maxNELEM_PP PARAMETER NOT LARGE ENOUGH IN PTMESH FOR&
  &      PROCESSOR MPI_RANK_P', 3I6)
206 FORMAT('PROGRAM STOPPED')
!
    RETURN
  END

END MODULE PTMESH_MODULE