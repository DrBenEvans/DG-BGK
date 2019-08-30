module ADVNCE_MODULE
contains
  SUBROUTINE ADVNCE(NNODE, NGEOM, NBNOI, &
                    NBOUN, NBNOR, COORD, &
                    BSIDO, NBOUN_PP, NPOIN, &
                    NPOIN_PP, NELEM_PP, NELEM, &
                    NTIME, BSIDO_PP, INTMA_PP, &
                    FLUXP_PP, FLUYP_PP, GEOME_PP, &
                    MMAT_PP, CMMAT_PP, RELEN, &
                    DTE, DELUN_PP, RSIDO_PP, &
                    UMEAN_PP, IMMAT, CSAFM, &
                    UX, UY, NSIDE_PP, &
                    ISIDE_PP, NX_PP, NY_PP, &
                    EL_PP, VNPNT, SUMWEIGHT, &
                    DISNF_PP, rv, VCORD, &
                    RHO_PP, UVEL_PP, VVEL_PP, &
                    PS_PP, TEMP_PP, ALPHA, &
                    CINF, GCOMM, IPCOM_PP, &
                    NCOMM_PP, maxNELEM_PP, RORDER, &
                    TORDER, maxNPOIN_PP, SDCOM_PP, &
                    IVD, FORCEOUT, d, &
                    RGas, M, MPI_RANK_P, &
                    MPI_SIZE_P, MPI_COMM_P, MPI_RANK_V, &
                    MPI_COMM_V, VSPACE_FIRST, &
                    VSPACE_LAST)

    USE INPUT_VARIABLE_MODULE, only: INPUTVARIABLES
    USE ADTHEM_MODULE, ONLY: ADTHEM
    USE CMMINC_MODULE, ONLY: CMMINC
    USE COMPMU_MODULE, ONLY: COMPMU
    USE FLUCON2_MODULE, ONLY: FLUCON2
    USE GETELC_MODULE, ONLY: GETELC
    USE GETFLA_MODULE, ONLY: GETFLA
    USE GETINC_MODULE, ONLY: GETINC
    USE GETMAC_MODULE, ONLY: GETMAC
    USE GETRES_MODULE, ONLY: GETRES
    USE MACROS_MODULE, ONLY: MACROS
    USE EDGFLXOPT_MODULE, ONLY: EDGFLXOPT
    USE INFLOW_MODULE, ONLY: INFLOW
    IMPLICIT NONE
    INCLUDE 'mpif.h'

    INTEGER NNODE, NGEOM, NSIDE_PP, NBNOI, NBNOR, NBOUN_PP
    INTEGER RORDER, TORDER, NCOMM_PP, FORCEOUT
    INTEGER NTIME, ITIME, VNPNT, NPOIN_PP, NELEM_PP
    INTEGER NGRPS, NPOIN, NBOUN
    INTEGER maxNPOIN_PP
    INTEGER GCOMM(MPI_SIZE_P, MPI_SIZE_P), BSIDO(NBNOI, NBOUN)
    INTEGER maxNELEM_PP
    ! For EDGFLX(OPT) communication logic
    INTEGER SENDRECV_TOT_LENGTH_MAX
    INTEGER SENDRECV_MAX_LENGTHS(MPI_SIZE_P)! MAX LENGTHS OF DATA TO SEND/RECV TO/FROM OTHER RANKS
    INTEGER SENDRECV_START_OFFSETS(MPI_SIZE_P) ! OFFSET OF FIRST ELEMENT TO SEND TO ANY GIVEN MPI_RANK_P
    REAL, ALLOCATABLE :: SEND_EDGE_DATA(:)!
    REAL, ALLOCATABLE :: RECV_EDGE_DATA(:)!
    INTEGER, ALLOCATABLE :: SEND_EDGE_DATA_IDX(:)! LOCAL EDGE INDICES ON OPPOSITE MPI_RANK_P
    INTEGER, ALLOCATABLE :: RECV_EDGE_DATA_IDX(:)! LOCAL EDGE INDICES ON CURRENT MPI_RANK_P

    INTEGER MPI_IERR, NELEM, IPCOM_PP(NPOIN_PP)
!
    REAL RSIDO_PP(NBNOR, NBOUN_PP), R, THETA, ETA, ZETA, rv, PI
    REAL GEOME_PP(NGEOM, maxNELEM_PP), COORD(2, NPOIN)
    REAL VCORD(3, VNPNT), CINF(4), NFO_PP(NNODE, NELEM_PP)
    REAL MMAT_PP(NNODE, maxNELEM_PP), CMMAT_PP(3, 3, maxNELEM_PP)
    REAL UX, UY, CSAFM, ALPHA, SPEED, C0, DTE, DTEc, DTEb, DTET
    REAL NX_PP(NSIDE_PP), NY_PP(NSIDE_PP), EL_PP(NSIDE_PP)
    REAL RHO_PP(NPOIN_PP), UVEL_PP(NPOIN_PP), VVEL_PP(NPOIN_PP)
    REAL TEMP_PP(NPOIN_PP), PS_PP(NPOIN_PP), ND_PP(NPOIN_PP)
    REAL RHO(NPOIN), UVEL(NPOIN), VVEL(NPOIN), ND(NPOIN), PS(NPOIN)
    REAL TEMP(NPOIN), RELEN(NELEM), MU_PP(NNODE, NELEM_PP)
    INTEGER VSPACE_FIRST, VSPACE_LAST, PRINT_EVERY
    REAL DISNF_PP(NNODE, VSPACE_FIRST:VSPACE_LAST, maxNELEM_PP)
    REAL RHS_PP(1, 3, NELEM_PP)
    REAL FLUXP_PP(1, 3, NELEM_PP), FLUYP_PP(1, 3, NELEM_PP)
    REAL UMEAN_PP(NELEM_PP), DELUN_PP(1, 3, NELEM_PP), ETA_PP(NBOUN_PP)
    REAL DISND_PP(NNODE, NELEM_PP), DISUX_PP(NNODE, NELEM_PP)
    REAL DISUY_PP(NNODE, NELEM_PP), DISPS_PP(NNODE, NELEM_PP)
    REAL SUMWEIGHT, DTETEST(2), SRITM1
    REAL SUMRES, SUMRESG

    REAL RESIDUAL(VSPACE_FIRST:VSPACE_LAST)
    REAL RESIDUAL_PP(VSPACE_FIRST:VSPACE_LAST)
    REAL LIFT, DRAG, d, RGas, M
!
    INTEGER INTMA_PP(NNODE, NELEM_PP), BSIDO_PP(NBNOI, NBOUN_PP)
    INTEGER ISIDE_PP(8, NSIDE_PP), IMMAT, SDCOM_PP(3, NCOMM_PP)
    ! mpi-stuff for position space partitioning
    INTEGER MPI_RANK_P, MPI_SIZE_P, MPI_COMM_P
    ! mpi-stuff for velocity space partitioning
    INTEGER MPI_RANK_V, MPI_COMM_V

    TYPE(InputVariables) :: IVD
!
    PARAMETER(PI=3.1416)
    DOUBLE PRECISION ITERTIME1
    DOUBLE PRECISION ITERTIME2
    INTEGER :: IE, IG, IP, IV
    !PRINT_EVERY = (VSPACE_LAST+1-VSPACE_FIRST)/5
    PRINT_EVERY = VNPNT + 1  ! No printing
!
! *** OPEN CHANNEL 16 TO WRITE TO THE RESIDUAL FILE
!
    IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN ! write only on master rank

!
      PRINT *, 'OPENING RESIDUAL FILE = ', IVD%ResidualFile
      OPEN (16, file=IVD%ResidualFile, status='UNKNOWN')
      WRITE (*, *)
!
! ***   SET UP CHANNEL 17 & 18 TO WRITE TO RESULTS FILES FOR GiD POST PROCESSING
!
      PRINT *, 'OPENING 1st RESULTS FILE = ', IVD%ResultsFile1
      OPEN (17, file=IVD%ResultsFile1, status='UNKNOWN')
      WRITE (*, *)
      WRITE (*, *)

      PRINT *, 'OPENING 2nd RESULTS FILE = ', IVD%ResultsFile2
      OPEN (18, file=IVD%ResultsFile2, status='UNKNOWN')
      WRITE (*, *)
    ENDIF !dskljdasda
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! *** CODE TO DO THE DISTRIBUTION FUNCTION PRINTOUTS AT VARIOUS MESH POINTS
!

!
    C0 = 0.0

!
! *** GET THE COURANT CONDITION MAX ALLOWABLE TIMESTEP
!
    IF (MPI_RANK_P .EQ. 0) THEN
      CALL ALOTIM(NELEM, CSAFM, RELEN, DTEc, rv)
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! *** CALL MACROS TO CONVERT nf'S INTO MACROSCOPIC QUANTITIES
!
    CALL MACROS(rv=rv, VNPNT=VNPNT, SUMWEIGHT=SUMWEIGHT, &
                NPOIN_PP=NPOIN_PP, NNODE=NNODE, DISNF_PP=DISNF_PP, &
                VCORD=VCORD, INTMA_PP=INTMA_PP, NELEM_PP=NELEM_PP, &
                GEOME_PP=GEOME_PP, NGEOM=NGEOM, ND_PP=ND_PP, &
                RHO_PP=RHO_PP, UVEL_PP=UVEL_PP, VVEL_PP=VVEL_PP, &
                PS_PP=PS_PP, TEMP_PP=TEMP_PP, DISND_PP=DISND_PP, &
                DISUX_PP=DISUX_PP, DISUY_PP=DISUY_PP, DISPS_PP=DISPS_PP, &
                M=M, R=RGas, &
                MPI_COMM_V=MPI_COMM_V, VSPACE_FIRST=VSPACE_FIRST, &
                VSPACE_LAST=VSPACE_LAST)
!
! *** CONVERT MACRO_PPs to MACROS
!
    NGRPS = MPI_SIZE_P

    CALL GETMAC(maxNPOIN_PP=maxNPOIN_PP, NPOIN_PP=NPOIN_PP, IPCOM_PP=IPCOM_PP, &
                ND_PP=ND_PP, RHO_PP=RHO_PP, UVEL_PP=UVEL_PP, &
                VVEL_PP=VVEL_PP, PS_PP=PS_PP, TEMP_PP=TEMP_PP, &
                NPOIN=NPOIN, ND=ND, RHO=RHO, &
                UVEL=UVEL, VVEL=VVEL, PS=PS, &
                TEMP=TEMP, NGRPS=NGRPS, MPI_RANK_P=MPI_RANK_P, &
                MPI_COMM_P=MPI_COMM_P) !
! *** WRITE OUTPUT DATA TO RESULTS FILE
!
    IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN ! write only on master rank
      WRITE (17, 500) 0
      WRITE (18, 500) 0
      DO IP = 1, NPOIN
        WRITE (17, 600) IP, ND(IP), UVEL(IP), VVEL(IP)
        WRITE (18, 600) IP, RHO(IP), PS(IP), TEMP(IP)
      ENDDO
!
! *** LOOP OVER THE PRESCRIBED NR. OF TIMESTEPS
!
      PRINT *, 'NTIME= ', NTIME
      WRITE (*, *)
    ENDIF           !LINKS WITH IF STATEMENT ON LINE 131
!
! *** BROADCAST NTIME TO ALL PROCESSORS
!
    CALL MPI_BCAST(NTIME, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, MPI_IERR)

! *** EDGE FLUX COMMUNICATION LOGIC
    DO IG = 1, NGRPS
      SENDRECV_MAX_LENGTHS(IG) = GCOMM(MPI_RANK_P + 1, IG)
    ENDDO

    SENDRECV_START_OFFSETS(1) = 0
    DO IG = 2, NGRPS
      SENDRECV_START_OFFSETS(IG) = SENDRECV_START_OFFSETS(IG - 1) + &
 &                                   SENDRECV_MAX_LENGTHS(IG - 1)
    ENDDO
    SENDRECV_TOT_LENGTH_MAX = SENDRECV_START_OFFSETS(NGRPS) + &
   &                         SENDRECV_MAX_LENGTHS(NGRPS)

    ALLOCATE (SEND_EDGE_DATA(2*SENDRECV_TOT_LENGTH_MAX))
    ALLOCATE (RECV_EDGE_DATA(2*SENDRECV_TOT_LENGTH_MAX))
    ALLOCATE (SEND_EDGE_DATA_IDX(SENDRECV_TOT_LENGTH_MAX))
    ALLOCATE (RECV_EDGE_DATA_IDX(SENDRECV_TOT_LENGTH_MAX))

    ! DEBUG
    !DO IG=1,NGRPS
    !  CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
    !  IF(MPI_RANK_P.EQ.(IG-1))THEN
    !    WRITE(*,"(3I3,I5)") MPI_RANK_P,MPI_RANK_V,NGRPS, &
    !                              SENDRECV_TOT_LENGTH_MAX
    !    WRITE(*,"(8I5)") SENDRECV_MAX_LENGTHS
    !    WRITE(*,"(8I5)") SENDRECV_START_OFFSETS
    !  ENDIF
    !  CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
    !ENDDO
    !DEBUG

    SEND_EDGE_DATA_IDX = 0
    RECV_EDGE_DATA_IDX = 0
!
! *** BEGIN THE TIMESTEPPING !!!
!

    ITERTIME1 = MPI_WTIME()
    SRITM1 = 0.0
    DO ITIME = 1, NTIME
      IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN
        PRINT *, 'TIMESTEP = ', ITIME ! PRINT TIMESTEP NUMBER TO SCREEN
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)

!
! ***   CALCULATE THE FLUX CONSERVATION PARAMETER (ETA) AT EACH BOUNDARY NODE
!
      CALL FLUCON2(NELEM_PP=NELEM_PP, VNPNT=VNPNT, &
                   DISNF_PP=DISNF_PP, NBNOI=NBNOI, &
                   NBOUN_PP=NBOUN_PP, BSIDO_PP=BSIDO_PP, &
                   NBNOR=NBNOR, RSIDO_PP=RSIDO_PP, &
                   VCORD=VCORD, rv=rv, &
                   ETA=ETA_PP, NSIDE_PP=NSIDE_PP, &
                   ISIDE_PP=ISIDE_PP, SUMWEIGHT=SUMWEIGHT, &
                   R=RGas, MPI_COMM_V=MPI_COMM_V, &
                   VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST, &
                   MPI_RANK_P=MPI_RANK_P)
!
! ***   COMPUTE THE VALUE OF MU FOR CONSTRUCTION OF THE BGK RHS SOURCE TERM
!
      CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
      CALL COMPMU(rv=rv, VNPNT=VNPNT, NNODE=NNODE, &
                  DISNF_PP=DISNF_PP, VCORD=VCORD, NELEM_PP=NELEM_PP, &
                  MU_PP=MU_PP, DISND_PP=DISND_PP, DISUX_PP=DISUX_PP, &
                  DISUY_PP=DISUY_PP, SUMWEIGHT=SUMWEIGHT, d=d, &
                  MPI_COMM_V=MPI_COMM_V, VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST)
!
! ***   COMPUTE THE MAX TIMESTEP BASED ON THE BGK TERM
!
      CALL ALOTIM2(NNODE, NELEM_PP, MU_PP, DTEb)
      CALL MPI_ALLREDUCE(DTEb, DTE, 1, MPI_REAL, MPI_MIN,&
   &                MPI_COMM_P, MPI_IERR)
!
! ***   LOOP OVER THE PROCESSORS TO FIND THE LIMITING BGK TIMESTEP
!

      IF ((MPI_RANK_P .EQ. 0)) THEN
        DTEb = DTE
        DTETEST(1) = DTEb
        DTETEST(2) = DTEc
        DTET = MINVAL(DTETEST)
        CALL MPI_ALLREDUCE(DTET, DTE, 1, MPI_REAL, MPI_MIN,&
   &                  MPI_COMM_V, MPI_IERR)
      ENDIF
!
      CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
      CALL MPI_BCAST(DTE, 1, MPI_REAL, 0, MPI_COMM_P, MPI_IERR)
!
!
! ***   LOOP OVER THE NODES IN VELOCITY SPACE
!
      DO IV = VSPACE_FIRST, VSPACE_LAST
        IF ((MPI_RANK_P .EQ. 0) .AND. ((MOD(IV, PRINT_EVERY)) .EQ. 0)) THEN
          WRITE (*, "(A3,I2,A25,I5)") 'VSR', MPI_RANK_V,&
   &               ": V-space iteration no.", IV
        ENDIF

! ***     SYNCHRONISE
!
        CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! ***     IDENTIFY THE POSITION IN THE VELOCITY SPACE MESH BY UX,UV
!

        ETA = VCORD(1, IV)
        ZETA = VCORD(2, IV)
        R = (rv/2)*(ETA + 1)
        THETA = ZETA*PI
        UX = R*COS(THETA)
        UY = R*SIN(THETA)
        SPEED = SQRT(UX*UX + UY*UY)

! ***     COMPUTE THE VALUE OF THE EQUILIBRIUM DIST FUNC FOR THIS POINT IN V-SPACE
! ***     BASED ON THE BULK CONDITIONS
!
        CALL GTEQNF(NFO_PP, NNODE, NELEM_PP, DISND_PP,&
   &        DISUX_PP, DISUY_PP, DISPS_PP, UX, UY, M)
!
!
! ***     OBTAIN THE FLUXES AT THE POINTS
!
        CALL GETFLA(NELEM_PP=NELEM_PP, NNODE=NNODE, DISNF_PP=DISNF_PP, &
                    FLUXP_PP=FLUXP_PP, FLUYP_PP=FLUYP_PP, UX=UX, &
                    UY=UY, IV=IV, VSPACE_FIRST=VSPACE_FIRST, &
                    VSPACE_LAST=VSPACE_LAST)
!
! ***     ELEMENT CONTRIBUTIONS (STEP 1)
!
        CALL GETELC(NELEM_PP=NELEM_PP, GEOME_PP=GEOME_PP, NGEOM=NGEOM, &
                    NNODE=NNODE, IV=IV, DISNF_PP=DISNF_PP, &
                    FLUXP_PP=FLUXP_PP, FLUYP_PP=FLUYP_PP, UMEAN_PP=UMEAN_PP, &
                    RHS_PP=RHS_PP, DTE=DTE, UX=UX, &
                    UY=UY, MU=MU_PP, NFO_PP=NFO_PP, &
                    VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST)
!
! ***     FILL DELUN(NAMAT=1,NPOIN) WITH THE VALUE SET IN C00
        CALL RFILLA(DELUN_PP, 1, NNODE, NELEM_PP, C0)
        CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)

        CALL EDGFLXOPT(NELEM=NELEM_PP, NSIDE_PP=NSIDE_PP, &
                       ISIDE=ISIDE_PP, RHS=RHS_PP, &
                       NX=NX_PP, NY=NY_PP, &
                       EL=EL_PP, UX=UX, &
                       UY=UY, RSIDO=RSIDO_PP, &
                       BSIDO=BSIDO_PP, NBOUN=NBOUN_PP, &
                       NBNOR=NBNOR, ALPHA=ALPHA, &
                       ETA=ETA_PP, VNPNT=VNPNT, &
                       IV=IV, DISNF=DISNF_PP, &
                       UMEAN=UMEAN_PP, rv=rv, &
                       NGRPS=NGRPS, MPI_RANK_P=MPI_RANK_P, &
                       NCOMM_PP=NCOMM_PP, VCORD=VCORD, &
                       RORDER=RORDER, TORDER=TORDER, &
                       SDCOM_PP=SDCOM_PP, R=RGas, &
                       MPI_COMM_P=MPI_COMM_P, VSPACE_FIRST=VSPACE_FIRST, &
                       VSPACE_LAST=VSPACE_LAST, SENDRECV_TOT_LENGTH_MAX=SENDRECV_TOT_LENGTH_MAX, &
                       SENDRECV_START_OFFSETS=SENDRECV_START_OFFSETS, SEND_EDGE_DATA=SEND_EDGE_DATA, &
                       RECV_EDGE_DATA=RECV_EDGE_DATA, SEND_EDGE_DATA_IDX=SEND_EDGE_DATA_IDX, &
                       RECV_EDGE_DATA_IDX=RECV_EDGE_DATA_IDX)

! ***       MULTIPLY BY DELTA-TIME-POINT
        DO IE = 1, NELEM_PP
          DO IP = 1, NNODE
            RHS_PP(1, IP, IE) = DTE*RHS_PP(1, IP, IE)
          ENDDO
        ENDDO
!
! ***     OBTAIN THE INCREMENTS -PREDICTIONS
!
        IF (IMMAT .EQ. 1) THEN
          CALL GETINC(NNODE=NNODE, MMAT=MMAT_PP, RHS=RHS_PP, &
                      DELUN=DELUN_PP, NELEM_PP=NELEM_PP, RESIDUAL=RESIDUAL_PP, &
                      maxNELEM_PP=maxNELEM_PP, DISNF=DISNF_PP, NFO=NFO_PP, &
                      MU=MU_PP, IV=IV, DTE=DTE, &
                      VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST)
        ELSE
          CALL CMMINC(NNODE=NNODE, CMMAT=CMMAT_PP, RHS=RHS_PP, &
                      DELUN=DELUN_PP, NELEM_PP=NELEM_PP, RESIDUAL=RESIDUAL_PP, &
                      maxNELEM_PP=maxNELEM_PP, DISNF=DISNF_PP, NFO=NFO_PP, &
                      MU=MU_PP, IV=IV, DTE=DTE, &
                      VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST)
        ENDIF
!
! ***     ADD INCREMENTS
!
        CALL ADTHEM(IV=IV, DISNF=DISNF_PP, DELUN=DELUN_PP, &
                    NELEM_PP=NELEM_PP, VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST)
!
! ***     RESET THE INFLOW BOUNDARY VALUES
!
        CALL INFLOW(IV=IV, UX=UX, UY=UY, &
                    NNODE=NNODE, NELEM_PP=NELEM_PP, &
                    INTMA=INTMA_PP, CINF=CINF, &
                    NBOUN=NBOUN_PP, NBNOI=NBNOI, BSIDO=BSIDO_PP, &
                    DISNF=DISNF_PP, &
                    R=RGas, M=M, VSPACE_FIRST=VSPACE_FIRST, &
                    VSPACE_LAST=VSPACE_LAST)

        CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! *** END LOOP OVER VELOCITY SPACE NODES IN THE CURRENT PARTITION
!
      ENDDO

      CALL GETRES(RESIDUAL=RESIDUAL, RESIDUAL_PP=RESIDUAL_PP, VSPACE_FIRST=VSPACE_FIRST, &
                  VSPACE_LAST=VSPACE_LAST, MPI_COMM_P=MPI_COMM_P)

!
! *** CALL MACROS
!
      CALL MACROS(rv=rv, VNPNT=VNPNT, SUMWEIGHT=SUMWEIGHT, &
                  NPOIN_PP=NPOIN_PP, NNODE=NNODE, DISNF_PP=DISNF_PP, &
                  VCORD=VCORD, INTMA_PP=INTMA_PP, NELEM_PP=NELEM_PP, &
                  GEOME_PP=GEOME_PP, NGEOM=NGEOM, ND_PP=ND_PP, &
                  RHO_PP=RHO_PP, UVEL_PP=UVEL_PP, VVEL_PP=VVEL_PP, &
                  PS_PP=PS_PP, TEMP_PP=TEMP_PP, DISND_PP=DISND_PP, &
                  DISUX_PP=DISUX_PP, DISUY_PP=DISUY_PP, DISPS_PP=DISPS_PP, &
                  M=M, R=RGas, &
                  MPI_COMM_V=MPI_COMM_V, VSPACE_FIRST=VSPACE_FIRST, &
                  VSPACE_LAST=VSPACE_LAST)

      IF (FORCEOUT .EQ. 0) THEN !NOT WRITING OUT FORCES TO RESIDUAL FILE
! *** CONSTRUCT GLOBAL MACRO VECTORS
        IF (MOD(ITIME, 1000) .EQ. 0) THEN
          CALL GETMAC(maxNPOIN_PP, NPOIN_PP, IPCOM_PP, ND_PP, RHO_PP,&
   &        UVEL_PP, VVEL_PP, PS_PP, TEMP_PP, NPOIN, ND, RHO, UVEL, VVEL,&
   &        PS, TEMP, NGRPS, MPI_RANK_P, MPI_COMM_P)
! ***       WRITE OUTPUT DATA TO RESULTS FILE
          IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN
! ***         CONSTRUCT GLOBAL MACRO VARIABLE VECTORS
            WRITE (17, 500) ITIME
            WRITE (18, 500) ITIME
            DO IP = 1, NPOIN
              WRITE (17, 600) IP, ND(IP), UVEL(IP), VVEL(IP)
              WRITE (18, 600) IP, RHO(IP), PS(IP), TEMP(IP)
            ENDDO
          ENDIF
        ENDIF
      ELSE  !WRITING OUT FORCES TO RESIDUAL FILE
! ***     CONSTRUCT GLOBAL MACRO VECTORS
        IF (MOD(ITIME, FORCEOUT) .EQ. 0) THEN
          CALL GETMAC(maxNPOIN_PP, NPOIN_PP, IPCOM_PP, ND_PP, RHO_PP,&
   &          UVEL_PP, VVEL_PP, PS_PP, TEMP_PP, NPOIN, ND, RHO, UVEL, VVEL,&
   &          PS, TEMP, NGRPS, MPI_RANK_P, MPI_COMM_P)
          IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN
            WRITE (17, 500) ITIME
            WRITE (18, 500) ITIME
            DO IP = 1, NPOIN
              WRITE (17, 600) IP, ND(IP), UVEL(IP), VVEL(IP)
              WRITE (18, 600) IP, RHO(IP), PS(IP), TEMP(IP)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
!
! *** WRITE TO MAX RESIDUAL FILE
!
      IF (MPI_RANK_P .EQ. 0) THEN
        SUMRES = MAXVAL(RESIDUAL)
        CALL MPI_REDUCE(SUMRES, SUMRESG, 1, MPI_REAL, MPI_MAX,&
   &                             0, MPI_COMM_V, MPI_IERR)
        IF (MPI_RANK_V .EQ. 0) THEN
          SUMRES = SUMRESG
          if (ITIME .EQ. 1) SRITM1 = SUMRES
          IF (IVD%FORCEOUT .EQ. 0) THEN     !IF NOT OUTPUTTING FORCES
            WRITE (16, *) ITIME, LOG(SUMRES/SRITM1)
          ELSE             !IF OUTPUTTING FORCES IN RESIDUAL FILE
            CALL GETFOR(LIFT, DRAG, NBOUN, NPOIN, NBNOI, PS, COORD, BSIDO)
            WRITE (16, *) ITIME, LOG(SUMRES/SRITM1), LIFT, DRAG
          ENDIF
        ENDIF
      ENDIF
!
! ***   SYNCHRONISE
!
      CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! ***   END OF TIMESTEPPING LOOP
!
      ITERTIME2 = MPI_WTIME()
      IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN
        WRITE (*, "(A17,F6.2)") "Iteration time:", (ITERTIME2 - ITERTIME1)
      ENDIF
      ITERTIME1 = ITERTIME2
!
    ENDDO ! DO ITIME=1,NTIME !
    IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN
      WRITE (*, *) "END OF TIME LOOP"
    ENDIF

    DEALLOCATE (SEND_EDGE_DATA)
    DEALLOCATE (RECV_EDGE_DATA)
    DEALLOCATE (SEND_EDGE_DATA_IDX)
    DEALLOCATE (RECV_EDGE_DATA_IDX)

!
! *** CONSTRACT GLOBAL MACRO VECTORS
!
!        CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,UVEL_PP,VVEL_PP,&
!    &         PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,&
!    &           NGRPS,MPI_RANK_P,MPI_COMM_P)
!
! *** WRITE OUTPUT DATA TO RESULTS FILE
!
    IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN
! *** CONSTRUCT GLOBAL MACRO VARIABLE VECTORS
      WRITE (17, 500) ITIME
      WRITE (18, 500) ITIME
      DO IP = 1, NPOIN
        WRITE (17, 600) IP, ND(IP), UVEL(IP), VVEL(IP)
        WRITE (18, 600) IP, RHO(IP), PS(IP), TEMP(IP)
      ENDDO
!
! *** CLOSE CHANNELS 16 AND 17 AND 18
!
      CLOSE (16)
      CLOSE (17)
      CLOSE (18)
    ENDIF

    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! *** FORMAT STATMENTS
!
500 FORMAT('RESULTS        1   ', I5, '    2      1      0')
600 FORMAT(I5, 4X, 3(1pe12.5, 1X))
!
    RETURN
  END

end module ADVNCE_MODULE
