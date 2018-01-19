         SUBROUTINE ADVNCE(NNODE,NGEOM,NBNOI,NBOUN,& 
     &       NBNOR,COORD,BSIDO,NBOUN_PP,NPOIN,NPOIN_PP,NELEM_PP,NELEM,& 
     &       NTIME,BSIDO_PP,INTMA_PP,& 
     &       FLUXP_PP,FLUYP_PP,& 
     &      GEOME_PP ,MMAT_PP  ,CMMAT_PP, RELEN ,DTE ,&         
     &      DELUN_PP ,RSIDO_PP ,UMEAN_PP ,IMMAT ,CSAFM ,UX,UY,& 
     &         NSIDE_PP,ISIDE_PP,NX_PP,NY_PP,EL_PP,UNKO_PP,&          
     &     VNPNT,SUMWEIGHT,DISNF_PP,rv,VCORD,&      
     &     RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,ALPHA,CINF,&   
     &     LCOMM_PP,GCOMM,NPGRP,IPCOM_PP,&   
     &   ISCOM_PP,MXCOM_PP,NCOMM_PP,maxNELEM_PP,RORDER,TORDER,&
     &   maxNPOIN_PP,SDCOM_PP,IVD,FORCEOUT,d,RGas,M,&
     &                     MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
     &                     MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &                     VSPACE_FIRST,VSPACE_LAST)
! 
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
      INTEGER NNODE,NGEOM,NSIDE_PP,NBNOI,NBNOR,NBOUN_PP 
      INTEGER RORDER,TORDER,TAG,NCOMM_PP,FORCEOUT
      INTEGER NTIME,ITIME,VNPNT,NPOIN_PP,NELEM_PP
      INTEGER MXCOM_PP,NGRPS,NPOIN,NBOUN 
      INTEGER LCOMM_PP(NSIDE_PP),NPGRP(MPI_SIZE_P-1),maxNPOIN_PP 
      INTEGER GCOMM(MPI_SIZE_P-1,MPI_SIZE_P-1),ISCOM_PP(NSIDE_PP),BSIDO(NBNOI,NBOUN) 
      INTEGER MPI_IERR,NELEM,IPCOM_PP(NPOIN_PP),MPI_STATUS(MPI_STATUS_SIZE) 
! 
      REAL RSIDO_PP(NBNOR,NBOUN_PP),R,THETA,ETA,ZETA,rv,PI 
      REAL GEOME_PP(NGEOM,maxNELEM_PP),COORD(2,NPOIN) 
      REAL VCORD(3,VNPNT),CINF(4),NFO_PP(NNODE,NELEM_PP) 
      REAL MMAT_PP(NNODE,maxNELEM_PP),CMMAT_PP(3,3,maxNELEM_PP)      
      REAL UX, UY, CSAFM,ALPHA,SPEED,C0,DTE,DTEc,DTEb,DTET
      REAL NX_PP(NSIDE_PP),NY_PP(NSIDE_PP),EL_PP(NSIDE_PP) 
      REAL RHO_PP(NPOIN_PP),UVEL_PP(NPOIN_PP),VVEL_PP(NPOIN_PP) 
      REAL TEMP_PP(NPOIN_PP),PS_PP(NPOIN_PP),ND_PP(NPOIN_PP) 
      REAL RHO(NPOIN),UVEL(NPOIN),VVEL(NPOIN),ND(NPOIN),PS(NPOIN) 
      REAL TEMP(NPOIN),RELEN(NELEM),MU_PP(NNODE,NELEM_PP)
      INTEGER VNPNT_PART,VSPACE_FIRST,VSPACE_LAST,PRINT_EVERY
      REAL DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,maxNELEM_PP)
      REAL RHS_PP(1,3,NELEM_PP) 
      REAL FLUXP_PP(1,3,NELEM_PP),FLUYP_PP(1,3,NELEM_PP) 
      REAL UMEAN_PP(NELEM_PP),DELUN_PP(1,3,NELEM_PP),ETA_PP(NBOUN_PP)
      REAL DISND_PP(NNODE,NELEM_PP),DISUX_PP(NNODE,NELEM_PP)
      REAL DISUY_PP(NNODE,NELEM_PP),DISPS_PP(NNODE,NELEM_PP)
      REAL SUMWEIGHT,DIFF,DTETEST(2),SRITM1
      REAL SUMRES,SUMRESG

      REAL RESIDUAL(VSPACE_FIRST:VSPACE_LAST) 
      REAL RESIDUAL_PP(VSPACE_FIRST:VSPACE_LAST) 
      REAL LIFT,DRAG,d,RGas,M
      REAL RHS_PP_CP(1,3,NELEM_PP), SUPR,L2R,L1R, DENOM ! TESTING
      REAL RHS_PP_D(1,3,NELEM_PP), TMP!TESTING
      INTEGER I,J,NONZERODIFF1,NONZERODIFF2 ! TESTING
! 
      INTEGER INTMA_PP(NNODE,NELEM_PP),BSIDO_PP(NBNOI,NBOUN_PP) 
      INTEGER ISIDE_PP(8,NSIDE_PP),IMMAT,SDCOM_PP(3,NCOMM_PP) 
      ! mpi-stuff for position space partitioning
      INTEGER MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,GROUP_P
      ! mpi-stuff for velocity space partitioning
      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,GROUP_V
      INTEGER :: MPI_COMM_SLAVES ! SLAVES COMMUNICATOR
      INTEGER :: COLOR,SLAVERANK 

      CHARACTER filename*80
!      CHARACTER dbgfilename*80 !TESTING-DEBUG
!
      TYPE(InputVariables) :: IVD  
!
      PARAMETER (PI=3.1416)
      PRINT_EVERY = (VSPACE_LAST+1-VSPACE_FIRST)/5
! 
! *** OPEN CHANNEL 16 TO WRITE TO THE RESIDUAL FILE 
! 
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN ! write only on master rank
        filename = IVD%ResidualFile
!
        PRINT*,'OPENING RESIDUAL FILE = ',IVD%ResidualFile
        OPEN  (16, file=filename, status='UNKNOWN') 
        WRITE(*,*) 
!  
! ***   SET UP CHANNEL 17 & 18 TO WRITE TO RESULTS FILES FOR GiD POST PROCESSING 
! 
        filename = IVD%ResultsFile1
        PRINT*,'OPENING 1st RESULTS FILE = ',IVD%ResultsFile1
        OPEN  (17, file=filename, status='UNKNOWN')  
        WRITE(*,*) 
        WRITE(*,*) 
        filename = IVD%ResultsFile2
        PRINT*,'OPENING 2nd RESULTS FILE = ',IVD%ResultsFile2
        OPEN (18, file=filename,status='UNKNOWN')
        WRITE(*,*) 
      ENDIF !dskljdasda
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
!
! *** CODE TO DO THE DISTRIBUTION FUNCTION PRINTOUTS AT VARIOUS MESH POINTS
!

! 
      C0=0.0 
! 
! *** GET THE COURANT CONDITION MAX ALLOWABLE TIMESTEP 
!       
      IF(MPI_RANK_P.EQ.0) CALL ALOTIM(NELEM ,CSAFM ,RELEN ,DTEc,rv ) 
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
! 
! *** CALL MACROS TO CONVERT nf'S INTO MACROSCOPIC QUANTITIES 
!
      IF(MPI_RANK_P.NE.0)THEN
        CALL MACROS(rv,VNPNT,SUMWEIGHT,NPOIN_PP,NNODE,DISNF_PP,&
     &           VCORD,INTMA_PP,NELEM_PP,GEOME_PP,NGEOM,&
     &           ND_PP,RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,&
     &           DISND_PP,DISUX_PP,DISUY_PP,DISPS_PP,M,RGas,&
     &           MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &           VSPACE_FIRST,VSPACE_LAST)
      ENDIF
! 
! *** CONVERT MACRO_PPs to MACROS 
! 
      NGRPS=MPI_SIZE_P-1 

      CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,&
     &          UVEL_PP,VVEL_PP,& 
     &          PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,& 
     &          NGRPS,MPI_RANK_P,MPI_COMM_P)
! 
! *** WRITE OUTPUT DATA TO RESULTS FILE 
!  
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN ! write only on master rank
        WRITE(17,500) 0 
        WRITE(18,500) 0 
          DO IP=1,NPOIN 
            WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)   
            WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP) 
          ENDDO
! 
! *** LOOP OVER THE PRESCRIBED NR. OF TIMESTEPS 
! 
        PRINT*,'NTIME= ',NTIME 
        WRITE(*,*) 
      ENDIF           !LINKS WITH IF STATEMENT ON LINE 131 
! 
! *** BROADCAST NTIME TO ALL PROCESSORS 
! 
      CALL MPI_BCAST(NTIME,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR) 
! 
! *** BEGIN THE TIMESTEPPING !!! 
! 

      IF (MPI_RANK_P.NE.0) THEN ! adadadfadfad
       ! To create new communicator 
            COLOR = 1
            SLAVERANK = MPI_RANK_P-1
      ELSE  ! IF (MPI_RANK_P.NE.0) ! adadadfadfad
            COLOR = MPI_UNDEFINED
            SLAVERANK = 0
      ENDIF ! IF (MPI_RANK_P.NE.0) ! adadadfadfad
 
      CALL MPI_COMM_SPLIT(MPI_COMM_P,COLOR,SLAVERANK,&
     &                          MPI_COMM_SLAVES,MPI_IERR)
      IF(MPI_RANK_P.NE.0)THEN !kdjfhsewew
        CALL MPI_COMM_RANK(MPI_COMM_SLAVES,SLAVERANK,MPI_IERR)
      ENDIF

  

      ITERTIME = MPI_WTIME()
      DO ITIME=1,NTIME !fsajfdclaksjdnckajlndaa
        IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0)) THEN
          PRINT*,'TIMESTEP = ',ITIME ! PRINT TIMESTEP NUMBER TO SCREEN
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 

! 
! ***   CALCULATE THE FLUX CONSERVATION PARAMETER (ETA) AT EACH BOUNDARY NODE 
!
        IF(MPI_RANK_P.NE.0)THEN
          CALL FLUCON2(NELEM_PP,VNPNT,DISNF_PP,NBNOI,NBOUN_PP,& 
     &          BSIDO_PP,NBNOR,RSIDO_PP,VCORD,rv,& 
     &          ETA_PP,NSIDE_PP,ISIDE_PP,RORDER,TORDER,SUMWEIGHT,RGas,&
     &          MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &          VSPACE_FIRST,VSPACE_LAST)

          ENDIF
!
! ***   COMPUTE THE VALUE OF MU FOR CONSTRUCTION OF THE BGK RHS SOURCE TERM
!
        CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
        IF(MPI_RANK_P.NE.0)THEN
          CALL COMPMU(rv,VNPNT,NPOIN_PP,NNODE,DISNF_PP,& 
     &            VCORD,NELEM_PP,GEOME_PP,NGEOM,&
     &          MU_PP,UX,UY,DISND_PP,DISUX_PP,DISUY_PP,SUMWEIGHT,d,&
     &          MPI_COMM_V,VSPACE_FIRST,VSPACE_LAST)
!
! ***   COMPUTE THE MAX TIMESTEP BASED ON THE BGK TERM
!
          CALL ALOTIM2(NNODE,NELEM_PP,MU_PP,DTEb)
        ENDIF
!
! ***   LOOP OVER THE PROCESSORS TO FIND THE LIMITING BGK TIMESTEP
!
        CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)

        CALL MPI_ALLREDUCE(DTEb,DTE,1,MPI_REAL,MPI_MIN,&
     &                  MPI_COMM_P,MPI_IERR)
        IF((MPI_RANK_P.EQ.0)) THEN
          DTEb = DTE
          DTETEST(1)=DTEb
          DTETEST(2)=DTEc
          DTET=MINVAL(DTETEST)
          PRINT*,'RANK_V=',MPI_RANK_V,&
     &            'DTEb=',DTEb,'DTEc=',DTEc,'DTE=',DTE

          CALL MPI_ALLREDUCE(DTET,DTE,1,MPI_REAL,MPI_MIN,&
     &                  MPI_COMM_V,MPI_IERR)
          PRINT*,'AFTER REDUCE RANK_V=',MPI_RANK_V,&
     &            'DTEb=',DTEb,'DTEc=',DTEc,'DTE=',DTE
        ENDIF  
!
        CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
        CALL MPI_BCAST(DTE,1,MPI_REAL,0,MPI_COMM_P,MPI_IERR) 
!
!
! ***   LOOP OVER THE NODES IN VELOCITY SPACE 
! 
        DO 7000 IV=VSPACE_FIRST,VSPACE_LAST
          IF((MPI_RANK_P.EQ.0).AND.((MOD(IV,PRINT_EVERY)).EQ.0))THEN 
            WRITE(*,"(A3,I2,A25,I5)"),'VSR',MPI_RANK_V,&
     &               ": V-space iteration no.", IV
          ENDIF

! ***     SYNCHRONISE 
! 
          CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
! 
! ***     IDENTIFY THE POSITION IN THE VELOCITY SPACE MESH BY UX,UV 
! 

          ETA=VCORD(1,IV)
          ZETA=VCORD(2,IV)
          R=(rv/2)*(ETA+1)
          THETA=ZETA*PI
          UX=R*COS(THETA)
          UY=R*SIN(THETA)
          SPEED=SQRT(UX*UX+UY*UY)

! ***     COMPUTE THE VALUE OF THE EQUILIBRIUM DIST FUNC FOR THIS POINT IN V-SPACE
! ***     BASED ON THE BULK CONDITIONS
!
          IF(MPI_RANK_P.NE.0)THEN
            CALL GTEQNF(NFO_PP,NNODE,NELEM_PP,DISND_PP,&
     &          DISUX_PP,DISUY_PP,DISPS_PP,UX,UY,MPI_RANK_P,M)
          ENDIF
!
          IF(MPI_RANK_P.NE.0)THEN !dvjkaapxxmxkayy
! 
! ***       OBTAIN THE FLUXES AT THE POINTS 
! 
            CALL GETFLA(NELEM_PP,NNODE,VNPNT,& 
     &         DISNF_PP,FLUXP_PP,FLUYP_PP,&
     &         UX,UY,IV,VSPACE_FIRST,VSPACE_LAST)
! 
! ***       ELEMENT CONTRIBUTIONS (STEP 1) 
! 
            CALL GETELC(NELEM_PP,GEOME_PP,NGEOM,NNODE,VNPNT,IV,DISNF_PP,& 
     &        FLUXP_PP,FLUYP_PP,UMEAN_PP,RHS_PP,DTE,UX,UY,MU_PP,NFO_PP,&
     &        ITIME,VSPACE_FIRST,VSPACE_LAST)
! 
! ***       FILL DELUN(NAMAT=1,NPOIN) WITH THE VALUE SET IN C00 
            CALL RFILLA(DELUN_PP,1,NNODE,NELEM_PP,C0) 
          ENDIF  !  IF(MPI_RANK_P.NE.0)THEN !dvjkaapxxmxkayy
          CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)


          CALL EDGFLXOPT(NELEM_PP,NSIDE_PP,ISIDE_PP,RHS_PP,& 
     &       NX_PP,NY_PP,EL_PP,UX,UY,RSIDO_PP,BSIDO_PP,NBOUN_PP,& 
     &       NBNOR,ALPHA,ETA_PP,VNPNT,IV,DISNF_PP,UMEAN_PP,& 
     &       CINF,rv,LCOMM_PP,NGRPS,ISCOM_PP,MPI_RANK_P,NCOMM_PP,VCORD,&
     &       RORDER,TORDER,SDCOM_PP,RGas,M,ITIME,GCOMM,MPI_COMM_P,&
     &       VSPACE_FIRST,VSPACE_LAST,MPI_COMM_SLAVES,SLAVERANK)

          IF(MPI_RANK_P.NE.0)THEN !dssserwfjvvfskjfs
! ***       MULTIPLY BY DELTA-TIME-POINT 
            DO 4500 IE=1,NELEM_PP 
              DO 4501 IP=1,NNODE      
                RHS_PP(1,IP,IE)=DTE*RHS_PP(1,IP,IE) 
 4501         CONTINUE 
 4500       CONTINUE
! 
! *** OBTAIN THE INCREMENTS -PREDICTIONS 
! 
            IF(IMMAT.EQ.1)THEN 
              CALL GETINC(NNODE ,MMAT_PP,RHS_PP,DELUN_PP,& 
     &           NELEM_PP,RESIDUAL_PP,maxNELEM_PP,VNPNT,DISNF_PP,&
     &           UMEAN_PP,NFO_PP,MU_PP,IV,DTE,&
     &           VSPACE_FIRST,VSPACE_LAST)
            ELSE 
              CALL CMMINC(NNODE,CMMAT_PP,RHS_PP,DELUN_PP,& 
     &           NELEM_PP,RESIDUAL_PP,maxNELEM_PP,VNPNT,DISNF_PP,&
     &           UMEAN_PP,NFO_PP,MU_PP,IV,DTE,&
     &           VSPACE_FIRST,VSPACE_LAST) 
            ENDIF
! 
! *** ADD INCREMENTS  
! 
            CALL ADTHEM(IV,VNPNT,DISNF_PP,DELUN_PP,NELEM_PP,&
     &        VSPACE_FIRST,VSPACE_LAST)
!
! *** RESET THE INFLOW BOUNDARY VALUES
!
            CALL INFLOW(IV,UX,UY,NNODE,VNPNT,NELEM_PP,NPOIN_PP,&
     &        INTMA_PP,CINF,NBOUN_PP,NBNOI,&
     &        BSIDO_PP,DISNF_PP,ITIME,PS_PP,RGas,M,&
     &        VSPACE_FIRST,VSPACE_LAST)

!
          ENDIF ! IF(MPI_RANK_P.NE.0)THEN !dssserwfjvvfskjfs
! 
          CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
! 
! *** END LOOP OVER VELOCITY SPACE NODES IN THE CURRENT PARTITION
! 
 7000   CONTINUE

        CALL GETRES(VNPNT,RESIDUAL,RESIDUAL_PP,MPI_RANK_P,MPI_SIZE_P,&
     &           VSPACE_FIRST,VSPACE_LAST,MPI_COMM_P) 

! 
! *** CALL MACROS 
!
        IF(MPI_RANK_P.NE.0)THEN 
          CALL MACROS(rv,VNPNT,SUMWEIGHT,NPOIN_PP,NNODE,DISNF_PP,&
     &           VCORD,INTMA_PP,NELEM_PP,GEOME_PP,NGEOM,&
     &           ND_PP,RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,&
     &           DISND_PP,DISUX_PP,DISUY_PP,DISPS_PP,M,RGas,&
     &           MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &           VSPACE_FIRST,VSPACE_LAST)
        ENDIF

        IF(FORCEOUT.EQ.0)THEN !NOT WRITING OUT FORCES TO RESIDUAL FILE
! *** CONSTRUCT GLOBAL MACRO VECTORS
          IF (MOD(ITIME,1000).EQ.0) THEN
            CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,&
     &        UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,&
     &        PS,TEMP,NGRPS,MPI_RANK_P,MPI_COMM_P)
! ***       WRITE OUTPUT DATA TO RESULTS FILE 
            IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN
! ***         CONSTRUCT GLOBAL MACRO VARIABLE VECTORS 
              WRITE(17,500) ITIME 
              WRITE(18,500) ITIME 
              DO 3000 IP=1,NPOIN 
                WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)   
                WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP) 
 3000         CONTINUE
            ENDIF
          ENDIF
        ELSE  !WRITING OUT FORCES TO RESIDUAL FILE
! ***     CONSTRUCT GLOBAL MACRO VECTORS
          CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,&
     &        UVEL_PP,VVEL_PP, PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,&
     &        PS,TEMP,NGRPS,MPI_RANK_P,MPI_COMM_P)

          IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0)&
     &         .AND.(MOD(ITIME,1000).EQ.0))THEN
            WRITE(17,500) ITIME
            WRITE(18,500) ITIME
            DO 3011 IP=1,NPOIN
              WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)
              WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP)
 3011       CONTINUE
          ENDIF
        ENDIF
! 
! *** WRITE TO MAX RESIDUAL FILE 
! 
        IF(MPI_RANK_P.EQ.0)THEN
          SUMRES=MAXVAL(RESIDUAL) 
          CALL MPI_REDUCE(SUMRES,SUMRESG,1,MPI_REAL,MPI_MAX,&
     &                             0,MPI_COMM_V,MPI_IERR)
          SUMRE=SUMRESG
          if(ITIME.EQ.1)SRITM1=SUMRES
          IF(MPI_RANK_V.EQ.0) THEN
            IF(IVD%FORCEOUT.EQ.0)THEN     !IF NOT OUTPUTTING FORCES
              WRITE(16,*) ITIME,LOG(SUMRES/SRITM1)
            ELSE             !IF OUTPUTTING FORCES IN RESIDUAL FILE
              CALL GETFOR(LIFT,DRAG,NBOUN,NPOIN,NBNOI,PS,COORD,BSIDO)
              WRITE(16,*) ITIME,LOG(SUMRES/SRITM1),LIFT,DRAG
            ENDIF
          ENDIF
        ENDIF 
! 
! ***   SYNCHRONISE 
! 
        CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
! 
! ***   END IF TIMESTEPPING LOOP 
! 
      ENDDO ! DO ITIME=1,NTIME !fsajfdclaksjdnckajlndaa
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN 
        WRITE(*,*) "END OF TIME LOOP"
      ENDIF
!
! *** CONSTRACT GLOBAL MACRO VECTORS
!
!        CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,UVEL_PP,VVEL_PP,& 
!    &         PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,& 
!    &           NGRPS,MPI_RANK_P,MPI_COMM_P)
! 
! *** WRITE OUTPUT DATA TO RESULTS FILE 
!  
      IF((MPI_RANK_P.EQ.0).AND.(MPI_RANK_V.EQ.0))THEN
! *** CONSTRUCT GLOBAL MACRO VARIABLE VECTORS 
        WRITE(17,500) ITIME 
        WRITE(18,500) ITIME 
        DO 3002 IP=1,NPOIN 
          WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)   
          WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP) 
 3002   CONTINUE
! 
! *** CLOSE CHANNELS 16 AND 17 AND 18 
! 
        CLOSE(16) 
        CLOSE(17) 
        CLOSE(18) 
      ENDIF 


      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
!  
! *** FORMAT STATMENTS 
! 
  300 FORMAT('WHAT WOULD YOU LIKE THE RESIDUAL FILE TO BE CALLED?& 
     & (____________.RES  :    ',$) 
  400 FORMAT('WHAT WOULD YOU LIKE THE RESULTS FILE 1 TO BE CALLED?& 
     & (____________.RES  :    ', $) 
  410 FORMAT('WHAT WOULD YOU LIKE THE RESULTS FILE 2 TO BE CALLED?& 
     & (____________.RES  :    ', $) 
  500 FORMAT('RESULTS        1   ',I5,'    2      1      0') 
  600 FORMAT(I5,4X,3(1pe12.5,1X)) 
! 
      RETURN 
      END 
