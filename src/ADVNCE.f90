         SUBROUTINE ADVNCE(NNODE,NGEOM,NBNOI,NBOUN,& 
     &       NBNOR,COORD,BSIDO,NBOUN_PP,NPOIN,NPOIN_PP,NELEM_PP,NELEM,& 
     &       NTIME,BSIDO_PP,INTMA_PP,& 
     &       FLUXP_PP,FLUYP_PP,& 
     &      GEOME_PP ,MMAT_PP  ,CMMAT_PP, RELEN ,DTE ,&         
     &      DELUN_PP ,RSIDO_PP ,UMEAN_PP ,IMMAT ,CSAFM ,UX,UY,& 
     &         NSIDE_PP,ISIDE_PP,NX_PP,NY_PP,EL_PP,UNKO_PP,&          
     &     VNPNT,SUMWEIGHT,DISNF_PP,rv,VCORD,&      
     &     RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,ALPHA,CINF,&   
     &   RANK,NPROC, LCOMM_PP,GCOMM,NPGRP,IPCOM_PP,&   
     &   ISCOM_PP,MXCOM_PP,NCOMM_PP,maxNELEM_PP,RORDER,TORDER,maxNPOIN_PP,&
     &     SDCOM_PP,IVD,FORCEOUT,d,RGas,M)
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
      INTEGER RANK,NPROC,MXCOM_PP,NGRPS,NPOIN,NBOUN 
      INTEGER LCOMM_PP(NSIDE_PP),NPGRP(NPROC-1),maxNPOIN_PP 
      INTEGER GCOMM(NPROC-1,NPROC-1),ISCOM_PP(NSIDE_PP),BSIDO(NBNOI,NBOUN) 
      INTEGER MPI_IERR,NELEM,IPCOM_PP(NPOIN_PP),MPI_STATUS(MPI_STATUS_SIZE) 
! 
      REAL RSIDO_PP(NBNOR,NBOUN_PP),R,THETA,ETA,ZETA,rv,PI 
      REAL GEOME_PP(NGEOM,maxNELEM_PP),COORD(2,NPOIN) 
      REAL VCORD(3,VNPNT),CINF(4),NFO_PP(NNODE,NELEM_PP) 
      REAL MMAT_PP(NNODE,maxNELEM_PP),CMMAT_PP(3,3,maxNELEM_PP)      
      REAL UX, UY, CSAFM,ALPHA,SPEED,C0,DTE,DTEc,DTEb,DTEbnew,RESIDUAL_PP 
      REAL NX_PP(NSIDE_PP),NY_PP(NSIDE_PP),EL_PP(NSIDE_PP) 
      REAL RHO_PP(NPOIN_PP),UVEL_PP(NPOIN_PP),VVEL_PP(NPOIN_PP) 
      REAL TEMP_PP(NPOIN_PP),PS_PP(NPOIN_PP),ND_PP(NPOIN_PP) 
      REAL RHO(NPOIN),UVEL(NPOIN),VVEL(NPOIN),ND(NPOIN),PS(NPOIN) 
      REAL TEMP(NPOIN),RELEN(NELEM),MU_PP(NNODE,NELEM_PP)
      REAL DISNF_PP(NNODE,VNPNT,NELEM_PP),RHS_PP(1,3,NELEM_PP) 
      REAL FLUXP_PP(1,3,NELEM_PP),FLUYP_PP(1,3,NELEM_PP) 
      REAL UMEAN_PP(NELEM_PP),DELUN_PP(1,3,NELEM_PP),ETA_PP(NBOUN_PP)
      REAL DISND_PP(NNODE,NELEM_PP),DISUX_PP(NNODE,NELEM_PP)
      REAL DISUY_PP(NNODE,NELEM_PP),DISPS_PP(NNODE,NELEM_PP)
      REAL SUMWEIGHT,DIFF,DTETEST(2),RESIDUAL(VNPNT),SUMRES,SRITM1 
      REAL LIFT,DRAG,d,RGas,M
! 
      INTEGER INTMA_PP(NNODE,NELEM_PP),BSIDO_PP(NBNOI,NBOUN_PP) 
      INTEGER ISIDE_PP(8,NSIDE_PP),IMMAT,SDCOM_PP(3,NCOMM_PP) 
! 
      CHARACTER filename*80
!
      TYPE(InputVariables) :: IVD  
!
      PARAMETER (PI=3.1416)
! 
! *** OPEN CHANNEL 16 TO WRITE TO THE RESIDUAL FILE 
! 
      IF(RANK.EQ.0)THEN 
      filename = IVD%ResidualFile
!
      PRINT*,'OPENING RESIDUAL FILE = ',IVD%ResidualFile
      OPEN  (16, file=filename, status='UNKNOWN') 
      WRITE(*,*) 
!  
! *** SET UP CHANNEL 17 & 18 TO WRITE TO RESULTS FILES FOR GiD POST PROCESSING 
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
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
!
! *** CODE TO DO THE DISTRIBUTION FUNCTION PRINTOUTS AT VARIOUS MESH POINTS
!
!      IF(RANK.EQ.3)THEN
!      filename='nfplot1.dat'
!      OPEN(15,file=filename,status='NEW')
!      filename='nfplot2.dat'
!      OPEN(25,file=filename,status='NEW')
!      filename='nfplot3.dat'
!      OPEN(35,file=filename,status='NEW')
!      filename='nfplot4.dat'
!      OPEN(45,file=filename,status='NEW')
!      filename='nfplot5.dat'
!      OPEN(55,file=filename,status='NEW')
!      filename='nfplot6.dat'
!      OPEN(65,file=filename,status='NEW')
!      filename='nfplot7.dat'
!      OPEN(75,file=filename,status='NEW')
!      filename='nfplot8.dat'
!      OPEN(85,file=filename,status='NEW')
!      ENDIF

! 
      C0=0.0 
! 
! *** GET THE COURANT CONDITION MAX ALLOWABLE TIMESTEP 
!       
      IF(RANK.EQ.0)THEN 
      CALL ALOTIM(NELEM ,CSAFM ,RELEN ,DTEc,rv ) 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
! 
! *** CALL MACROS TO CONVERT nf'S INTO MACROSCOPIC QUANTITIES 
!
      IF(RANK.NE.0)THEN
      CALL MACROS(rv,VNPNT,SUMWEIGHT,NPOIN_PP,NNODE,DISNF_PP,&
     &           VCORD,INTMA_PP,NELEM_PP,GEOME_PP,NGEOM,&
     &           ND_PP,RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,RANK,&
     &           DISND_PP,DISUX_PP,DISUY_PP,DISPS_PP,M,RGas)
      ENDIF
! 
! *** CONVERT MACRO_PPs to MACROS 
! 
       NGRPS=NPROC-1 
       CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,UVEL_PP,VVEL_PP,& 
     &          PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,& 
     &          NGRPS,RANK)
! 
! *** WRITE OUTPUT DATA TO RESULTS FILE 
!  
      IF(RANK.EQ.0)THEN 
      WRITE(17,500) 0 
      WRITE(18,500) 0 
      DO 3001 IP=1,NPOIN 
      WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)   
      WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP) 
 3001 CONTINUE
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
      DO 30000 ITIME=1,NTIME
! 
         IF(RANK.EQ.0) THEN
           PRINT*,'TIMESTEP = ',ITIME ! PRINT TIMESTEP NUMBER TO SCREEN
         ENDIF
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 

! 
! *** CALCULATE THE FLUX CONSERVATION PARAMETER (ETA) AT EACH BOUNDARY NODE 
!
         IF(RANK.NE.0)THEN       !SLAVE PROCESSOR ONLY 
           CALL FLUCON2(NELEM_PP,VNPNT,DISNF_PP,NBNOI,NBOUN_PP,& 
     &               BSIDO_PP,NBNOR,RSIDO_PP,VCORD,rv,& 
     &   ETA_PP,NSIDE_PP,ISIDE_PP,RANK,RORDER,TORDER,SUMWEIGHT,RGas)
         ENDIF         !LINKS WITH IF STATEMENT ON LINE 154
!
! *** COMPUTE THE VALUE OF MU FOR CONSTRUCTION OF THE BGK RHS SOURCE TERM
!
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
         IF(RANK.NE.0)THEN
           CALL COMPMU(rv,VNPNT,NPOIN_PP,NNODE,DISNF_PP,& 
     &             VCORD,INTMA_PP,NELEM_PP,GEOME_PP,NGEOM,&
     &           MU_PP,UX,UY,DISND_PP,DISUX_PP,DISUY_PP,RANK,SUMWEIGHT,d)
!
! *** COMPUTE THE MAX TIMESTEP BASED ON THE BGK TERM
!
           CALL ALOTIM2(NNODE,NELEM_PP,MU_PP,DTEb)
         ENDIF
!
! *** LOOP OVER THE PROCESSORS TO FIND THE LIMITING BGK TIMESTEP
!
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
        IF(RANK.EQ.0)DTEb=100000
          DO 7050 IG=1,NGRPS
           TAG=10*IG
           IF(RANK.EQ.IG)THEN
             CALL MPI_SEND(DTEb,1,MPI_REAL,0,TAG,MPI_COMM_WORLD,MPI_IERR)
           ENDIF
           IF(RANK.EQ.0)THEN
             CALL MPI_RECV(DTEbnew,1,MPI_REAL,IG,TAG,MPI_COMM_WORLD,MPI_STATUS,MPI_IERR)
             IF(DTEbnew.LT.DTEb)DTEb=DTEbnew
           ENDIF
 7050 CONTINUE
        IF(RANK.EQ.0)THEN
          DTETEST(1)=DTEb
          DTETEST(2)=DTEc
          DTE=MINVAL(DTETEST)
          PRINT*,'DTEb=',DTEb,'DTEc=',DTEc,'DTE=',DTE
        ENDIF         !LINKS WITH IF STATEMENT ON LINE 169
!
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
        CALL MPI_BCAST(DTE,1,MPI_REAL,0,MPI_COMM_WORLD,MPI_IERR) 
!
!
! *** LOOP OVER THE NODES IN VELOCITY SPACE 
! 
        DO 7000 IV=1,VNPNT
! 
! *** SYNCHRONISE 
! 
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
!         PRINT *,'4-GOT HERE',RANK,VNPNT
!         FLUSH(6)
! 
! *** IDENTIFY THE POSITION IN THE VELOCITY SPACE MESH BY UX,UV 
! 
         ETA=VCORD(1,IV)
         ZETA=VCORD(2,IV)
         R=(rv/2)*(ETA+1)
         THETA=ZETA*PI
         UX=R*COS(THETA)
         UY=R*SIN(THETA)
         SPEED=SQRT(UX*UX+UY*UY)
! *** PERFORM THE DISTFUNC.PLOT TEST PRINTOUT
!
!        IF((ITIME.EQ.1).OR.(ITIME.EQ.100).OR.(ITIME.EQ.200).OR.(ITIME.EQ.300).OR.(ITIME.EQ.400).OR.(ITIME.EQ.500).OR.(ITIME.EQ.600).OR.(ITIME.EQ.700).OR.(ITIME.EQ.800).OR.(ITIME.EQ.900).OR.(ITIME.EQ.1000))THEN
!         IF(RANK.EQ.3)THEN
!           WRITE(15,*) UX,UY,DISNF_PP(1,IV,1616)
!           WRITE(25,*) UX,UY,DISNF_PP(2,IV,2017)
!           WRITE(35,*) UX,UY,DISNF_PP(2,IV,2640)
!           WRITE(45,*) UX,UY,DISNF_PP(2,IV,2964)
!           WRITE(55,*) UX,UY,DISNF_PP(2,IV,2126)
!           WRITE(65,*) UX,UY,DISNF_PP(2,IV,2147)
!           WRITE(75,*) UX,UY,DISNF_PP(2,IV,2549)
!           WRITE(85,*) UX,UY,DISNF_PP(2,IV,3374)
!         ENDIF
!         ENDIF

! *** COMPUTE THE VALUE OF THE EQUILIBRIUM DIST FUNC FOR THIS POINT IN V-SPACE
! *** BASED ON THE BULK CONDITIONS
!
        IF(RANK.NE.0)THEN
          CALL GTEQNF(NFO_PP,NNODE,NELEM_PP,DISND_PP,DISUX_PP,DISUY_PP,&
     &           DISPS_PP,VNPNT,VCORD,UX,UY,RANK,M)
        ENDIF
!
        IF(RANK.NE.0)THEN             !SLAVE PROCESSORS ONLY 
! 
! *** OBTAIN THE FLUXES AT THE POINTS 
! 
        CALL GETFLA(NELEM_PP,NNODE,VNPNT,& 
     &         DISNF_PP,FLUXP_PP,FLUYP_PP,UX,UY,IV)
! 
! *** ELEMENT CONTRIBUTIONS (STEP 1) 
! 
      CALL GETELC(NELEM_PP,GEOME_PP,NGEOM,NNODE,VNPNT,IV,DISNF_PP,& 
     &      FLUXP_PP,FLUYP_PP,UMEAN_PP,RHS_PP,DTE,UX,UY,MU_PP,NFO_PP,&
     &      RANK,ITIME)
! 
! *** FILL DELUN(NAMAT=1,NPOIN) WITH THE VALUE SET IN C00 
      CALL RFILLA(DELUN_PP,1,NNODE,NELEM_PP,C0) 
       ENDIF !LINKS WITH IF STATEMENT ON LINE ?
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
         PRINT *,'5-GOT HERE',RANK,ITIME,IV, UX, UY
         FLUSH(6)

      CALL EDGFLXOPT(NELEM_PP,NSIDE_PP,ISIDE_PP,RHS_PP,& 
     &           NX_PP,NY_PP,EL_PP,UX,UY,RSIDO_PP,BSIDO_PP,NBOUN_PP,& 
     &           NBNOR,ALPHA,ETA_PP,VNPNT,IV,DISNF_PP,UMEAN_PP,& 
     &        CINF,rv,LCOMM_PP,NGRPS,ISCOM_PP,RANK,NCOMM_PP,VCORD,RORDER,&
     &         TORDER,SDCOM_PP,RGas,M,ITIME,GCOMM)

!      CALL EDGFLX(NELEM_PP,NSIDE_PP,ISIDE_PP,RHS_PP,& 
!     &           NX_PP,NY_PP,EL_PP,UX,UY,RSIDO_PP,BSIDO_PP,NBOUN_PP,& 
!     &           NBNOR,ALPHA,ETA_PP,VNPNT,IV,DISNF_PP,UMEAN_PP,& 
!     &        CINF,rv,LCOMM_PP,NGRPS,ISCOM_PP,RANK,NCOMM_PP,VCORD,RORDER,&
!     &         TORDER,SDCOM_PP,RGas,M,ITIME)
         PRINT *,'6-GOT HERE',RANK
         FLUSH(6)

       IF(RANK.NE.0)THEN
! 
! *** MULTIPLY BY DELTA-TIME-POINT 
! 
        DO 4500 IE=1,NELEM_PP 
        DO 4501 IP=1,NNODE      
        RHS_PP(1,IP,IE)=DTE*RHS_PP(1,IP,IE) 
 4501 CONTINUE 
 4500 CONTINUE
! 
! *** OBTAIN THE INCREMENTS -PREDICTIONS 
! 
        IF(IMMAT.EQ.1)THEN 
        CALL GETINC(NNODE ,MMAT_PP,RHS_PP,DELUN_PP,& 
     &         NELEM_PP,RESIDUAL_PP,maxNELEM_PP,VNPNT,DISNF_PP,&
     &          UMEAN_PP,NFO_PP,MU_PP,IV,DTE)
        ELSE 
        CALL CMMINC(NNODE,CMMAT_PP,RHS_PP,DELUN_PP,& 
     &         NELEM_PP,RESIDUAL_PP,maxNELEM_PP,VNPNT,DISNF_PP,&
     &          UMEAN_PP,NFO_PP,MU_PP,IV,DTE) 
        ENDIF
! 
! *** ADD INCREMENTS  
! 
      CALL ADTHEM(IV,VNPNT,DISNF_PP,DELUN_PP,NELEM_PP)
!
! *** RESET THE INFLOW BOUNDARY VALUES
!
      CALL INFLOW(IV,UX,UY,NNODE,VNPNT,NELEM_PP,NPOIN_PP,INTMA_PP,&
     &    CINF,NBOUN_PP,NBNOI,BSIDO_PP,DISNF_PP,ITIME,PS_PP,RGas,M)
!
      ENDIF !LINK WITH IF STATEMENT ON LINE 241
! 
!         PRINT *,'7-GOT HERE',RANK
!         FLUSH(6)
! *** GETRES SUBROUTINE HERE TO GET THE MAX RESIDUAL FROM ALL PROCS 
! 
       CALL GETRES(IV,VNPNT,RESIDUAL,RESIDUAL_PP,RANK,NPROC) 
! 
!         PRINT *,'8-GOT HERE',RANK
!         FLUSH(6)
! *** SYNCHRONISE 
! 
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
!         PRINT *,'9-GOT HERE',RANK
!         FLUSH(6)
! 
! *** END LOOP OVER VELOCITY SPACE NODES 
! 
 7000 CONTINUE
! 
! *** CALL MACROS 
!
       IF(RANK.NE.0)THEN 
        CALL MACROS(rv,VNPNT,SUMWEIGHT,NPOIN_PP,NNODE,DISNF_PP,&
     &           VCORD,INTMA_PP,NELEM_PP,GEOME_PP,NGEOM,&
     &           ND_PP,RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,RANK,&
     &           DISND_PP,DISUX_PP,DISUY_PP,DISPS_PP,M,RGas)
       ENDIF

      IF(FORCEOUT.EQ.0)THEN !NOT WRITING OUT FORCES TO RESIDUAL FILE
        IF((ITIME.EQ.1000).OR.(ITIME.EQ.2000).OR.(ITIME.EQ.3000).OR.(ITIME.EQ.4000)&
     &  .OR.(ITIME.EQ.5000).OR.(ITIME.EQ.6000).OR.(ITIME.EQ.7000)&
     &  .OR.(ITIME.EQ.8000).OR.(ITIME.EQ.9000).OR.(ITIME.EQ.10000)&
     &  .OR.(ITIME.EQ.11000).OR.(ITIME.EQ.12000).OR.(ITIME.EQ.13000)&
     &  .OR.(ITIME.EQ.14000).OR.(ITIME.EQ.15000).OR.(ITIME.EQ.20000)&
     &  .OR.(ITIME.EQ.25000).OR.(ITIME.EQ.30000).OR.(ITIME.EQ.35000)&
     &   .OR.(ITIME.EQ.40000).OR.(ITIME.EQ.45000).OR.(ITIME.EQ.50000))THEN
             PRINT *,'WRITING1',RANK,ITIME
             FLUSH(6)
!
! *** CONSTRACT GLOBAL MACRO VECTORS
!
          CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,UVEL_PP,VVEL_PP,& 
     &         PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,& 
     &           NGRPS,RANK)
! 
! *** WRITE OUTPUT DATA TO RESULTS FILE 
!  
           IF(RANK.EQ.0)THEN
! *** CONSTRUCT GLOBAL MACRO VARIABLE VECTORS 
             WRITE(17,500) ITIME 
             WRITE(18,500) ITIME 
             DO 3000 IP=1,NPOIN 
               WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)   
               WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP) 
 3000        CONTINUE
           ENDIF
        ENDIF
      ELSE  !WRITING OUT FORCES TO RESIDUAL FILE
!
! *** CONSTRACT GLOBAL MACRO VECTORS
!
          CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,UVEL_PP,VVEL_PP,&
     &         PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,&
     &           NGRPS,RANK)

        IF(RANK.EQ.0)THEN
          IF((ITIME.EQ.1000).OR.(ITIME.EQ.2000).OR.(ITIME.EQ.3000).OR.(ITIME.EQ.4000)&
     &  .OR.(ITIME.EQ.5000).OR.(ITIME.EQ.6000).OR.(ITIME.EQ.7000)&
     &  .OR.(ITIME.EQ.8000).OR.(ITIME.EQ.9000).OR.(ITIME.EQ.10000)&
     &  .OR.(ITIME.EQ.11000).OR.(ITIME.EQ.12000).OR.(ITIME.EQ.13000)&
     &  .OR.(ITIME.EQ.14000).OR.(ITIME.EQ.15000).OR.(ITIME.EQ.20000)&
     &  .OR.(ITIME.EQ.25000).OR.(ITIME.EQ.30000).OR.(ITIME.EQ.35000)&
     &  .OR.(ITIME.EQ.40000).OR.(ITIME.EQ.45000).OR.(ITIME.EQ.50000))THEN
             PRINT *,'WRITING2',RANK,ITIME
             FLUSH(6)
             WRITE(17,500) ITIME
             WRITE(18,500) ITIME
             DO 3011 IP=1,NPOIN
               WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)
               WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP)
 3011        CONTINUE
          ENDIF
        ENDIF
      ENDIF
! 
! *** WRITE TO MAX RESIDUAL FILE 
! 
     IF(RANK.EQ.0)THEN
        IF(IVD%FORCEOUT.EQ.0)THEN     !IF NOT OUTPUTTING FORCES
          SUMRES=MAXVAL(RESIDUAL) 
          if(ITIME.EQ.1)SRITM1=SUMRES
          WRITE(16,*) ITIME,LOG(SUMRES/SRITM1)
        ELSE             !IF OUTPUTTING FORCES IN RESIDUAL FILE
          SUMRES=MAXVAL(RESIDUAL)
          if(ITIME.EQ.1)SRITM1=SUMRES
          CALL GETFOR(LIFT,DRAG,NBOUN,NPOIN,NBNOI,PS,COORD,BSIDO)
          WRITE(16,*) ITIME,LOG(SUMRES/SRITM1),LIFT,DRAG
        ENDIF
     ENDIF 
! 
! *** SYNCHRONISE 
! 
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR) 
! 
! *** END IF TIMESTEPPING LOOP 
! 
30000 CONTINUE 
!
! *** CONSTRACT GLOBAL MACRO VECTORS
!
 !        CALL GETMAC(maxNPOIN_PP,NPOIN_PP,IPCOM_PP,ND_PP,RHO_PP,UVEL_PP,VVEL_PP,& 
 !    &         PS_PP,TEMP_PP,NPOIN,ND,RHO,UVEL,VVEL,PS,TEMP,& 
 !    &           NGRPS,RANK)
! 
! *** WRITE OUTPUT DATA TO RESULTS FILE 
!  
       IF(RANK.EQ.0)THEN
! *** CONSTRUCT GLOBAL MACRO VARIABLE VECTORS 
      WRITE(17,500) ITIME 
      WRITE(18,500) ITIME 
      DO 3002 IP=1,NPOIN 
      WRITE(17,600)IP,ND(IP),UVEL(IP),VVEL(IP)   
      WRITE(18,600)IP,RHO(IP),PS(IP),TEMP(IP) 
 3002 CONTINUE
      ENDIF
! 
! *** CLOSE CHANNELS 16 AND 17 AND 18 
! 
      IF(RANK.EQ.0)THEN 
      CLOSE(16) 
      CLOSE(17) 
      CLOSE(18) 
      ENDIF 
!      IF(RANK.EQ.3)THEN
!      CLOSE(15)
!      CLOSE(25)
!      CLOSE(35)
!      CLOSE(45)
!      CLOSE(55)
!      CLOSE(65)
!      CLOSE(75)
!      CLOSE(85)
!      ENDIF

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
