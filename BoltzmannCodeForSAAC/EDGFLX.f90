           SUBROUTINE EDGFLX(NELEM,NSIDE,ISIDE,RHS,& 
     &             NX,NY,EL,UX,UY,RSIDO,BSIDO,NBOUN,& 
     &                       NBNOR,ALPHA,ETA,VNPNT,& 
     &               IV,DISNF,UMEAN,CINF,rv,LCOMM_PP,& 
     &             NGRPS,ISCOM_PP,RANK,NCOMM_PP,VCORD,RORDER,&
     &              TORDER,SDCOM_PP,R,M,ITIME)
! 
! *** SUBROUTINE TO CALCULATE THE FLUXES TRANSFERRED BETWEEN ELEMENTS AT EDGES 
! 
      IMPLICIT NONE 
      INCLUDE 'mpif.h' 
! 
      INTEGER NELEM,NSIDE,NBOUN,NBNOR,MXCOM_PP,NGRPS,ISG,FLAG 
      INTEGER IS_PP,ISC_PP,IG,TEST,IGT,NSIDETMP,IST,ISGT,NCOMM_PP,ISC 
      INTEGER LCOMM_PP(NSIDE),ISCOM_PP(NSIDE),MPI_IERR,RANK,IE,ISL,ISR 
      INTEGER ISIDE(8,NSIDE),VNPNT,IV,RANKL,RANKR,RORDER,TORDER 
      INTEGER MPI_STATUS(MPI_STATUS_SIZE),BSIDO(6,NBOUN),INL1,INL2 
      INTEGER IS,IP1,IEL,IER,IN,JN,INR1, INR2,SDCOM_PP(3,NCOMM_PP) 
      INTEGER ITIME
! 
      REAL NX(NSIDE),NY(NSIDE),UX,UY,ALEN,RSIDO(NBNOR,NBOUN) 
      REAL MMAT(2,2),FN(2),EL(NSIDE),CM,ALPHA,CINF(4) 
      REAL RHS(1,3,NELEM),RHSI(2),DISNF(3,VNPNT,NELEM) 
      REAL UMEAN(NELEM),EDGUN(2),SCPR,RV,VCORD(3,VNPNT) 
      REAL FLUXN(4),FLUYN(4),ETA(NBOUN),CO,R,M
! 
      DATA MMAT/2.0,1.0,1.0,2.0/ 
      FLAG=0 
      CO=0.0
! 
! *** BEGIN LOOP OVER EACH ELEMENT EDGE 
!
       IF(RANK.NE.0)THEN 
         DO 1000 IS=1,NSIDE 
! 
! *** NODES AND ELEMENTS ASSOCIATED WITH EDGE 
! 
           IP1=ISIDE(1,IS)                 !FIRST NODE NUMBER 
           IEL=ISIDE(3,IS)                 !LHS ELEMENT NUMBER 
           IER=ISIDE(4,IS)                 !RHS ELEMENT NUMBER 
           INL1=ISIDE(5,IS)        !LHS LOCAL NODE 1 
           INL2=ISIDE(6,IS)        !LHS LOCAL NODE 2 
           INR1=ISIDE(7,IS)        !RHS LOCAL NODE 1 
           INR2=ISIDE(8,IS)        !RHS LOCAL NODE 2
! 
! *** SKIP TO SEPARATE SUBROUTINE FOR BOUNDARY EDGES 
!	 
           IF(IER.EQ.0)THEN 
             CALL GETBOU(NBOUN,NBNOR,NELEM,BSIDO,IEL,RSIDO,IP1,& 
     &            RHS,INL1,INL2,UX,UY,ALPHA,ETA,VNPNT,&
     &            IV,DISNF,UMEAN,CINF,rv,RANK,VCORD,RORDER,TORDER&
     &            ,R,M) 
             GOTO 1001 
           ENDIF
! ------------------------------------------------------------------------------------------------- 
!  
! *** FOR INTERNAL SIDES: 
!
           IF((IER.NE.-1).AND.(IEL.NE.-1))THEN 
!  
! *** STORE THE UPSTREAM EDGE UNKNOWNS IN UNK1 AND UNK2
! *** AND COMPUTE THE BGK COLLISION TERM 
! 
             SCPR=UX*NX(IS)+UY*NY(IS) 
             IF(SCPR.GT.0.0)THEN 
               EDGUN(1)=DISNF(INL1,IV,IEL)-UMEAN(IEL) 
               EDGUN(2)=DISNF(INL2,IV,IEL)-UMEAN(IEL)
             ELSE 
               EDGUN(1)=DISNF(INR1,IV,IER)-UMEAN(IER) 
               EDGUN(2)=DISNF(INR2,IV,IER)-UMEAN(IER)
             ENDIF 
! 
! *** GET THE FLUXES 
! 
             CALL GETFLU(2,EDGUN,FLUXN,FLUYN,UX,UY) 
! 
! *** SIDE LENGTH x (1/6) 
             ALEN=(1.0/6.0)*EL(IS) 
!	 
! *** MULTIPLY THE !UPSTREAM! FLUXES BY THE APPROPRIATE NORMALS 
! 
             FN(1)=NX(IS)*FLUXN(1)+NY(IS)*FLUYN(1) 
             FN(2)=NX(IS)*FLUXN(2)+NY(IS)*FLUYN(2)
! 
        CALL RFILLV(RHSI,2,CO) 
! 
             DO 401 IN=1,2 
             DO 402 JN=1,2 
               CM=FN(JN)*MMAT(IN,JN)*ALEN 
               RHSI(IN)=RHSI(IN)-CM 
  402 CONTINUE 
  401 CONTINUE 
! 
! *** UPDATE THE UPSTREAM AND DOWNSTREAM ELEMENT RHS VALUES 
! 
             RHS(1,INR1,IER)=RHS(1,INR1,IER)-RHSI(1) 
             RHS(1,INR2,IER)=RHS(1,INR2,IER)-RHSI(2) 
             RHS(1,INL1,IEL)=RHS(1,INL1,IEL)+RHSI(1) 
             RHS(1,INL2,IEL)=RHS(1,INL2,IEL)+RHSI(2)  
! 
           ENDIF
 1001 CONTINUE 
!        
! 
!	   
 1000 CONTINUE
       ENDIF 
! 
! ----------------------------------------------------------------------------- 
! 
! 
! *** NOW DEAL WITH PROCESSOR BOUNDARY EDGES HERE!!!!! 
! 
! *** LOOP OVER THE GROUPS 
!
       DO 4000 IG=1,NGRPS 
         IF(RANK.EQ.IG)THEN 
! 
! *** SET UP A TEMPORARY NSIDE VARIABLE THAT ALL PROCS MUST LOOP OVER 
! *** AND COMMUNICATE TO ALL PROCS
! 
           NSIDETMP=NCOMM_PP
         ENDIF
         CALL MPI_BCAST(NSIDETMP,1,MPI_INTEGER,IG,MPI_COMM_WORLD,MPI_IERR)
! 
! ***  ALL PROCS LOOP OVER PROC IG'S COMMUNICATING EDGES 
! 
         DO 4001 ISC_PP=1,NSIDETMP
           IF(RANK.EQ.IG)THEN
             IS_PP=SDCOM_PP(1,ISC_PP)
             IEL=ISIDE(3,IS_PP) 
             IER=ISIDE(4,IS_PP) 
             SCPR=UX*NX(IS_PP)+UY*NY(IS_PP)
           ENDIF
           CALL MPI_BCAST(IEL,1,MPI_INTEGER,IG,MPI_COMM_WORLD,MPI_IERR)
           CALL MPI_BCAST(IER,1,MPI_INTEGER,IG,MPI_COMM_WORLD,MPI_IERR)
           CALL MPI_BCAST(SCPR,1,MPI_REAL,IG,MPI_COMM_WORLD,MPI_IERR)
           CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
           IF((IEL.EQ.-1).AND.(SCPR.LT.0.0))THEN 
             IF(RANK.EQ.IG)THEN
               RANKL=SDCOM_PP(3,ISC_PP) 
               ISL=SDCOM_PP(2,ISC_PP)
!
! *** SEND THE OPPOSITE PROCESSOR THE DETAILS REQUIRED TO CONSTRUCT THE INTER-
! *** ELEMENT FLUX
!
               INR1=ISIDE(7,IS_PP)
               INR2=ISIDE(8,IS_PP)
               EDGUN(1)=DISNF(INR1,IV,IER)-UMEAN(IER)
               EDGUN(2)=DISNF(INR2,IV,IER)-UMEAN(IER)
!
! *** UPDATE RHS
!
               CALL GETFLU(2,EDGUN,FLUXN,FLUYN,UX,UY)
!
               ALEN=(1.0/6.0)*EL(IS_PP)
!
!  *** MULTPLY UPSTREAM FLUXES BY THE APPROPRIATE NORMALS
!
               FN(1)=NX(IS_PP)*FLUXN(1)+NY(IS_PP)*FLUYN(1)
               FN(2)=NX(IS_PP)*FLUXN(2)+NY(IS_PP)*FLUYN(2)
!
               CALL RFILLV(RHSI,2,CO)
!
               DO 501 IN=1,2
                 DO 502 JN=1,2
                   CM=FN(JN)*MMAT(IN,JN)*ALEN
                   RHSI(IN)=RHSI(IN)-CM
 502 CONTINUE
 501 CONTINUE
!
! *** UPDATE THE DOWNSTREAM ELEMENT RHS VALUES
!
               RHS(1,INR1,IER)=RHS(1,INR1,IER)-RHSI(1)
               RHS(1,INR2,IER)=RHS(1,INR2,IER)-RHSI(2)
             ENDIF
             CALL MPI_BCAST(RANKL,1,MPI_INTEGER,IG,MPI_COMM_WORLD,MPI_IERR)
             IF(RANK.EQ.IG)THEN
               CALL MPI_SEND(ISL,1,MPI_INTEGER,RANKL,1,&
     &                         MPI_COMM_WORLD,MPI_IERR)
               CALL MPI_SEND(RHSI,2,MPI_REAL,RANKL,2,&
                               MPI_COMM_WORLD,MPI_IERR)
             ENDIF
             IF(RANK.EQ.RANKL)THEN
               CALL MPI_RECV(ISL,1,MPI_INTEGER,IG,1,&
     &                         MPI_COMM_WORLD,MPI_STATUS,MPI_IERR)
               CALL MPI_RECV(RHSI,2,MPI_REAL,IG,2,&
     &                         MPI_COMM_WORLD,MPI_STATUS,MPI_IERR)
!
! *** GET THE RELEVENT DISCONTINUOUS NODE NUMBERS AND ELEMENT NUMBER
!
               IEL=ISIDE(3,ISL)
               INL1=ISIDE(5,ISL)
               INL2=ISIDE(6,ISL)
!
! *** UPDATE THE DOWNSTREAM ELEMENT RHS VALUES
!
               RHS(1,INL1,IEL)=RHS(1,INL1,IEL)+RHSI(1)
               RHS(1,INL2,IEL)=RHS(1,INL2,IEL)+RHSI(2)
             ENDIF
             CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
           ELSEIF((IER.EQ.-1).AND.(SCPR.GT.0.0))THEN
             IF(RANK.EQ.IG)THEN
               RANKR=SDCOM_PP(3,ISC_PP)
               ISR=SDCOM_PP(2,ISC_PP)
 
!
! *** SEND THE OPPOSITE PROCESSOR THE DETAILS REQUIRED TO CONSTRUCT THE INTER-
! *** ELEMENT FLUX
!
               INL1=ISIDE(5,IS_PP)
               INL2=ISIDE(6,IS_PP)
               EDGUN(1)=DISNF(INL1,IV,IEL)-UMEAN(IEL)
               EDGUN(2)=DISNF(INL2,IV,IEL)-UMEAN(IEL)
!
! *** UPDATE RHS
!
               CALL GETFLU(2,EDGUN,FLUXN,FLUYN,UX,UY)
!
               ALEN=(1.0/6.0)*EL(IS_PP)
!
!  *** MULTPLY UPSTREAM FLUXES BY THE APPROPRIATE NORMALS
!
               FN(1)=NX(IS_PP)*FLUXN(1)+NY(IS_PP)*FLUYN(1)
               FN(2)=NX(IS_PP)*FLUXN(2)+NY(IS_PP)*FLUYN(2)
!
               CALL RFILLV(RHSI,2,CO)
!
               DO 701 IN=1,2
                 DO 702 JN=1,2
                   CM=FN(JN)*MMAT(IN,JN)*ALEN
                   RHSI(IN)=RHSI(IN)-CM
 702 CONTINUE
 701 CONTINUE
!
! *** UPDATE THE UPSTREAM ELEMENT RHS VALUES
!
               RHS(1,INL1,IEL)=RHS(1,INL1,IEL)+RHSI(1)
               RHS(1,INL2,IEL)=RHS(1,INL2,IEL)+RHSI(2)
             ENDIF
             CALL MPI_BCAST(RANKR,1,MPI_INTEGER,IG,MPI_COMM_WORLD,MPI_IERR)
             IF(RANK.EQ.IG)THEN
               CALL MPI_SEND(ISR,1,MPI_INTEGER,RANKR,1,&
     &                     MPI_COMM_WORLD,MPI_IERR)
               CALL MPI_SEND(RHSI,2,MPI_REAL,RANKR,2,&
                            MPI_COMM_WORLD,MPI_IERR)
             ENDIF
             IF(RANK.EQ.RANKR)THEN
               CALL MPI_RECV(ISR,1,MPI_INTEGER,IG,1,&
     &                      MPI_COMM_WORLD,MPI_STATUS,MPI_IERR)
               CALL MPI_RECV(RHSI,2,MPI_REAL,IG,2,&
     &                         MPI_COMM_WORLD,MPI_STATUS,MPI_IERR)
!
! *** GET THE RELEVENT DISCONTINUOUS NODE NUMBERS AND ELEMENT NUMBER
!
               IER=ISIDE(4,ISR)
               INR1=ISIDE(7,ISR)
               INR2=ISIDE(8,ISR)
!
! *** UPDATE THE DOWNSTREAM ELEMENT RHS VALUES
!
               RHS(1,INR1,IER)=RHS(1,INR1,IER)-RHSI(1)
               RHS(1,INR2,IER)=RHS(1,INR2,IER)-RHSI(2)
             ENDIF
             CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
           ENDIF
!
 4001 CONTINUE
!
 4000 CONTINUE
! 
! *** FORMAT STATEMENTS 
! 
 3000 FORMAT('ERROR IN LCOMM_PP MATRIX!!!') 
 3001 FORMAT('INCREASE MXCOM PARAMETER IN EDGFLX SUBROUTINE!!!') 
! 
      RETURN 
! 
 2000 PRINT*,'ERROR IN DATA CONTAINED IN ISIDE!! PROGRAM STOPPED'     
      STOP 
      END 
