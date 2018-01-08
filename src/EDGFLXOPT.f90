           SUBROUTINE EDGFLXOPT(NELEM,NSIDE_PP,ISIDE,RHS,& 
     &             NX,NY,EL,UX,UY,RSIDO,BSIDO,NBOUN,& 
     &                       NBNOR,ALPHA,ETA,VNPNT,& 
     &               IV,DISNF,UMEAN,CINF,rv,LCOMM_PP,& 
     &             NGRPS,ISCOM_PP,MPI_RANK_P,NCOMM_PP,VCORD,RORDER,&
     &              TORDER,SDCOM_PP,R,M,ITIME,GCOMM,MPI_COMM_P,&
     &                     VSPACE_FIRST,VSPACE_LAST)
! 
! *** SUBROUTINE TO CALCULATE THE FLUXES TRANSFERRED BETWEEN ELEMENTS AT EDGES 
! 
      IMPLICIT NONE 
      INCLUDE 'mpif.h' 
! 
      INTEGER NELEM,NSIDE_PP,NBOUN,NBNOR,MXCOM_PP,NGRPS,ISG,FLAG 
      INTEGER IS_PP,ISC_PP,IG,TEST,IGT,NSIDETMP,IST,ISGT,NCOMM_PP,ISC 
      INTEGER LCOMM_PP(NSIDE_PP),ISCOM_PP(NSIDE_PP),MPI_IERR,MPI_RANK_P
      INTEGER MPI_COMM_P,IE,ISL,ISR 
      INTEGER ISIDE(8,NSIDE_PP),VNPNT,IV,OPPSLAVERANK,RORDER,TORDER 
      INTEGER MPI_STATUS(MPI_STATUS_SIZE),BSIDO(6,NBOUN),INL1,INL2 
      INTEGER IS,IP1,IEL,IER,IN,JN,INR1, INR2,SDCOM_PP(3,NCOMM_PP) 
      INTEGER ITIME, IS_OTHER, IELR, INLR1,INLR2
      INTEGER VSPACE_FIRST,VSPACE_LAST
! 
      REAL NX(NSIDE_PP),NY(NSIDE_PP),UX,UY,ALEN,RSIDO(NBNOR,NBOUN) 
      REAL MMAT(2,2),FN(2),EL(NSIDE_PP),CM,ALPHA,CINF(4) 
      REAL RHS(1,3,NELEM),RHSI(2)
      REAL DISNF(3,VSPACE_FIRST:VSPACE_LAST,NELEM) 
      REAL UMEAN(NELEM),EDGUN(2),SCPR,RV,VCORD(3,VNPNT) 
      REAL FLUXN(4),FLUYN(4),ETA(NBOUN),CO,R,M
      INTEGER GCOMM(NGRPS,NGRPS)
 
!     COMMUNICATION STUFF (OPTIMIZATION)
      INTEGER :: SENDRECV_TOT_LENGTH_MAX, OFFSET 
      INTEGER SENDRECV_MAX_LENGTHS(NGRPS)! MAX LENGTHS OF DATA TO SEND/RECV TO/FROM OTHER RANKS
      INTEGER SEND_LENGTHS(NGRPS) ! ACTUAL LENGTHS (SEND) FOR EACH MPI_RANK_P
      INTEGER RECV_LENGTHS(NGRPS) ! ACTUAL LENGTHS (RECV) FOR EACH MPI_RANK_P
      INTEGER SENDRECV_START_OFFSETS(NGRPS) ! OFFSET OF FIRST ELEMENT TO SEND TO ANY GIVEN MPI_RANK_P 
      REAL :: SEND_EDGE_DATA(5000)!
      REAL :: RECV_EDGE_DATA(5000)!
      INTEGER :: SEND_EDGE_DATA_IDX(5000)! LOCAL EDGE INDICES ON OPPOSITE MPI_RANK_P    
      INTEGER :: RECV_EDGE_DATA_IDX(5000)! LOCAL EDGE INDICES ON CURRENT MPI_RANK_P
      INTEGER :: MPI_COMM_SLAVES ! SLAVES COMMUNICATOR
      INTEGER :: COLOR,SLAVERANK 
      INTEGER :: EDGCOUNT

      DATA MMAT/2.0,1.0,1.0,2.0/ 
      FLAG=0 
      EDGCOUNT = 0

      CO=0.0



!     SETUP COMMUNICATION STUFF (OPTIMIZATION)

      IF (MPI_RANK_P.NE.0) THEN ! adadadfadfad
     
        DO IG=1,NGRPS
            SENDRECV_MAX_LENGTHS(IG) = GCOMM(MPI_RANK_P,IG)
        ENDDO

        SENDRECV_START_OFFSETS(1) = 0
        DO IG=2,NGRPS
            SENDRECV_START_OFFSETS(IG) = SENDRECV_START_OFFSETS(IG-1)+ &
     &                                     SENDRECV_MAX_LENGTHS(IG-1)
        ENDDO
        SENDRECV_TOT_LENGTH_MAX = SENDRECV_START_OFFSETS(NGRPS) + &
     &                           SENDRECV_MAX_LENGTHS(NGRPS)

        DO IG=1,NGRPS
            SEND_LENGTHS(IG) = 0
        ENDDO
        
        DO OFFSET=1,5000                   !DEBUG
          SEND_EDGE_DATA_IDX(OFFSET) = 0!DEBUG
          RECV_EDGE_DATA_IDX(OFFSET) = 0!DEBUG
        ENDDO!DEBUG

   
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
! 
! *** BEGIN LOOP OVER EACH ELEMENT EDGE 
!
         DO 1000 IS=1,NSIDE_PP 
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
           SCPR=UX*NX(IS)+UY*NY(IS) 
! 
! *** SKIP TO SEPARATE SUBROUTINE FOR BOUNDARY EDGES 
!	 
           IF(IER.EQ.0)THEN  ! dlkfsdaklfasklda 
             CALL GETBOU(NBOUN,NBNOR,NELEM,BSIDO,IEL,RSIDO,IP1,& 
     &            RHS,INL1,INL2,UX,UY,ALPHA,ETA,VNPNT,&
     &            IV,DISNF,UMEAN,CINF,rv,MPI_RANK_P,VCORD,RORDER,TORDER&
     &            ,R,M,VSPACE_FIRST,VSPACE_LAST) 

           ELSEIF((IER.NE.-1).AND.(IEL.NE.-1))THEN !dlkfsdaklfasklda 
! *** FOR INTERNAL SIDES: 
!  
! *** STORE THE UPSTREAM EDGE UNKNOWNS IN UNK1 AND UNK2
! *** AND COMPUTE THE BGK COLLISION TERM 
! 
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

           ELSE   !dlkfsdaklfasklda 


! *** NOW DEAL WITH PROCESSOR BOUNDARY EDGES HERE!!!!! 
               EDGCOUNT = EDGCOUNT + 1 ! index for edges on partition border
               
               IS_PP = SDCOM_PP(1,EDGCOUNT)
               IF(IS_PP.NE.IS) THEN  ! DEBUG
                 WRITE(*,*) "Rank, IS_PP != IS", MPI_RANK_P,IS_PP,IS
               STOP
               ENDIF                 !DEBUG


               OPPSLAVERANK = SDCOM_PP(3,EDGCOUNT) - 1 ! ranks in the 
                                                 ! MPI_COMM_SLAVES communicator
                                                 ! are the rank in
                                                 ! MPI_COMM_P communicator
                                                 ! decreased by 1
               IS_OTHER = SDCOM_PP(2,EDGCOUNT)
              
           IF((IEL.EQ.-1).AND.(SCPR.LT.0.0))THEN 
               INLR1 = INR1
               INLR2 = INR2
               IELR = IER
               FLAG = -1
           ELSEIF((IER.EQ.-1).AND.(SCPR.GT.0.0))THEN
               INLR1 = INL1
               INLR2 = INL2
               IELR = IEL
               FLAG = 1
           ELSE
           CYCLE ! Sometimes SCPR == 0, we skip this iteration.
           ENDIF
               EDGUN(1)=DISNF(INLR1,IV,IELR)-UMEAN(IELR)
               EDGUN(2)=DISNF(INLR2,IV,IELR)-UMEAN(IELR)
!
! *** UPDATE RHS
!
               CALL GETFLU(2,EDGUN,FLUXN,FLUYN,UX,UY)
!
               ALEN=(1.0/6.0)*EL(IS_PP)
!
!  *** MULTPLY UP(DOWN)STREAM FLUXES BY THE APPROPRIATE NORMALS
!
               FN(1)=NX(IS_PP)*FLUXN(1)+NY(IS_PP)*FLUYN(1)
               FN(2)=NX(IS_PP)*FLUXN(2)+NY(IS_PP)*FLUYN(2)
!
               CALL RFILLV(RHSI,2,CO)
!
               DO IN=1,2
                 DO JN=1,2
                   CM=FN(JN)*MMAT(IN,JN)*ALEN
                   RHSI(IN)=RHSI(IN)-CM
                 ENDDO
               ENDDO
               
!
! *** UPDATE THE DOWN(UP)STREAM ELEMENT RHS VALUES
!
               RHS(1,INLR1,IELR)=RHS(1,INLR1,IELR)+FLAG*RHSI(1)
               RHS(1,INLR2,IELR)=RHS(1,INLR2,IELR)+FLAG*RHSI(2)
                  


               OFFSET=SENDRECV_START_OFFSETS(OPPSLAVERANK+1)+&
     &                                SEND_LENGTHS(OPPSLAVERANK+1)
               SEND_EDGE_DATA(2*OFFSET+1)=RHSI(1)
               SEND_EDGE_DATA(2*OFFSET+2)=RHSI(2)
               SEND_EDGE_DATA_IDX(OFFSET+1)=FLAG*IS_OTHER ! STORING FLAG IN THIS WAY

               SEND_LENGTHS(OPPSLAVERANK+1)= &
     &                     SEND_LENGTHS(OPPSLAVERANK+1)+1 ! will be
                                                        ! multiplied by 2 when necessary

!               WRITE(50+MPI_RANK_P,"(A9,3I4.1,I5.1,2ES13.4E2)"),"NEWSEND",&!DEBUG
!     &                       MPI_RANK_P,MPI_RANK_P,OPPSLAVERANK+1,IS_OTHER,RHSI   !DEBUG !+1 to compare with old

           ENDIF !dlkfsdaklfasklda 

 1000 CONTINUE ! END OF CYCLE OVER SIDES
     
      CALL MPI_ALLTOALL(SEND_LENGTHS,1,MPI_INTEGER,&
     &             RECV_LENGTHS,1,MPI_INTEGER,&
     &             MPI_COMM_SLAVES,MPI_IERR)
!      CALL MPI_BARRIER(MPI_COMM_SLAVES,MPI_IERR)! DEBUG
!      WRITE(*,"(A10,8I5.1)") "NEWSEND", MPI_RANK_P, SEND_LENGTHS !DEBUG
!      WRITE(*,"(A10,8I5.1)") "NEWRECV", MPI_RANK_P, RECV_LENGTHS !DEBUG
!      CALL MPI_BARRIER(MPI_COMM_SLAVES,MPI_IERR) !DEBUG
!
      CALL MPI_ALLTOALLV(SEND_EDGE_DATA,2*SEND_LENGTHS,&
     &      2*SENDRECV_START_OFFSETS,MPI_REAL,&
     &      RECV_EDGE_DATA,2*RECV_LENGTHS,&
     &      2*SENDRECV_START_OFFSETS,MPI_REAL,&
     &      MPI_COMM_SLAVES,MPI_IERR)

      CALL MPI_ALLTOALLV(SEND_EDGE_DATA_IDX,SEND_LENGTHS,&
     &      SENDRECV_START_OFFSETS,MPI_INTEGER,&
     &      RECV_EDGE_DATA_IDX,RECV_LENGTHS,&
     &      SENDRECV_START_OFFSETS,MPI_INTEGER,&
     &      MPI_COMM_SLAVES,MPI_IERR)

      DO IG=1,NGRPS
         DO OFFSET=SENDRECV_START_OFFSETS(IG),&
     &              SENDRECV_START_OFFSETS(IG)+RECV_LENGTHS(IG)-1
            IS = ABS(RECV_EDGE_DATA_IDX(OFFSET+1))
            IF(RECV_EDGE_DATA_IDX(OFFSET+1).LT.0) THEN 
               ! WE ARE ON THE LEFT OF THE EDGE
               IELR =  ISIDE(3,IS)
               INLR1 = ISIDE(5,IS) 
               INLR2 = ISIDE(6,IS) 
               FLAG = -1
            ELSEIF(RECV_EDGE_DATA_IDX(OFFSET+1).GT.0) THEN 
               ! WE ARE ON THE RIGHT OF THE EDGE
               IELR =  ISIDE(4,IS)
               INLR1 = ISIDE(7,IS) 
               INLR2 = ISIDE(8,IS) 
               FLAG = 1
            ELSE
             WRITE(*,*) "THIS SHOULD NOT HAPPEN."
            ENDIF 
            
            RHSI(1) = RECV_EDGE_DATA(2*OFFSET+1)
            RHSI(2) = RECV_EDGE_DATA(2*OFFSET+2)
             
            RHS(1,INLR1,IELR)=RHS(1,INLR1,IELR)-FLAG*RHSI(1)
            RHS(1,INLR2,IELR)=RHS(1,INLR2,IELR)-FLAG*RHSI(2)
!            WRITE(50+MPI_RANK_P,"(A9,3I4.1,4I5.1,2ES13.4E2)"),"NEWRECV",&!DEBUG
!     &         MPI_RANK_P,IG,MPI_RANK_P,INLR1,INLR2,IELR,IS,RHSI   !DEBUG !+1 to compare with old


         ENDDO
      ENDDO

      ENDIF  ! IF(MPI_RANK_P.NE.0)  !kdjfhsewew
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
