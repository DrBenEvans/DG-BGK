      SUBROUTINE MACROS(rv,VNPNT,SUMWEIGHT,NPOIN_PP,NNODE,& 
     &           DISNF_PP,VCORD,INTMA_PP,NELEM_PP,GEOME_PP,NGEOM,& 
     &           ND,RHO_PP,UVEL_PP,VVEL_PP,PS_PP,TEMP_PP,&
     &           DISND_PP,DISUX_PP,DISUY_PP,DISPS,M,R,&
     &           MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
     &           VSPACE_FIRST,VSPACE_LAST)
! 
! *** SUBROUTINE CALCULATES MACROSCOPIC PROPERTIES FROM THE DISTRIBUTION 
! *** FUNCTION 
! 
      IMPLICIT NONE 
      INCLUDE 'mpif.h'
! 
      INTEGER VNPNT,NPOIN_PP,NNODE,NELEM_PP,NGEOM 
      INTEGER IN,IV,IP,IE,MPI_IERR 
      INTEGER INTMA_PP(NNODE,NELEM_PP),TMP 
      REAL rv,R,AE,RT,PI,JAC,ETA,ZETA,N 
      REAL WEIGHT,THETA,BETA,M,MM 
      REAL ND(NPOIN_PP),VCORD(3,VNPNT) 
      REAL RHO_PP(NPOIN_PP),UVEL_PP(NPOIN_PP),VVEL_PP(NPOIN_PP)
      REAL TEMP_PP(NPOIN_PP),PS_PP(NPOIN_PP),FLAG(NPOIN_PP)
      REAL CX,CY,CXDASH,CYDASH,SUMWEIGHT,C4,FRAC
      REAL GEOME_PP(NGEOM,NELEM_PP)
      REAL DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,NELEM_PP)
      REAL DISND_PP(NNODE,NELEM_PP)
      REAL DISUX_PP(NNODE,NELEM_PP)
      REAL DISUY_PP(NNODE,NELEM_PP)

      REAL DISPS(NNODE,NELEM_PP)

      INTEGER MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,GROUP_V
      INTEGER VNPNT_PART,VSPACE_FIRST,VSPACE_LAST

      REAL  INT_ARR(3*NELEM_PP), INTG_ARR(3*NELEM_PP)
      REAL INTU_ARR(3*NELEM_PP),INTUG_ARR(3*NELEM_PP)
      REAL INTV_ARR(3*NELEM_PP),INTVG_ARR(3*NELEM_PP)
      REAL INTP_ARR(3*NELEM_PP),INTPG_ARR(3*NELEM_PP)
      INTEGER IDX


! 
! *** MOLECULAR MASS OF O2 AND GAS CONSTANT 
! 
      PARAMETER(PI=3.1416,C4=4.0) 
!
! *** COMPUTE THE WEIGHTING FACTOR FOR THE FULL V-SPACE QUADRATURE
!
      FRAC=C4/SUMWEIGHT
!
! *** COMPUTE THE MOLECULAR MASS
!
      MM = M/6.022E26
! 
! *** INITIALISE ND(NPOIN_PP),RHO_PP(NPOIN_PP),UVEL_PP(NPOIN_PP),VVEL_PP(NPOIN_PP),TEMP_PP(NPOINT),PS_PP(NPOIN_PP) 
! 
      CALL RFILLV(RHO_PP,NPOIN_PP,0.0) 
      CALL RFILLV(UVEL_PP,NPOIN_PP,0.0)
      CALL RFILLV(VVEL_PP,NPOIN_PP,0.0) 
      CALL RFILLV(TEMP_PP,NPOIN_PP,0.0) 
      CALL RFILLV(PS_PP,NPOIN_PP,0.0) 
      CALL RFILLV(ND,NPOIN_PP,0.0) 
      CALL RFILLV(FLAG,NPOIN_PP,0.0)
!
! *** INITIALISE THE DISCONTINUOUS BULK VARIABLE ARRAYS
!
      CALL RFILLM(DISND_PP,NNODE,NELEM_PP,0.0)
      CALL RFILLM(DISUX_PP,NNODE,NELEM_PP,0.0)
      CALL RFILLM(DISUY_PP,NNODE,NELEM_PP,0.0)
      CALL RFILLM(DISPS,NNODE,NELEM_PP,0.0)
! 
! *** LOOP OVER EACH DISCONTINUOUS NODE IN PHYSICAL SPACE 
! 
      DO IE=1,NELEM_PP !fdgjshfvsfkskjfvs 
        AE=GEOME_PP(7,IE) 
        DO IN=1,NNODE  !fjhsdfkgjhccccfflkjdfgs
          IP=INTMA_PP(IN,IE) 
          FLAG(IP)=FLAG(IP)+AE 
          IDX = (IE-1)*3+IN
! 
! *** LOOP OVER ALL VSPACE NODES TO CALCULATE NUMBER DENSITY & BULK VELOCITIES 
! *** LOOP OVER VSPACE NODES TO CALCULATE STATIC PRESSURE 
!	
          INT_ARR (IDX)=0.0 
          INTU_ARR(IDX)=0.0 
          INTV_ARR(IDX)=0.0 
          INTP_ARR(IDX)=0.0 
          DO IV=VSPACE_FIRST,VSPACE_LAST!dfjhasdfa
            ETA=VCORD(1,IV) 
            ZETA=VCORD(2,IV) 
            WEIGHT=VCORD(3,IV)
            RT=ETA*(rv/2)+(rv/2)!MAP BACK TO POLAR COORDINATES 
            THETA=ZETA*PI
            CX=RT*COS(THETA)!CONVERT TO CARTESIANS 
            CY=RT*SIN(THETA) 
            JAC=PI*rv*rv*0.25*(ETA+1)!CALCULATE THE JACOBIAN OF THE TRANSFORMATION 
            CXDASH=CX-DISUX_PP(IN,IE)
            CYDASH=CY-DISUY_PP(IN,IE) 
            TMP=0.5*(CXDASH*CXDASH+CYDASH*CYDASH) 
             INT_ARR(IDX)= INT_ARR(IDX)+&
     &                      DISNF_PP(IN,IV,IE)*WEIGHT*JAC*FRAC
            INTU_ARR(IDX)=INTU_ARR(IDX)+&
     &                      CX*DISNF_PP(IN,IV,IE)*WEIGHT*JAC*FRAC
            INTV_ARR(IDX)=INTV_ARR(IDX)+&
     &                      CY*DISNF_PP(IN,IV,IE)*WEIGHT*JAC*FRAC
            INTP_ARR(IDX)=INTP_ARR(IDX)+&
     &                      TMP*DISNF_PP(IN,IV,IE)*WEIGHT*JAC*FRAC
          ENDDO ! DO IV=VSPACE_FIRST,VSPACE_LAST!dfjhasdfa
        ENDDO !  DO IN=1,NNODE  !fjhsdfkgjhccccfflkjdfgs
      ENDDO !  DO IE=1,NELEM_PP !fdgjshfvsfkskjfvs 

      MPI_ALLREDUCE( INT_ARR, INTG_ARR,3*NELEM_PP,MPI_REAL,MPI_SUM,&
     &               MPI_COMM_V,MPI_IERR)
      MPI_ALLREDUCE(INTU_ARR,INTUG_ARR,3*NELEM_PP,MPI_REAL,MPI_SUM,&
     &               MPI_COMM_V,MPI_IERR)
      MPI_ALLREDUCE(INTV_ARR,INTVG_ARR,3*NELEM_PP,MPI_REAL,MPI_SUM,&
     &               MPI_COMM_V,MPI_IERR)
      MPI_ALLREDUCE(INTP_ARR,INTPG_ARR,3*NELEM_PP,MPI_REAL,MPI_SUM,&
     &               MPI_COMM_V,MPI_IERR)

      DO IE=1,NELEM_PP !fdgjshfvsfkskjfvs2
        AE=GEOME_PP(7,IE) 
        DO IN=1,NNODE  !fjhsdfkgjhccccfflkjdfgs2
          IDX = (IE-1)*3+IN
          IP=INTMA_PP(IN,IE) 
          ND(IP)=ND(IP)+INTG_ARR(IDX)*AE
          IF(INTG_ARR(IDX).GE.1e-20) THEN
            UVEL_PP(IP)=UVEL_PP(IP)+AE*INTUG_ARR(IDX)/INTG_ARR(IDX) 
            VVEL_PP(IP)=VVEL_PP(IP)+AE*INTVG_ARR(IDX)/INTG_ARR(IDX)
! ***  STORE THE DISCONTINUOUS ND & VELOCITIES TO BE USED IN CONSTRUCTING THE EQUILIBRIUM
! ***  ND FOR THE BGK RHS
            DISND_PP(IN,IE)= INTG_ARR(IDX)
            DISUX_PP(IN,IE)=INTUG_ARR(IDX)/INTG_ARR(IDX)
            DISUY_PP(IN,IE)=INTVG_ARR(IDX)/INTG_ARR(IDX)
          ENDIF

          PS_PP(IP)=PS_PP(IP)+MM*INTPG_ARR*AE
!         DISPS(IN,IE)=INTP*MM-&
!     &     .5*DISND_PP(IN,IE)*MM*(DISUX_PP(IN,IE)*DISUX_PP(IN,IE)+&
!     &        DISUY_PP(IN,IE)*DISUY_PP(IN,IE))
          DISPS(IN,IE)=INTPG_ARR*MM

        ENDDO !  DO IN=1,NNODE  !fjhsdfkgjhccccfflkjdfgs2
      ENDDO !  DO IE=1,NELEM_PP !fdgjshfvsfkskjfvs2


! 
! 
 
! 
! *** NORMALISE THE MACRO VARIABLES 
! 
      DO IP=1,NPOIN_PP  !kdjhfaccdfdaskcjkchda
        IF(ND(IP).LT.1e-20)THEN	!TO AVOID THE 'DIVIDE BY ZERO' 
          RHO_PP(IP)=0.0
          UVEL_PP(IP)=0.0
          VVEL_PP(IP)=0.0
          PS_PP(IP)=0.0
          TEMP_PP(IP)=0.0 
        ELSE
          ND(IP)=ND(IP)/FLAG(IP) 
          RHO_PP(IP)=ND(IP)*MM 
          UVEL_PP(IP)=UVEL_PP(IP)/FLAG(IP) 
          VVEL_PP(IP)=VVEL_PP(IP)/FLAG(IP) 
!          PS_PP(IP)=(PS_PP(IP)/FLAG(IP))-&
!     &   0.5*RHO_PP(IP)*(UVEL_PP(IP)*UVEL_PP(IP)+VVEL_PP(IP)*VVEL_PP(IP)) 
          PS_PP(IP)=PS_PP(IP)/FLAG(IP)
          TEMP_PP(IP)=PS_PP(IP)/(R*RHO_PP(IP)) 
        ENDIF 
      ENDDO !   DO IP=1,NPOIN_PP  !kdjhfaccdfdaskcjkchda
! 
      RETURN 
      END 
