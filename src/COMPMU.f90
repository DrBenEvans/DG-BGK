      SUBROUTINE COMPMU(rv,VNPNT,NPOIN_PP,NNODE,&
     &          DISNF_PP,VCORD,NELEM_PP,GEOME_PP,NGEOM,&
     &          MU_PP,UX,UY,DISND_PP,DISUX_PP,DISUY_PP,SUMWEIGHT,d,&
     &           MPI_COMM_V,VSPACE_FIRST,VSPACE_LAST)
!
! *** THIS SUBROUTINE COMPUTES THE MU_PP VALUE REQUIRED TO CONSTRUCT THE BGK RHS
!
      IMPLICIT NONE
      include 'mpif.h' 
!
      INTEGER VNPNT,NPOIN_PP,NNODE,NELEM_PP,NGEOM,IE,IP,IV,IN,IVT
!
      REAL rv,d,PI,CO,ETA,ZETA,RT,THETA,CX,CY
      REAL WEIGHT,JAC,SUMWEIGHT,C1,C2
      REAL FRAC,C4,DISUX_PP(NNODE,NELEM_PP),DISUY_PP(NNODE,NELEM_PP)
      REAL CR,UX,UY 
      REAL CX0,CY0,CDASH,SIGT,INTMU,SUMMUV,MUV(VNPNT),VCORD(3,VNPNT) 
      REAL GEOME_PP(NGEOM,NELEM_PP),DISND_PP(NNODE,NELEM_PP),ND 
      REAL MU_PP(NNODE,NELEM_PP)
      REAL DISNF_PP(NNODE,VSPACE_FIRST:VSPACE_LAST,NELEM_PP)

      INTEGER MPI_COMM_V,VSPACE_FIRST,VSPACE_LAST

      REAL CT(3*NELEM_PP)
      REAL CTG(3*NELEM_PP)
      INTEGER IDX



!  ! *** PARAMETERS !
      PARAMETER(PI=3.1416,C4=4.0)
!
! *** CALCULATE THE SCALING FRACTIN FOR V-SPACE INTEGRATION
!
      FRAC=C4/SUMWEIGHT
!
! *** INITIALISE THE MU_PP ARRAY
!
      CO=0.0
      CALL RFILLM(MU_PP,NNODE,NELEM_PP,CO)
!
! *** CALCULATE THE TOTAL COLLISION CROSS SECTION FOR AN OXYGEN MOLECULE
!
      SIGT=PI*d*d
!
! *** LOOP OVER EACH DISCONTINUOUS NODE IN P-SPACE
!
      DO IE=1,NELEM_PP !ljfkajdfkajdfkljvv
        DO IN=1,NNODE !rlfdkfvsdfvvvfdd
          IDX=3*IE+IN
!
! ***     NUMBER DENSITY AT THIS NODE
!
          ND=DISND_PP(IN,IE)
          IF(ND.GE.1e-20)
          THEN
! ***       BULK VELOCITY COMPONENTS
            CX0=DISUX_PP(IN,IE)
            CY0=DISUY_PP(IN,IE)
! ***       LOOP OVER ALL V-SPACE POINTS TO DETERMINE THE MEAN THERMAL VELOCITY
            CT(IDX)=0.0
            DO 1003 IV=VSPACE_FIRST,VSPACE_LAST
              ETA=VCORD(1,IV)
              ZETA=VCORD(2,IV)
              RT=(rv/2)*(ETA+1)
              THETA=ZETA*PI
              CX=RT*COS(THETA)
              CY=RT*SIN(THETA)
              WEIGHT=VCORD(3,IV)
              JAC=PI*rv*rv*0.25*(ETA+1)     !JACOBIAN OF THE TRANSFORMATION
              C1=(CX0-CX)*(CX0-CX)+(CY0-CY)*(CY0-CY)
              CDASH=SQRT(C1)     !RELATIVE VELOCITY
              CT(IDX)=CT(IDX)+&
     &                 CDASH*(DISNF_PP(IN,IV,IE)/ND)*WEIGHT*JAC*FRAC
            ENDDO
          ENDIF
        ENDDO ! DO IN=1,NNODE !rlfdkfvsdfvvvfdd
      ENDDO ! DO IE=1,NELEM_PP !ljfkajdfkajdfkljvv


      MPI_ALLREDUCE(CT,CTG,3*NELEM_PP,MPI_REAL,MPI_SUM,&
     &      MPI_COMM_V,MPI_IERR)

!
! ***   NOW CALCULATE THE MEAN COLLISION FREQUENCY
!
      DO IE=1,NELEM_PP !ljfkajdfkajdfkljvv
        DO IN=1,NNODE !rlfdkfvsdfvvvfdd
          IDX=3*IE+IN
!
          C2=2.0
          CR=SQRT(C2)*CTG(IDX)
          MU_PP(IN,IE)=CR*SIGT*DISND_PP(IN,IE)
!
! ***   END LOOP OVER THE DISCONTINUOUS NODES
!
        ENDDO ! DO IN=1,NNODE !rlfdkfvsdfvvvfdd
      ENDDO ! DO IE=1,NELEM_PP !ljfkajdfkajdfkljvv
!
! *** RETURN TO ADVNCE
!
      RETURN
      END
