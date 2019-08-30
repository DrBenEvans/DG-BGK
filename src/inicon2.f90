MODULE INICON2_MODULE

CONTAINS

  SUBROUTINE INICON2(NPOIN, NNODE, NELEM, &
                     NELEM_PP, ELGRP, VNPNT, &
                     VCORD, DISNF_PP, CINF, &
                     rv, IVD, GC, &
                     M, MPI_RANK_P, MPI_COMM_P, &
                     VSPACE_FIRST, VSPACE_LAST)
!
! *** THIS SUBROUTINE ESTABLISHES THE INITIAL CONDITIONS OF THE PROBLEM
!
    USE INPUT_VARIABLE_MODULE, ONLY: INPUTVARIABLES
    IMPLICIT NONE
    INCLUDE 'mpif.h'
!
! CREATE STRUCTURE FOR INPUT VARIABLES
!

! *** DECLARE VARIABLES
!
    TYPE(InputVariables) :: IVD
    INTEGER IN, INV, RS, IV, I, NELEM
    INTEGER VNPNT, NPOIN, NNODE, NELEM_PP
    INTEGER TEST1, TEST2, IE_PP, IE
    INTEGER ELGRP(NELEM, 2), MPI_IERR
    ! mpi-stuff for position space partitioning
    INTEGER MPI_RANK_P, MPI_COMM_P
    ! mpi-stuff for velocity space partitioning
    INTEGER VSPACE_FIRST, VSPACE_LAST
    INTEGER NPART_P

    REAL VCORD(3, VNPNT)
    REAL T1, P1, RHO1, n1, NA
    REAL M, BETA1, TMP1, U0, V0
    REAL UX, UY, TMP, SPEED, CO, PI, F0
    REAL ETA, ZETA, R, THETA, rv, GC
    REAL CINF(4)
!
    REAL DISNF_PP(NNODE, VSPACE_FIRST:VSPACE_LAST, NELEM_PP)
    REAL DISNFPARCEL(NNODE)

!
    CHARACTER filename*80, TEXT*80
!
! *** SET PI, AVAGADRO'S NUMBER AS PARAMETERS
!
    PARAMETER(NA=6.022E+026, PI=3.1416)
!
! *** INITIALISE DISNF_PP
!
    DISNF_PP = 0.0
!
! *** ASK THE USER IF THEY ARE USING A RESTART FILE AS INITIAL CONDITIONS
!
    RS = IVD%RS
!
! *** IF USING A RESTART FILE ASK FOR THE NAME OF THE FILE
!
    IF (RS .EQ. 1) THEN !dlafcccnwweppsdc
      PRINT *, 'I AM READING FROM THE RESTART FILE'
      IF (MPI_RANK_P .EQ. 0) THEN !dlssprrertepi
        filename = IVD%RestartInFile
        PRINT *, 'READING RESTART FILE = ', IVD%RestartInFile
        OPEN (11, file=filename, status='old')
        WRITE (*, *)
!
! ***     READ IN DATA FROM RESTRT FILE
!
        READ (11, *) TEXT
        READ (11, *) TEST1, TEST2
!
! ***     CHECK DATA
!
        IF ((TEST1 .NE. NPOIN) .OR. (TEST2 .NE. VNPNT)) THEN !kjlkfjsdklalaa
          WRITE (*, 208)
          WRITE (*, 209)
          STOP
        ENDIF !IF((TEST1.NE.NPOIN).OR.(TEST2.NE.VNPNT))THEN !kjlkfjsdklalaa
!
      ENDIF ! IF(MPI_RANK_P.EQ.0)THEN !dlssprrertepi
!
! ***   READ DISNF ONTO MASTER AND DISTRIBUTE
!
      ! THIS IS SUPER SLOW BUT FAST TO IMPLEMENT
      DO IV = 1, VNPNT
        DO IE = 1, NELEM
          IF (MPI_RANK_P .EQ. 0) THEN !flkjdssssssyy
            READ (11, *) (DISNFPARCEL(I), I=1, NNODE)
! ***         CHECK THE SLAVE ALLOCATION OF THIS ELEMENT
            NPART_P = ELGRP(IE, 1)
            IE_PP = ELGRP(IE, 2)
          ENDIF !flkjdssssssyy
! ***       BROADCAST
          IF ((IV .GE. VSPACE_FIRST) .AND. (IV .LE. VSPACE_LAST)) THEN !dsfkjshd
            CALL MPI_BCAST(NPART_P, 1, MPI_INTEGER, 0,&
   &                    MPI_COMM_P, MPI_IERR)
            CALL MPI_BCAST(IE_PP, 1, MPI_INTEGER, 0,&
   &                    MPI_COMM_P, MPI_IERR)
            CALL MPI_BCAST(DISNFPARCEL, NNODE, MPI_REAL, 0,&
   &                    MPI_COMM_P, MPI_IERR)
            IF ((MPI_RANK_P + 1) .EQ. NPART_P) THEN !fldjssmfslms
              DO IN = 1, NNODE
                DISNF_PP(IN, IV, IE_PP) = DISNFPARCEL(IN)
              ENDDO
            ENDIF  ! fldjssmfslms
          ENDIF !dsfkjshd
        ENDDO
      ENDDO
!
! ***   CLOSE RESTART FILE
!
      IF (MPI_RANK_P .EQ. 0) CLOSE (11)
!
    ELSE ! IF(RS.EQ.1)THEN !dlafcccnwweppsdc
!
! ***   GET THE INITIAL CONDITIONS:
!
      IF (MPI_RANK_P .EQ. 0) THEN !difadocaodmcadc
        T1 = IVD%T1
        P1 = IVD%P1
        U0 = IVD%U0
        V0 = IVD%V0
!
! ***     CALCULATE GAS PROPERTIES USED FOR MAXWELL DISTRIBUTION
!
        RHO1 = (P1*(10**05))/(GC*T1)        !DENSITIES
        n1 = RHO1*NA/M
        TMP1 = RHO1/(2*P1*(10**05))
        BETA1 = SQRT(TMP1)                !BETA VARIABLE
      ENDIF !IF(MPI_RANK_P.EQ.0)THEN !difadocaodmcadc
      CALL MPI_BCAST(BETA1, 1, MPI_REAL, 0, MPI_COMM_P, MPI_IERR)
      CALL MPI_BCAST(n1, 1, MPI_REAL, 0, MPI_COMM_P, MPI_IERR)
      CALL MPI_BCAST(U0, 1, MPI_REAL, 0, MPI_COMM_P, MPI_IERR)
      CALL MPI_BCAST(V0, 1, MPI_REAL, 0, MPI_COMM_P, MPI_IERR)
!
! ***   BEGIN LOOP OVER THE PHYSICAL SPACE DISCONTINUOUS NODES
!
!
! ***   IF VACCUUM FOR INTITIAL CONDITION
!
      IF (n1 .GE. 1e-20) THEN !djfhadcpaaaargekj
        DO IE = 1, NELEM_PP
          DO IN = 1, NNODE
!
! ***         LOOP OVER ALL VELOCITY SPACE NODES
!
            CO = (BETA1**2)/(PI)            !COEFFICIENT OF THE MAXWELL DISTRIBUTION FUNCTION
            DO INV = VSPACE_FIRST, VSPACE_LAST
              ETA = VCORD(1, INV)
              ZETA = VCORD(2, INV)
!               TRANSFORM FROM ETA-ZETA COORDS TO R-THETA COORDS
              R = ETA*(rv/2) + (rv/2)
              THETA = ZETA*PI
!               TRANSORM TO CARTESIANS
              UX = R*COS(THETA)
              UY = R*SIN(THETA)
              TMP = (UX - U0)*(UX - U0) + (UY - V0)*(UY - V0)
              SPEED = SQRT(TMP)     !MOLECULAR SPEED AT THIS COORDINATE IN VELSPACE MESH
              TMP = -((BETA1**2)*(SPEED**2))
              F0 = CO*EXP(TMP)         !DISTRIBUTION FUNCTION
              DISNF_PP(IN, INV, IE) = n1*F0!INITIAL CONDITION (nf) MATRIX
!
! ***           END LOOP OVER VELOCITY SPACE NODES
!
            ENDDO
!
! ***         END LOOP OVER THE PHYSICAL SPACE DISCONTINUOUS NODES
!
          ENDDO
!
        ENDDO
      ENDIF  ! IF(n1.GE.1e-20)THEN !djfhadcpaaaargekj
!
! ***   SYNCHRONISE PROCESSORS
!
      CALL MPI_BARRIER(MPI_COMM_P, MPI_IERR)
!
    ENDIF  !   ELSE ! IF(RS.EQ.1)THEN !dlafcccnwweppsdc

!
! *** ASK WHETHER THERE IS AN INFLOW PRESENT?
!
    I = IVD%INF
    IF (I .EQ. 1) THEN
      CINF = IVD%CINF
    ELSE
      CINF(1) = 0.0
      CINF(2) = 0.0
      CINF(3) = 0.0
      CINF(4) = 0.0
    ENDIF
!
! *** FORMAT STATEMENTS
!
208 FORMAT('THIS RESTART FILE IS NOT COMPATIBLE WITH YOUR MESHES!')
209 FORMAT('PROGRAM STOPPED')
!
    RETURN
  END

END MODULE INICON2_MODULE
