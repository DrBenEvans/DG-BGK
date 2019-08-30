MODULE DG_LIN_CONVEC_MODULE
  IMPLICIT NONE

CONTAINS

  SUBROUTINE DG_LIN_CONVEC(NDIMN, NNODE, NPOIN, &
                           NELEM, NBOUN, INTMA, &
                           COORD, NTIME, CSAFM, &
                           IMMAT, NBNOI, BSIDO, &
                           NBNOR, NGEOM, RORDER, &
                           TORDER, VNPNT, SUMWEIGHT, &
                           VCORD, rv, ALPHA, &
                           CINF, MXSID, IVD, &
                           FORCEOUT, d, R, &
                           M, MPI_RANK_P, MPI_SIZE_P, &
                           MPI_COMM_P, MPI_RANK_V, MPI_SIZE_V, &
                           MPI_COMM_V, VSPACE_FIRST, VSPACE_LAST)
!
! *** KEYWORDS :
!  - DISCONTINUOUS GALERKIN
!  - SECOND ORDER TAYLOR GALERKIN
!  - TWO-STEP SCHEME
!  - LINEAR FEM
!  - GLOBAL TIMESTEPPING

    USE ADVNCE_MODULE, ONLY: ADVNCE
    USE GETNOR_MODULE, ONLY: GETNOR
    USE INICON2_MODULE, ONLY: INICON2
    USE INPUT_VARIABLE_MODULE, only: INPUTVARIABLES
    USE OUTPUT_MODULE, ONLY: OUTPUT
    USE PTMESH_MODULE, ONLY: PTMESH
    IMPLICIT NONE
    INCLUDE 'mpif.h'

!   DECLARE DYNAMIC ARRAYS
    INTEGER, ALLOCATABLE :: ELGRP(:, :)
    REAL, ALLOCATABLE :: MMAT(:, :)
    REAL, ALLOCATABLE :: DISNF_PP(:, :, :)
    INTEGER, ALLOCATABLE :: NPGRP(:)
    INTEGER, ALLOCATABLE :: NEGRP(:)
    INTEGER, ALLOCATABLE :: GCOMM(:, :)
!
    TYPE(InputVariables) :: IVD
!
! *** DECLARE INTEGER VARIABLES
    INTEGER NSIDE, NDIMN, NNODE, NBNOR, NBNOI, NCOMM_PP
    INTEGER NGEOM, NELEM, NPOIN, NBOUN, NTIME
    INTEGER MXSID, MPI_IERR, IMMAT, VNPNT, RORDER, TORDER
    INTEGER FORCEOUT
    INTEGER maxNELEM_PP, maxNBOUN_PP, maxNSIDE_PP, maxNPOIN_PP
    INTEGER NELEM_PP, NSIDE_PP, NPOIN_PP, NBOUN_PP
    ! mpi-stuff for position space partitioning
    INTEGER MPI_RANK_P, MPI_SIZE_P, MPI_COMM_P
    ! mpi-stuff for velocity space partitioning
    INTEGER MPI_RANK_V, MPI_SIZE_V, MPI_COMM_V
    INTEGER VNPNT_PART, VSPACE_FIRST, VSPACE_LAST

! *** PARAMETERS
    PARAMETER(maxNELEM_PP=30000, maxNBOUN_PP=2400, maxNPOIN_PP=16000)
    PARAMETER(maxNSIDE_PP=81000)
! *** DECLARE REAL VARIABLES MMAT (MASS MATRIX), UX,UY,CSAFM
    REAL UX, UY, CSAFM, ALPHA, DTE
    REAL COORD(2, NPOIN), GEOME(NGEOM, NELEM), RSIDO(3, NBOUN)
    REAL CMMAT(3, 3, NELEM), RELEN(NELEM)
    REAL rv, CINF(4), VCORD(3, VNPNT)
    REAL RHO_PP(maxNPOIN_pp), UVEL_PP(maxNPOIN_PP)
    REAL TEMP_PP(maxNPOIN_PP), NX_PP(maxNSIDE_PP)
    REAL EL_PP(maxNSIDE_PP)
    REAL COORD_PP(2, maxNPOIN_PP), PS_PP(maxNPOIN_PP)
    REAL RSIDO_PP(NBNOR, maxNBOUN_PP), MMAT_PP(NNODE, maxNELEM_PP)
    REAL CMMAT_PP(3, 3, maxNELEM_PP), GEOME_PP(NGEOM, maxNELEM_PP)
    REAL DELUN_PP(1, 3, maxNELEM_pp), UMEAN_PP(maxNELEM_PP)
    REAL NX(mxsid), NY(mxsid)
    REAL FLUXP_PP(1, 3, maxNELEM_pp), FLUYP_PP(1, 3, maxNELEM_pp)
    REAL EL(mxsid), VVEL_PP(maxNPOIN_pp), NY_PP(maxNSIDE_PP)
    REAL SUMWEIGHT, d, R, M
! *** DECLARE INTEGER ARRAYS
    INTEGER ISIDE(8, MXSID)
    INTEGER ISCOM_PP(maxNSIDE_pp)
    INTEGER INTMA(3, NELEM), BSIDO(NBNOI, NBOUN)
    INTEGER IBCOM_PP(maxNBOUN_pp)
    INTEGER LCOMM_PP(maxNSIDE_PP), IPCOM_PP(maxNPOIN_PP)
    INTEGER INTMA_PP(3, maxNELEM_PP), BSIDO_PP(NBNOI, maxNBOUN_PP)
    INTEGER ISIDE_PP(8, maxNSIDE_PP), SDCOM_PP(3, maxNSIDE_PP)
! *** VARIABLES FROM THE GETSID SUBROUTINE
    INTEGER LWHER(NPOIN), LHOWM(NPOIN), ICONE(3*NELEM)
    DOUBLE PRECISION OUTPUT_TIMING

! *** ALLOCATE MEMORY

    VNPNT_PART = VNPNT/MPI_SIZE_V
    ALLOCATE (ELGRP(NELEM, 2), MMAT(NNODE, NELEM))
    ALLOCATE (DISNF_PP(NNODE, VSPACE_FIRST:VSPACE_LAST, maxNELEM_PP))
    ALLOCATE (NPGRP(MPI_SIZE_P), NEGRP(MPI_SIZE_P), GCOMM(MPI_SIZE_P, MPI_SIZE_P))
!
    IF (MPI_RANK_P .EQ. 0) THEN
      ! GET THE ELEMENT SIDES DATA (ROWS 1-4 IN ISIDE)
      CALL GETSID(NELEM, NPOIN, NSIDE, INTMA, ISIDE, LWHER,&
  &                      LHOWM, ICONE, MXSID)
    ENDIF
    !BROADCAST NSIDE TO ALL PROCS
    CALL MPI_BCAST(NSIDE, 1, MPI_INTEGER, 0, MPI_COMM_P, MPI_IERR)

    IF (MPI_RANK_P .EQ. 0) THEN ! sfasdasfsda
! ***   GET THE NORMALS AND LENGTHS OF THE ELEMENT SIDES
      CALL NORMAL(ISIDE, COORD, NSIDE, NPOIN, NX, NY, EL)
! ***   OBTAIN THE GEOMETRY-ARRAYS
      CALL GETGEO(NELEM, NPOIN, NNODE, NDIMN, NGEOM, INTMA, COORD,&
   &              GEOME)
! ***   GET THE REPRESENTATIVE ELEMENT/POINT LENGTHS
      CALL GETLEN(NNODE, NELEM, RELEN,&
   &              NGEOM, GEOME)
! ***   FIND THE ELEMENT LUMPED MASS MATRICES
      IF (IMMAT .EQ. 1) THEN
        CALL GETMAT(NELEM, NNODE, NGEOM, GEOME, MMAT)
      ELSE
        CALL CONMAT(NELEM, NGEOM, GEOME, CMMAT)
      ENDIF
! ***   GET THE NORMALS OF THE MESH

      CALL GETNOR(NDIMN=NDIMN, NBNOI=NBNOI, NBNOR=NBNOR, &
                  NPOIN=NPOIN, NBOUN=NBOUN, BSIDO=BSIDO, &
                  RSIDO=RSIDO, COORD=COORD )

    ENDIF !  IF(MPI_RANK_P.EQ.0)THEN ! sfasdasfsda
!
! *** CONVERT FROM METIS NODAL PARTITION TO ELEMENT PARTITION
!
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)

!     The PTMESH routine does not care about velocity space and velocity
!     space partitioning.
    CALL PTMESH(NELEM=NELEM, INTMA=INTMA, NPOIN=NPOIN, &
                MPI_SIZE_P=MPI_SIZE_P, NSIDE=NSIDE, ISIDE=ISIDE, &
                NELEM_PP=NELEM_PP, NPOIN_PP=NPOIN_PP, NBOUN_PP=NBOUN_PP, &
                NSIDE_PP=NSIDE_PP, NEGRP=NEGRP, LCOMM_PP=LCOMM_PP, &
                GCOMM=GCOMM, ELGRP=ELGRP, &
                NGEOM=NGEOM, NBNOI=NBNOI, NBOUN=NBOUN, &
                NBNOR=NBNOR, NNODE=NNODE, NX=NX, &
                NY=NY, EL=EL, GEOME=GEOME, &
                MMAT=MMAT, CMMAT=CMMAT, COORD=COORD, &
                BSIDO=BSIDO, RSIDO=RSIDO, maxNSIDE_PP=maxNSIDE_PP, &
                maxNBOUN_PP=maxNBOUN_PP, maxNELEM_PP=maxNELEM_PP, maxNPOIN_PP=maxNPOIN_PP, &
                NX_PP=NX_PP, NY_PP=NY_PP, EL_PP=EL_PP, &
                GEOME_PP=GEOME_PP, MMAT_PP=MMAT_PP, CMMAT_PP=CMMAT_PP, &
                INTMA_PP=INTMA_PP, COORD_PP=COORD_PP, BSIDO_PP=BSIDO_PP, &
                RSIDO_PP=RSIDO_PP, IMMAT=IMMAT, MPI_RANK_P=MPI_RANK_P, &
                IPCOM_PP=IPCOM_PP, IBCOM_PP=IBCOM_PP, ISCOM_PP=ISCOM_PP, &
                NPGRP=NPGRP, ISIDE_PP=ISIDE_PP, NCOMM_PP=NCOMM_PP, &
                SDCOM_PP=SDCOM_PP, IVD=IVD, MPI_COMM_P=MPI_COMM_P)

!
! *** SET UP THE INITIAL CONDITIONS OF THE PROBLEM (FOR A GAS EXPANSION)
!
    IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN ! write only on master rank
      WRITE (*, *) "Initialization..."
    ENDIF
!      CALL INICON(NPOIN,NNODE,NELEM,NELEM_PP,INTMA_PP,ELGRP,&
!     &NPOIN_PP,COORD_PP,VNPNT,VCORD,DISNF_PP,IPCOM_PP,rv,&
!     &                 IVD,&
!     &                 MPI_RANK_P,MPI_SIZE_P,MPI_COMM_P,&
!     &                 MPI_RANK_V,MPI_SIZE_V,MPI_COMM_V,&
!     &                 VSPACE_FIRST,VSPACE_LAST)
!
!      CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERR)
!
! *** FOR A UNIFORM INITIAL CONDITION
!

    CALL INICON2(NPOIN=NPOIN, NNODE=NNODE, NELEM=NELEM, &
                 NELEM_PP=NELEM_PP, ELGRP=ELGRP, &
                 VNPNT=VNPNT, VCORD=VCORD, &
                 DISNF_PP=DISNF_PP, CINF=CINF, &
                 rv=rv, IVD=IVD, GC=R, &
                 M=M, MPI_RANK_P=MPI_RANK_P, &
                 MPI_COMM_P=MPI_COMM_P, &
                 VSPACE_FIRST=VSPACE_FIRST, VSPACE_LAST=VSPACE_LAST)

! *** SYNCHRONISE
!
    CALL MPI_BARRIER(MPI_COMM_WORLD, MPI_IERR)
!
! *** PERFORM THE TIME STEP ITERATIONS
!
    CALL ADVNCE(NNODE=NNODE, NGEOM=NGEOM, NBNOI=NBNOI,&
      NBOUN=NBOUN, NBNOR=NBNOR, COORD=COORD,&
      BSIDO=BSIDO, NBOUN_PP=NBOUN_PP, NPOIN=NPOIN,&
      NPOIN_PP=NPOIN_PP, NELEM_PP=NELEM_PP, NELEM=NELEM,&
      NTIME=NTIME, BSIDO_PP=BSIDO_PP, INTMA_PP=INTMA_PP,&
      FLUXP_PP=FLUXP_PP, FLUYP_PP=FLUYP_PP, GEOME_PP=GEOME_PP,&
      MMAT_PP=MMAT_PP, CMMAT_PP=CMMAT_PP, RELEN=RELEN,&
      DTE=DTE, DELUN_PP=DELUN_PP, RSIDO_PP=RSIDO_PP,&
      UMEAN_PP=UMEAN_PP, IMMAT=IMMAT, CSAFM=CSAFM,&
      UX=UX, UY=UY, NSIDE_PP=NSIDE_PP,&
      ISIDE_PP=ISIDE_PP, NX_PP=NX_PP, NY_PP=NY_PP,&
      EL_PP=EL_PP, DISNF_PP=DISNF_PP, VNPNT=VNPNT,&
      SUMWEIGHT=SUMWEIGHT, rv=rv, VCORD=VCORD,&
      RHO_PP=RHO_PP, UVEL_PP=UVEL_PP, VVEL_PP=VVEL_PP,&
      PS_PP=PS_PP, TEMP_PP=TEMP_PP, ALPHA=ALPHA,&
      CINF=CINF, GCOMM=GCOMM, IPCOM_PP=IPCOM_PP,&
      NCOMM_PP=NCOMM_PP, maxNELEM_PP=maxNELEM_PP, RORDER=RORDER,&
      TORDER=TORDER, maxNPOIN_PP=maxNPOIN_PP, SDCOM_PP=SDCOM_PP,&
      IVD=IVD, FORCEOUT=FORCEOUT, d=d,&
      RGas=R, M=M, MPI_RANK_P=MPI_RANK_P,&
      MPI_SIZE_P=MPI_SIZE_P, MPI_COMM_P=MPI_COMM_P, MPI_RANK_V=MPI_RANK_V,&
      MPI_COMM_V=MPI_COMM_V, VSPACE_FIRST=VSPACE_FIRST,  VSPACE_LAST=VSPACE_LAST)

! *** OUTPUT FOR GID AND RESTART FILE
!
    OUTPUT_TIMING = MPI_WTIME()
!    CALL OUTPUT(NNODE=NNODE, NPOIN=NPOIN, NELEM=NELEM, &
!                maxNELEM_PP=maxNELEM_PP, COORD=COORD, INTMA=INTMA, &
!                DISNF_PP=DISNF_PP, VNPNT=VNPNT, ELGRP=ELGRP, &
!                NEGRP=NEGRP, IVD=IVD, MPI_RANK_P=MPI_RANK_P, &
!                MPI_SIZE_P=MPI_SIZE_P, MPI_COMM_P=MPI_COMM_P, MPI_RANK_V=MPI_RANK_V, &
!                MPI_SIZE_V=MPI_SIZE_V, MPI_COMM_V=MPI_COMM_V, VSPACE_FIRST=VSPACE_FIRST, &
!                VSPACE_LAST=VSPACE_LAST, VCORD=VCORD, rv=rv)

    OUTPUT_TIMING = MPI_WTIME() - OUTPUT_TIMING
    IF ((MPI_RANK_P .EQ. 0) .AND. (MPI_RANK_V .EQ. 0)) THEN ! write only on master rank
      WRITE (*, "(A20,F6.2,A3)") "OUTPUT written in ", OUTPUT_TIMING,&
  &                             "sec"
    ENDIF
!
    RETURN
  END
END MODULE DG_LIN_CONVEC_MODULE
