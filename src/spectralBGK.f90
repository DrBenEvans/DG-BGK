! *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c 
!					2D DISCONTINUOUS GALERKIN FE		    c 
!						 SPECTRAL BGK                       c 
!				by Ben Evans, Swansea Unviersity		    c 
! *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c *** c 
! 
! 
! *** THIS IS THE SPECTRAL V-SPACE VERSION (APRIL 2006) 
! 
     PROGRAM MAIN 
! 
     IMPLICIT NONE 
     INCLUDE 'mpif.h' 
! 
! *** THE WORKING FILES ARE  
!			- LOBATTO WEIGHTING FILE 1   ...........CHANNEL 10
!			- READ IN RESTART FILE ...............CHANNEL 11 
!			- METIS P-SPACE PARTITION FILE........CHANNEL 12 
!	        	- CONFIGURATION FILE .................CHANNEL 13 
!			- PHYSICAL SPACE INPUT FILE ..........CHANNEL 14 
!                       - OUTPUT FILE ........................CHANNEL 15 
!                       - RESIDUAL FILE.......................CHANNEL 16 
!                       - RESULTS FILE 1......................CHANNEL 17 
!		       	- RESULTS FILE 2 .....................CHANNEL 18 
!			- GID MESH FILE ......................CHANNEL 19 
!			- NEW RESTART FILE ...................CHANNEL 20 
!			- INPUT FILE .........................CHANNEL 21 
!
! CREATE STRUCTURE FOR INPUT VARIABLES
!
        TYPE InputVariables
          INTEGER :: TORDER !Quadrature order in Theta for V-SPACE
          INTEGER :: NTIME !Number of timesteps
          INTEGER :: FORCEOUT !Force output at each timestep (1 = yes, 0 = no)
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
! *** DECLARE DYNAMIC ARRAYS 
! 
        INTEGER, ALLOCATABLE :: INTMA(:,:) 
        INTEGER, ALLOCATABLE :: BSIDO(:,:) 
        INTEGER, ALLOCATABLE :: IELSI(:,:) 
        REAL, ALLOCATABLE :: COORD(:,:) 
        REAL, ALLOCATABLE :: SIZE(:) 
        REAL, ALLOCATABLE :: RHO(:) 
        REAL, ALLOCATABLE :: UVEL(:) 
        REAL, ALLOCATABLE :: VVEL(:) 
        REAL, ALLOCATABLE :: PS(:) 
        REAL, ALLOCATABLE :: TEMP(:) 
        REAL, ALLOCATABLE :: VCORD(:,:) 
! 
! *** DECLARE INTEGER VARIABLES 
! 
        INTEGER NLINES,NDIMN,NNODE,NBNOR,NBNOI 
        INTEGER NGEOM,NELEM,NPOIN,NBOUN,MXSID ,IP
        INTEGER ILINE,NTIME,VNPNT,IV,ORDER2,VNPNT2,FORCEOUT
        INTEGER RORDER,TORDER !THESE ARE THE ORDERS OF THE R AND THETA DISCRETISATIONS RESPECTIVELY IN V-SPACE
        INTEGER IMMAT,maxpn_pp,maxbn_pp 
        INTEGER MPI_RANK,MPI_SIZE,MPI_IERR 
! 
! *** DECLARE REAL VARIABLES 
! 
        REAL CINF(4) 
        REAL CSAFM, rv,ALPHA,SUMWEIGHT 
        REAL d,R,M
! 
! *** DECLARE CHARACTER VARIABLES 
! 
        CHARACTER TEXT*80, filename*80 
!
! *** DECLARE INPUT VARIABLES DATASTRUCTURE
!
      TYPE(InputVariables) :: IVD
            NAMELIST / Input / IVD
! 
! *** INITIALISE MPI 
! 
      CALL MPI_INIT(MPI_IERR) !Initialize the MPI execution environment 
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MPI_RANK,MPI_IERR) !Determines the rank of the calling process in the communicator 
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,MPI_SIZE,MPI_IERR) !Determines the size of the group associated with the communictor  
! 
! *** DISPLAY PROGRAM HEADER TO THE SCREEN 
! 
          IF(MPI_RANK.EQ.0)THEN !PROCESSOR 0 ONLY 
            WRITE(*,*) 
            WRITE(*,*) 
            WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' 
            WRITE(*,*)'%     2D SPECTRAL BGK BOLTZMANN SOLVER      %' 
            WRITE(*,*)'%                      USING                %' 
            WRITE(*,*)'%           DISCONTINUOUS GALERKIN FE       %' 
            WRITE(*,*)'%          (PARALLEL VERSION)               %' 
            WRITE(*,*)'%      running on',MPI_SIZE,'procs          %' 
            WRITE(*,*)'%             written by                    %' 
            WRITE(*,*)'%               Ben Evans                   %' 
            WRITE(*,*)'%                                           %' 
            WRITE(*,*)'% Civil & Computational Engineering Centre  %' 
            WRITE(*,*)'%         University of Wales, Swansea      %' 
            WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' 
            WRITE(*,*) 
            WRITE(*,*) 
            WRITE(*,*) 
! 
! *** OPEN THE CONFIGURATION FILE (filename stored in INPUT FILE) 
! 
            filename='run.inp' 
! 
            OPEN(21,file=filename,status='old')
            READ(21,*) filename 
            CLOSE(21)
!
            PRINT*,'CONFIGURATION FILE = ',filename 
            OPEN (13,file=filename,status='old') 
            WRITE(*,*)   
!
! *** READ IN ALL THE INPUT VARIABLES AS A NAMELIST FROM THE CONFIGURATION FILE
!
            IVD%TORDER = 20
            IVD%NTIME = 10000
            IVD%FORCEOUT = 0
            IVD%IMMAT = 1 
            IVD%INF = 1
            IVD%RS = 0
            IVD%CSAFM = 0.5
            IVD%rv = 2000
            IVD%T1 = 293
            IVD%P1 = 0.1
            IVD%U0 = 100
            IVD%V0 = 0.0
            IVD%W = 1.0
            IVD%ALPHA = 0.9
            IVD%R = 287
            IVD%d = 250e-12
            IVD%M = 32
            IVD%CINF(1) = 294 !INFLOW TEMP
            IVD%CINF(2) = 0.001473 !INFLOW PRESSURE
            IVD%CINF(3) = 171.85 !INFLOW U-VEL
            IVD%CINF(4) = 0.0   !INFLOW V-VEL
            IVD%LobattoFile = 'Lobatto20.txt' !Lobatto Quadrature File
            IVD%PSpaceFile = 'AEROFOIL.RES' !PSpace Mesh File
            IVD%OutFile = 'FILE.OUT' !Output File
            IVD%PartitionFile = 'METIS3WAYPARTAEROFOIL13.RES' !METIS PARTITION FILE
            IVD%RestartInFile = 'RESTART.RES' !RESTART INPUT FILE
            IVD%ResidualFile = 'RESIDUAL.RES' !RESIDUAL FILE
            IVD%ResultsFile1 = 'RESULTS1.RES' !RESULTS FILE 1
            IVD%ResultsFile2 = 'RESULTS2.RES' !RESULTS FILE 2
            IVD%RestartOutFile = 'RESTART.RES' !RESTART OUTPUT FILE
            IVD%GIDMeshFile = 'GIDMESH.RES' !GID MESH FILE
!
            READ(13,Input)
            CLOSE(13)
! 
! *** OPEN THE LOBATTO WEIGHTINGS FILES 
! *** FIRST FOR FULL VMESH INTEGRATION
            filename = IVD%LobattoFile
            PRINT*,'READING V-SPACE LOBATTO FILE = ',IVD%LobattoFile
            OPEN (10, file=filename, status='old') 
            WRITE(*,*)
! 
! *** READ IN THE LOBATTO ORDER 
!	 
            READ (10,*) RORDER 
            TORDER = IVD%TORDER
! 
! *** CALCULATE THE NUMBER OF NODES IN VSPACE 
! 
            VNPNT=RORDER*TORDER
            ALLOCATE(VCORD(3,VNPNT))  
! 
! *** CALL VSPACE SUBROUTINE 
! 
            CALL VSPACE(RORDER,TORDER,VNPNT,VCORD,SUMWEIGHT) 
            rv = IVD%rv
            FORCEOUT = IVD%FORCEOUT
! 
          ENDIF
! 
! *** MPI BCASTS 
! 
          CALL MPI_BCAST(VNPNT,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR)
          IF(MPI_RANK.NE.0)THEN
            ALLOCATE(VCORD(3,VNPNT))
          ENDIF 
          CALL MPI_BCAST(rv,1,MPI_REAL,0,MPI_COMM_WORLD,MPI_IERR) 
          CALL MPI_BCAST(FORCEOUT,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR)
          CALL MPI_BCAST(SUMWEIGHT,1,MPI_REAL,0,MPI_COMM_WORLD,MPI_IERR)
          CALL MPI_BCAST(RORDER,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR)
          CALL MPI_BCAST(TORDER,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERR)
          CALL MPI_BCAST(VCORD,3*VNPNT,MPI_REAL,0,MPI_COMM_WORLD,MPI_IERR)
!
! *** OPEN PHYSICAL SPACE MESH DATA FILE AND READ GLOBAL PARAMETERS 
! 
      IF(MPI_RANK.EQ.0)THEN  !sdfhasdcabadfadaaa
! 
! *** OPEN THE INPUT FILE ON CHANNEL5 
! 
        filename = IVD%PSpaceFile
        PRINT*,'READING P-SPACE MESH FILE = ',IVD%PSpaceFile 
        OPEN  (14, file=filename, status='old')  
        WRITE(*,*) 

        WRITE(*,*) 
        filename = IVD%OutFile
        PRINT*,'OPENDING OUTPUT FILE = ',IVD%OutFile
        OPEN  (15, file=filename, status='UNKNOWN')    
        WRITE(*,*) 
        WRITE(*,*) 

        NTIME = IVD%NTIME 
        IMMAT = IVD%IMMAT
        ALPHA = IVD%ALPHA
        CSAFM = IVD%CSAFM
        d = IVD%d
        R = IVD%R
        M = IVD%M
! 
! *** READ (FROM INPUT FILE) AND WRITE (TO OUTPUT FILE) TITLE OF COMPUTATION 
!  
        READ(14,*) NLINES 
! 
        DO 400 ILINE=1,NLINES 
          READ(14,*) TEXT 
          WRITE(15,1)TEXT 
  400   CONTINUE 
! 
! *** SPECIFY DIMENSION 
! 
        READ(14,*) TEXT 
        WRITE(15,1)TEXT 
        READ(14,*) NDIMN 
        WRITE(15,5)NDIMN 
! 
! *** FROM THIS OBTAIN  
!          NNODE : NR. OF NODES PER ELEMENT 
!          NBNOR : ROWS OF DATA FOR REAL BOUNDARY ARRAY 
!          NBNOI : ROWS OF DATA FOR INT. BOUNDARY ARRAY 
!          NGEOM : N,X N,Y N,Z RJAC  AT ELEMENT LEVEL 
!     
        NNODE=NDIMN+1 
        NBNOR=NDIMN+1 
        NBNOI=NDIMN+4 
        NGEOM=1+NDIMN*NNODE 
! 
! 
! *** SPECIFY: NR. OF ELEMENTS , NR. OF NODAL POINTS ,NR. OF BOUNDARY NO 
! 
        READ(14,*) TEXT 
        WRITE(15,1)TEXT 
        READ(14,*) NELEM,NPOIN,NBOUN 
        WRITE(15,5)NELEM,NPOIN,NBOUN  
! *** SET VALUE OF MXSID 
        mxsid=2*(NELEM+NBOUN) ! (3*NELEM+NBOUN)/2 IS ENOUGH,
                              ! YOU MUPPET
! 
! *** WRITE V-SPACE INFO TO OUTPUT FILE 
!      
        WRITE(15,170) 
        WRITE(15,180) 
        WRITE(15,5) RORDER,TORDER,VNPNT 
! 
! *** ALLOCATE DIMENSIONS TO ARRAYS 
! 
        ALLOCATE (INTMA(3,NELEM),BSIDO(NBNOI,NBOUN),IELSI(2,NELEM)) 
        ALLOCATE (COORD(2,NPOIN),RHO(NPOIN)) 
        ALLOCATE (UVEL(NPOIN),VVEL(NPOIN),PS(NPOIN),TEMP(NPOIN)) 
! 
! *** READ AND WRITE ALL GLOBAL DATA 
! 
        CALL GTINPT( NDIMN ,NNODE ,NPOIN ,NELEM ,& 
     &                  NBOUN ,INTMA ,COORD , BSIDO ,& 
     &                  IELSI , NBNOI) 
      ENDIF ! IF(MPI_RANK.EQ.0)THEN  !sdfhasdcabadfadaaa
! 
! *** BROADCAST THE RELEVANT P-SPACE DATA TO THE SLAVE PROCESSORS 
! 
       CALL BCASTP(NDIMN,NNODE,NBNOR,NBNOI,NGEOM,CSAFM,& 
     &    IMMAT,ALPHA,NELEM,NPOIN,NBOUN,MXSID,d,R,M) 
!
      IF(MPI_RANK.NE.0)THEN
        ALLOCATE (INTMA(3,NELEM),BSIDO(NBNOI,NBOUN),IELSI(2,NELEM)) 
        ALLOCATE (COORD(2,NPOIN),RHO(NPOIN)) 
        ALLOCATE (UVEL(NPOIN),VVEL(NPOIN),PS(NPOIN),TEMP(NPOIN)) 
      ENDIF 
! 
! *** TELL THE USER WHETHER THE CODE IS USING LUMPER OR CONSISTENT MMAT 
! 
       IF(MPI_RANK.EQ.0)THEN 
         IF(IMMAT.EQ.0)THEN 
          PRINT*,'NOTE: YOU HAVE CHOSEN TO USE CONSISTENT MASS MATRICES' 
         ELSE 
           PRINT*,'NOTE: YOU HAVE CHOSEN TO USE LUMPED MASS MATRICES' 
         ENDIF 
         WRITE(*,*) 
         WRITE(*,*) 
! 
! *** INITIALISE RHO(NPOIN),UVEL(NPOIN),VVEL(NPOIN),TEMP(NPOINT),PS(NPOIN) 
! 
! 
         CALL RFILLV(RHO,NPOIN,0.0) 
         CALL RFILLV(UVEL,NPOIN,0.0) 
         CALL RFILLV(VVEL,NPOIN,0.0) 
         CALL RFILLV(TEMP,NPOIN,0.0) 
         CALL RFILLV(PS,NPOIN,0.0)   
! 
! *** CLOSE CHANNELS 14 and 15
! 
         CLOSE(14) 
         CLOSE(15) 
       ENDIF  
! 
! *** CALL DG_LIN_CONVEC FOR THE DISTRIBUTION FUNCTION CONVECTION PROCESS 
! *** (CALLED BY ALL PROCESSORS RANK 0 -> MPI_SIZE-1) 
       CALL DG_LIN_CONVEC(NDIMN ,NNODE ,NPOIN ,NELEM ,& 
     &                  NBOUN ,INTMA ,COORD ,& 
     &                  NTIME ,CSAFM ,& 
     &                  IMMAT ,NBNOI ,BSIDO ,& 
     &                  NBNOR,NGEOM,RORDER,TORDER,& 
     &                  VNPNT,SUMWEIGHT,VCORD,& 
     &                  rv,RHO,UVEL,VVEL,& 
     &                  PS,TEMP,ALPHA,CINF,MXSID,MPI_RANK,& 
     &                  MPI_SIZE,IVD,FORCEOUT,d,R,M) 
! 
! *** END OF PROGRAM 
! 
      IF(MPI_RANK.EQ.0)THEN 
        PRINT*,'THANK YOU FOR USING BENS BOLTZMANN SOLVER' 
        PRINT*,'  I HOPE YOU HAVE ENJOYED THE EXPERIENCE' 
        PRINT*,'      PLEASE RETURN SOON				' 
      ENDIF 
! 
! *** SHUT DOWN MPI ENVIRONMENT 
! 
      CALL MPI_FINALIZE(MPI_IERR)  
! 
! *** FORMAL STATEMENTS 
! 
    1 FORMAT(20A30) 
    5 FORMAT(I6,I6,I6) 
  150 FORMAT('WHAT IS THE NAME OF THE PHYSICAL SPACE MESH FILE?  ',$) 
  160 FORMAT('WHAT WOULD YOU LIKE THE OUTPUT FILE TO BE CALLED?& 
     & (____________.RES) :  ',$) 
  170 FORMAT('V-SPACE STATS:') 
  180 FORMAT('LOBATTO ORDER      TOTAL NUMBER OF NODES') 
  201  FORMAT('WHAT IS THE NAME OF THE VELOCITY SPACE MESH FILE?  ', $) 
  202  FORMAT(20A4) 
  203  FORMAT(10I6) 
! 
      STOP 
      END 
