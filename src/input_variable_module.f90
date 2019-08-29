module INPUT_VARIABLE_MODULE
  implicit none

  TYPE InputVariables
    INTEGER :: TORDER !Quadrature order in Theta for V-SPACE
    INTEGER :: NTIME !Number of timesteps
    INTEGER :: FORCEOUT !Force output at each timestep (1 = yes, 0 = no)
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
    CHARACTER(len=80) :: LobattoFile
    CHARACTER(len=80) :: PSpaceFile
    CHARACTER(len=80) :: OutFile
    CHARACTER(len=80) :: PartitionFile
    CHARACTER(len=80) :: RestartInFile
    CHARACTER(len=80) :: ResidualFile
    CHARACTER(len=80) :: ResultsFile1
    CHARACTER(len=80) :: ResultsFile2
    CHARACTER(len=80) :: RestartOutFile
    CHARACTER(len=80) :: GIDMeshFile
  END TYPE

end Module INPUT_VARIABLE_MODULE
