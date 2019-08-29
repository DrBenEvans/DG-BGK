SUBROUTINE NORMAL(ISIDE, COORD, NSIDE, NPOIN, NORMX,&
     &NORMY, EL)
  IMPLICIT NONE
!
! *** DECLARE VARIABLES
!
  INTEGER NSIDE, NPOIN, IP1, IP2, IS
  INTEGER ISIDE(8, NSIDE)
!
  REAL COORD(2, NPOIN), NORMX(NSIDE), NORMY(NSIDE)
  REAL DY(NSIDE), DX(NSIDE), EL(NSIDE)
  REAL C0
! *** ZERO NORMX, NORMY, EL, DY, DX
!
  C0 = 0.0
  CALL RFILLV(NORMX, NSIDE, C0)
  CALL RFILLV(NORMY, NSIDE, C0)
  CALL RFILLV(EL, NSIDE, C0)
  CALL RFILLV(DY, NSIDE, C0)
  CALL RFILLV(DX, NSIDE, C0)
!
! *** FIND THE NORMALS AND LENGTHS OF THE ELEMENT SIDES
!
  DO IS = 1, NSIDE
    IP1 = ISIDE(1, IS)
    IP2 = ISIDE(2, IS)
    DX(IS) = COORD(1, IP2) - COORD(1, IP1)
    DY(IS) = COORD(2, IP2) - COORD(2, IP1)
    EL(IS) = SQRT((DX(IS)*DX(IS)) + (DY(IS)*DY(IS)))
    NORMX(IS) = DY(IS)/EL(IS)
    NORMY(IS) = -DX(IS)/EL(IS)
!
  ENDDO

  RETURN
END
