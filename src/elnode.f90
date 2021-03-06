SUBROUTINE ELNODE(NELEM, NNODE, NPOIN, INTMA, UNKO, UNKNO)
! *** THIS SUBROUTINE SEPARATES THE INITIAL MESH NODE UNKNOWNS (UNKO) INTO THE
! *** DISCONTINUOUS FORM ELEMENT NODE UNKNOWNS (UNKNO)
  IMPLICIT NONE
!
! ***        DECLARE VARIABLES
!
  INTEGER NELEM, NNODE, NPOIN, INTMA(NNODE, NELEM), IE, IN, IP
  REAL UNKO(NPOIN), UNKNO(1, 3, NELEM)
!
  DO IE = 1, NELEM
    DO IN = 1, NNODE
      IP = INTMA(IN, IE)
      UNKNO(1, IN, IE) = UNKO(IP)
    ENDDO
  ENDDO
!
  RETURN
END
