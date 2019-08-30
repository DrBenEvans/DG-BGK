      SUBROUTINE GETRH0(NELEM, NPOIN, NNODE,&
     &                  INTMA, GEOME, DELUN, RHS)
        IMPLICIT NONE

        INTEGER :: NELEM, NPOIN, NNODE, IELEM, IN, INODE, JN
        REAL :: CM, RJAC, RJ
!
        REAL DELUN(NPOIN), RHS(NPOIN), GEOME(7, NELEM)
        REAL RHSEL(3), DELEL(3), MMAT(3, 3)
!
        INTEGER INTMA(NNODE, NELEM), NODEN(9)
!
        DATA MMAT/2.0, 1.0, 1.0,&
       &            1.0, 2.0, 1.0,&
       &            1.0, 1.0, 2.0/
!
! *** SET RHS=0
!
        CALL RFILLV(RHS, NPOIN, 0.0)
!
! *** LOOP OVER THE ELEMENTS
!
        DO IELEM = 1, NELEM
!
! *** PICK UP THE VALUES NEEDED
!
          DO INODE = 1, NNODE
            IN = INTMA(INODE, IELEM)
            NODEN(INODE) = IN
            ! *** NODAL VALUES OF DELUN NEEDED
            DELEL(INODE) = DELUN(IN)
          ENDDO
          ! *** JACOBIAN NEEDED
          RJAC = GEOME(7, IELEM)
          RJ = -RJAC/24.

! *** OBTAIN THE ELEMENT CONTRIBUTION
          CALL RFILLV(RHSEL, 3, 0.0)
!
          DO IN = 1, 3
          DO JN = 1, 3
            CM = RJ*MMAT(IN, JN)
            RHSEL(IN) = RHSEL(IN) + CM*DELEL(JN)
          ENDDO
          ENDDO
!
! *** ADD TO THE RHS :
!
          DO INODE = 1, NNODE
            IN = NODEN(INODE)
            RHS(IN) = RHS(IN) + RHSEL(INODE)
          ENDDO
!
! *** END OF LOOP OVER THE ELEMENTS
!
        ENDDO
!
        RETURN
      END
