SUBROUTINE RFILLA(NAME, A, B, C, VAL)
  IMPLICIT NONE
  INTEGER A, B, C, I, J, K
  REAL VAL
  REAL NAME(A, B, C)
!
  DO K = 1, C
    DO J = 1, B
      DO I = 1, A
        NAME(I, J, K) = VAL
      ENDDO
    ENDDO
  ENDDO
!
  RETURN
!
END
