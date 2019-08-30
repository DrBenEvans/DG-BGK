SUBROUTINE IFILLA(MA, NA, LA, KA, val)
  IMPLICIT NONE
  INTEGER NA, LA, KA, I, J, K
  INTEGER MA(NA, LA, KA), val
!
  DO K = 1, KA
    DO J = 1, LA
      DO I = 1, NA
        MA(I, J, K) = val
      ENDDO
    ENDDO
  ENDDO
!
  RETURN
END
