SUBROUTINE CGQF_F77(RULE,ORDER,X,W)
IMPLICIT NONE
INTEGER :: ORDER, RULE
REAL(8), DIMENSION(ORDER),INTENT(OUT) :: X, W
CALL CGQF(ORDER, RULE, 0.0D0, 0.0D0, -1.0D0, +1.0D0, X, W)
END SUBROUTINE CGQF_F77