MODULE decimal
  IMPLICIT NONE
  !
  ! Define precision
  !
  INTEGER, PARAMETER :: dp=KIND(1.0d0)
  REAL, PARAMETER :: precision = EPSILON(1.0_dp)
END MODULE decimal
