MODULE numeric_schemes
USE decimal
USE tipos
IMPLICIT NONE
CONTAINS
subroutine Engquist_Osher(time_scheme, uu, N, ntime, dx, dt, Flux, DiffMat)
  INTEGER                   :: time_scheme
  INTEGER                   :: N
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: dt
  REAL(kind = dp)           :: mu
  INTEGER                   :: ntime
  REAL(kind = dp), ALLOCATABLE    :: uu(:), uold(:)
  REAL(kind = dp), ALLOCATABLE    :: fplus(:), fminus(:), KK(:)
  INTEGER                   :: tt, i, j
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright, fplusleft, fminusright
  REAL(kind = dp), ALLOCATABLE    :: Kleft, Kright
  REAL(kind = dp),external        :: Flux, DiffMat
  
  uleft = uu(1); uright = uu(N)
  ALLOCATE(fplus(N), fminus(N), uold(N), KK(N))
  fplus = 0.0_dp;   fminus = 0.0_dp
  uold = 0.0_dp; KK = 0.0_dp
  DO tt = 1,ntime
    uold = uu
    DO j = 1,N
      IF (uold(j) > 0) THEN
        fplus(j) = Flux(uold(j)); fminus(j) = 0.0_dp
      ELSE
        fplus(j) = 0.0_dp; fminus(j) = Flux(uold(j))
      END IF
      KK(j) = DiffMat(uold(j))
    END DO
    fplusleft = Flux(uleft); fminusright = Flux(uright)
    Kleft = DiffMat(uleft); Kright = DiffMat(uright)

    !Engquist-Osher Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      j = 1
      uu(j) = uold(j) - dt/dx * (fplus(j) + fminus(j+1) - fplusleft-fminus(j)) +&
      dt/dx**2*(KK(j+1) - 2*KK(j) + Kleft)
      DO j = 2,(N-1)
        uu(j) = uold(j) - dt/dx * (fplus(j) + fminus(j+1) - fplus(j-1)-fminus(j)) +&
        dt/dx**2*(KK(j+1) - 2*KK(j) + KK(j-1))
      END DO
      j = N
      uu(j) = uold(j) - dt/dx * (fplus(j) + fminusright - fplus(j-1)-fminus(j)) +&
      dt/dx**2*(Kright - 2*KK(j) + KK(j-1))
    END IF
  END DO
end subroutine Engquist_Osher
END MODULE numeric_schemes
