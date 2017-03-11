MODULE numeric_schemes
USE decimal
USE tipos
IMPLICIT NONE
CONTAINS
subroutine Engquist_Osher(time_scheme, uu, N, ntime, dx, dt, Flux, DiffMat)
  INTEGER                   :: time_scheme
  INTEGER                   :: N
  REAL(kind = dp)           :: dx, dt
  INTEGER                   :: ntime
  REAL(kind = dp), ALLOCATABLE    :: uu(:), uold(:)
  REAL(kind = dp), ALLOCATABLE    :: fplus(:), fminus(:), KK(:)
  INTEGER                   :: tt, i, j
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright, fplusleft, fminusright
  REAL(kind = dp),external        :: Flux, DiffMat
  REAL(kind = dp), ALLOCATABLE    :: Kleft, Kright
  
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

subroutine Entropy_Conservative(time_scheme, uu, N, ntime, dx, dt, Flux, DiffMat)
  INTEGER                   :: time_scheme
  INTEGER                   :: N
  REAL(kind = dp)           :: dx, dt
  INTEGER                   :: ntime
  REAL(kind = dp), ALLOCATABLE    :: uu(:), uold(:)
  INTEGER                   :: tt, i, j
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright
  REAL(kind = dp),external        :: Flux, DiffMat
  REAL(kind = dp), ALLOCATABLE    :: KK(:)
  REAL(kind = dp), ALLOCATABLE    :: Kleft, Kright
  
  uleft = uu(1); uright = uu(N)
  Kleft = DiffMat(uleft); Kright = DiffMat(uright)
  ALLOCATE(uold(N), KK(N))
  uold = 0.0_dp; KK = 0.0_dp
  DO tt = 1,ntime
    uold = uu
    DO j = 1,N
      KK(j) = DiffMat(uold(j))
    END DO
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      j = 1
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uleft, uold(j))) +&
      dt/dx**2*(KK(j+1) - 2*KK(j) + Kleft)
      DO j = 2,(N-1)
        uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uold(j-1), uold(j))) +&
        dt/dx**2*(KK(j+1) - 2*KK(j) + KK(j-1))
      END DO
      j = N
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uright)-Flux(uold(j-1), uold(j))) +&
      dt/dx**2*(Kright - 2*KK(j) + KK(j-1))
    END IF
  END DO
end subroutine Entropy_Conservative

subroutine Entropy_NonConservative(time_scheme, uu, N, ntime, dx, dt, Flux, KK)
  INTEGER                   :: time_scheme
  INTEGER                   :: N
  REAL(kind = dp)           :: dx, dt
  INTEGER                   :: ntime
  REAL(kind = dp), ALLOCATABLE    :: uu(:), uold(:)
  INTEGER                   :: tt, i, j
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright
  REAL(kind = dp),external        :: Flux, KK
  
  uleft = uu(1); uright = uu(N)
  ALLOCATE(uold(N))
  uold = 0.0_dp
  DO tt = 1,ntime
    uold = uu
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      j = 1
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uleft, uold(j))) +&
      dt/dx**2*(KK(uold(j),uold(j+1))*(uold(j+1)-uold(j)) - KK(uleft,uold(j))*(uold(j)-uleft))
      DO j = 2,(N-1)
        uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uold(j-1), uold(j))) +&
        dt/dx**2*(KK(uold(j),uold(j+1))*(uold(j+1)-uold(j)) - KK(uold(j-1),uold(j))*(uold(j)-uold(j-1)))
      END DO
      j = N
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uright)-Flux(uold(j-1), uold(j))) +&
      dt/dx**2*(KK(uold(j),uright)*(uright-uold(j)) - KK(uold(j-1),uold(j))*(uold(j)-uold(j-1)))
    END IF
  END DO
end subroutine Entropy_NonConservative
END MODULE numeric_schemes
