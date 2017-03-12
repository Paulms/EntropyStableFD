MODULE numeric_schemes
USE decimal
USE tipos
IMPLICIT NONE
CONTAINS
subroutine Engquist_Osher(time_scheme, uu, N, ntime, dx, dt, Flux, DiffMat)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, percentage
  REAL(kind = dp)           :: dx, dt, limit
  INTEGER                   :: ntime
  REAL(kind = dp)           :: uu(:)
  REAL(kind = dp), ALLOCATABLE    :: uold(:), utemp(:), utemp2(:), fplus(:), fminus(:), KK(:)
  INTEGER                   :: tt
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright, fplusleft, fminusright
  REAL(kind = dp),external        :: Flux, DiffMat
  REAL(kind = dp)                 :: Kleft, Kright
  percentage = 0
  limit = ntime/5
  uleft = uu(1); uright = uu(N)
  ALLOCATE(fplus(N), fminus(N), uold(N), KK(N), utemp(N), utemp2(N))
  fplus = 0.0_dp;   fminus = 0.0_dp
  uold = 0.0_dp; KK = 0.0_dp; utemp = 0.0_dp; utemp2 = 0.0_dp
  fplusleft = Flux(uleft); fminusright = Flux(uright)
  Kleft = DiffMat(uleft); Kright = DiffMat(uright)
  Print *, "Starting computing with monotone scheme"
  DO tt = 1,ntime
    uold = uu
    CALL Compute_fluxes_EO(uold, N, fplus, fminus, KK, Flux, DiffMat)
    !Engquist-Osher Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_EO(uu, uold, N, dx, dt, fplus, fminus, fplusleft, fminusright, KK, Kleft, Kright)
    !Engquist-Osher Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_EO(utemp, uold, N, dx, dt, fplus, fminus, fplusleft, fminusright, KK, Kleft, Kright)
      !SECOND STEP
      CALL Compute_fluxes_EO(utemp, N, fplus, fminus, KK, Flux, DiffMat)
      CALL update_u_EO(utemp2, utemp, N, dx, dt, fplus, fminus, fplusleft, fminusright, KK, Kleft, Kright)
      uu = 0.5*(uold + utemp2)
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + ntime/5
      print *, percentage, "% completed"
    END IF
  END DO
  print *, "completed..."
  DEALLOCATE(fplus, fminus, uold, KK, utemp, utemp2)
end subroutine Engquist_Osher

subroutine Compute_fluxes_EO(uold, N, fplus, fminus, KK, Flux, DiffMat)
  REAL(kind = dp), intent(in)     :: uold(:)
  REAL(kind = dp)                 :: fplus(:), fminus(:), KK(:)
  INTEGER                         :: j, N
  REAL(kind = dp),external        :: Flux, DiffMat
  DO j = 1,N
    IF (uold(j) > 0.0_dp) THEN
      fplus(j) = Flux(uold(j)); fminus(j) = 0.0_dp
    ELSE
      fplus(j) = 0.0_dp; fminus(j) = Flux(uold(j))
    END IF
    KK(j) = DiffMat(uold(j))
  END DO
end subroutine Compute_fluxes_EO
subroutine update_u_EO(uu, uold, N, dx, dt, fplus, fminus, fplusleft, fminusright, KK, Kleft, Kright)
    REAL(kind = dp)           :: dx, dt
    REAL(kind = dp)           :: uu(:), uold(:), fplus(:), fminus(:), KK(:)
    INTEGER                   :: N, j
    REAL(kind = dp)           :: fplusleft, fminusright
    REAL(kind = dp)           :: Kleft, Kright
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
end subroutine update_u_EO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Entropy Stable Scheme, conservative diffusion
subroutine Entropy_Conservative(time_scheme, Extra_Viscosity, uu, N, ntime, dx, dt, Flux, DiffMat, epsilon)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, percentage
  REAL(kind = dp)           :: dx, dt, limit
  INTEGER                   :: ntime
  REAL(kind = dp)           :: uu(:)
  INTEGER                   :: tt, j
  REAL(kind = dp)           :: uleft, uright
  REAL(kind = dp), external        :: Flux, DiffMat
  REAL(kind = dp), ALLOCATABLE    :: uold(:), KK(:), utemp(:), utemp2(:)
  REAL(kind = dp)                 :: Kleft, Kright
  LOGICAL                         :: Extra_Viscosity
  REAL(kind = dp)                 :: epsilon
  percentage = 0
  limit = ntime/5
  uleft = uu(1); uright = uu(N)
  Kleft = DiffMat(uleft); Kright = DiffMat(uright)
  ALLOCATE(uold(N), KK(N), utemp(N), utemp2(N))
  uold = 0.0_dp; KK = 0.0_dp
  utemp = 0.0_dp; utemp2 = 0.0_dp
  Print *, "Starting computing with Entropy Stable scheme (conservative diffusion)"
  DO tt = 1,ntime
    uold = uu
    DO j = 1,N
      KK(j) = DiffMat(uold(j))
    END DO
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_EC(uu, uold, N, dx, dt, KK, Kleft, Kright, uleft, uright, Extra_Viscosity, epsilon, Flux)
    !Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_EC(utemp, uold, N, dx, dt, KK, Kleft, Kright, uleft, uright, Extra_Viscosity, epsilon, Flux)
      !SECOND STEP
      DO j = 1,N
        KK(j) = DiffMat(utemp(j))
      END DO
      CALL update_u_EC(utemp2, utemp, N, dx, dt, KK, Kleft, Kright, uleft, uright, Extra_Viscosity, epsilon, Flux)
      uu = 0.5*(uold + utemp2)
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + ntime/5
      print *, percentage, "% completed"
    END IF
  END DO
  print *, "completed..."
  DEALLOCATE(uold, KK, utemp, utemp2)
end subroutine Entropy_Conservative

subroutine update_u_EC(uu, uold, N, dx, dt, KK, Kleft, Kright, uleft, uright, Extra_Viscosity, epsilon, Flux)
    REAL(kind = dp)           :: epsilon, dx, dt
    REAL(kind = dp)           :: uu(:), uold(:), KK(:)
    INTEGER                   :: N, j
    REAL(kind = dp)           :: uleft, uright
    REAL(kind = dp)           :: Kleft, Kright
    LOGICAL                   :: Extra_Viscosity
    REAL(kind = dp), external :: Flux
    j = 1
    uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uleft, uold(j))) +&
    dt/dx**2*(KK(j+1) - 2*KK(j) + Kleft) +&
    epsilon*dt/dx**2*merge(uold(j+1)-2*uold(j)+uleft,0.0_dp,Extra_Viscosity)
    DO j = 2,(N-1)
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uold(j-1), uold(j))) +&
      dt/dx**2*(KK(j+1) - 2*KK(j) + KK(j-1))+&
    epsilon*dt/dx**2*merge(uold(j+1)-2*uold(j)+uold(j-1),0._dp,Extra_Viscosity)
    END DO
    j = N
    uu(j) = uold(j) - dt/dx * (Flux(uold(j), uright)-Flux(uold(j-1), uold(j))) +&
    dt/dx**2*(Kright - 2*KK(j) + KK(j-1))+&
    epsilon*dt/dx**2*merge(uright-2*uold(j)+uold(j-1),0.0_dp,Extra_Viscosity)
end subroutine update_u_EC
!!!!!!!!!!!!!!!!!!!!!!!!!11 Entropy Stable Scheme, non conservative Diffusion
subroutine Entropy_NonConservative(time_scheme, Extra_Viscosity, uu, N, ntime, dx, dt, Flux, KK, epsilon)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, percentage
  REAL(kind = dp)           :: dx, dt, limit
  INTEGER                   :: ntime
  REAL(kind = dp)           :: uu(:)
  REAL(kind = dp), ALLOCATABLE    :: uold(:), utemp(:), utemp2(:)
  INTEGER                   :: tt
  REAL(kind = dp)           :: uleft, uright
  REAL(kind = dp),external        :: Flux, KK
  LOGICAL                         :: Extra_Viscosity
  REAL(kind = dp)                 :: epsilon
  percentage = 0
  limit = ntime/5
  uleft = uu(1); uright = uu(N)
  ALLOCATE(uold(N), utemp(N), utemp2(N))
  uold = 0.0_dp; utemp = 0.0_dp; utemp2 = 0.0_dp
  Print *, "Starting computing with Entropy Stable scheme (non-conservative diffusion)"
  DO tt = 1,ntime
    uold = uu
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_NC(uu, uold, N, dx, dt, KK, uleft, uright, Extra_Viscosity, epsilon, Flux)
      !Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_NC(utemp, uold, N, dx, dt, KK, uleft, uright, Extra_Viscosity, epsilon, Flux)
      !SECOND STEP
      CALL update_u_NC(utemp2, utemp, N, dx, dt, KK, uleft, uright, Extra_Viscosity, epsilon, Flux)
      uu = 0.5*(uold + utemp2)
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + ntime/5
      print *, percentage, "% completed"
    END IF
  END DO
  print *, "completed..."
  DEALLOCATE(uold, utemp, utemp2)
end subroutine Entropy_NonConservative

subroutine update_u_NC(uu, uold, N, dx, dt, KK, uleft, uright, Extra_Viscosity, epsilon, Flux)
    REAL(kind = dp)           :: epsilon, dx, dt
    REAL(kind = dp)           :: uu(:), uold(:)
    INTEGER                   :: N, j
    REAL(kind = dp)           :: uleft, uright
    LOGICAL                   :: Extra_Viscosity
    REAL(kind = dp), external :: Flux, KK
    j = 1
    uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uleft, uold(j))) +&
    dt/dx**2*(KK(uold(j),uold(j+1))*(uold(j+1)-uold(j)) - KK(uleft,uold(j))*(uold(j)-uleft)) +&
    epsilon*dt/dx**2*merge(uold(j+1)-2*uold(j)+uleft,0.0_dp,Extra_Viscosity)
    DO j = 2,(N-1)
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uold(j-1), uold(j))) +&
      dt/dx**2*(KK(uold(j),uold(j+1))*(uold(j+1)-uold(j)) - KK(uold(j-1),uold(j))*(uold(j)-uold(j-1)))+&
      epsilon*dt/dx**2*merge(uold(j+1)-2*uold(j)+uold(j-1),0._dp,Extra_Viscosity)
    END DO
    j = N
    uu(j) = uold(j) - dt/dx * (Flux(uold(j), uright)-Flux(uold(j-1), uold(j))) +&
    dt/dx**2*(KK(uold(j),uright)*(uright-uold(j)) - KK(uold(j-1),uold(j))*(uold(j)-uold(j-1)))+&
    epsilon*dt/dx**2*merge(uright-2*uold(j)+uold(j-1),0.0_dp,Extra_Viscosity)
end subroutine update_u_NC
END MODULE numeric_schemes
