MODULE numeric_schemes
USE decimal
USE tipos
IMPLICIT NONE
  abstract interface
    subroutine compute_dt(uu, CFL, dx, dt)
      REAL(kind=8), INTENT(in)  :: uu(:), CFL, dx 
      REAL(kind=8)              :: dt
    end subroutine
  end interface
CONTAINS
subroutine Engquist_Osher(time_scheme, uu, N, Tend, dx, CFL, Flux, DiffMat, Cdt)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, percentage
  REAL(kind = dp)           :: dx, dt, CFL, limit, Tend
  REAL(kind = dp)           :: uu(:), tt
  REAL(kind = dp), ALLOCATABLE    :: uold(:), utemp(:), utemp2(:)
  REAL(kind = dp),external        :: Flux, DiffMat
  procedure(compute_dt)           :: Cdt
  percentage = 0
  limit = Tend/5
  ALLOCATE(uold(N), utemp(N), utemp2(N))
  uold = 0.0_dp; utemp = 0.0_dp; utemp2 = 0.0_dp
  Print *, "Starting computing with monotone scheme"
  tt = 0.0_dp
  DO WHILE (tt <= Tend)
    uold = uu
    CALL cdt(uold, CFL, dx, dt)
    !Engquist-Osher Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_EO(uu, uold, N, dx, dt, Flux, DiffMat)
    !Engquist-Osher Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_EO(utemp, uold, N, dx, dt, Flux, DiffMat)
      !SECOND STEP
      CALL update_u_EO(utemp2, utemp, N, dx, dt, Flux, DiffMat)
      uu = 0.5*(uold + utemp2)
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + Tend/5
      print *, percentage, "% completed"
    END IF
    tt = tt + dt
  END DO
  print *, "completed..."
  DEALLOCATE(uold, utemp, utemp2)
end subroutine Engquist_Osher

subroutine update_u_EO(uu, uold, N, dx, dt, Flux, DiffMat)
    REAL(kind = dp), INTENT(IN)   :: dx, dt
    REAL(kind = dp)               :: uu(:), uold(:)
    REAL(kind = dp), ALLOCATABLE  :: fplus(:), fminus(:), KK(:)
    INTEGER                       :: N, j
    REAL(kind = dp)               :: fplusleft, fminusright
    REAL(kind = dp)               :: Kleft, Kright, uleft, uright
    REAL(kind = dp),external      :: Flux, DiffMat
    !Allocate vectores
    ALLOCATE(fplus(N), fminus(N), KK(N))
    fplus = 0.0_dp;   fminus = 0.0_dp; KK = 0.0_dp
    
    DO j = 1,N
      IF (uold(j) > 0.0_dp) THEN
        fplus(j) = Flux(uold(j)); fminus(j) = 0.0_dp
      ELSE
        fplus(j) = 0.0_dp; fminus(j) = Flux(uold(j))
      END IF
      KK(j) = DiffMat(uold(j))
    END DO
    uleft = uold(1); uright = uold(N)
    fplusleft = merge(Flux(uleft),0.0_dp,uleft > 0.0_dp)
    fminusright = merge(0.0_dp,Flux(uright),uright > 0.0_dp)
    Kleft = DiffMat(uleft); Kright = DiffMat(uright)

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
    DEALLOCATE(fplus, fminus, KK)
end subroutine update_u_EO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Entropy Stable Scheme, conservative diffusion
subroutine Entropy_Conservative(time_scheme, Extra_Viscosity, uu, N, Tend, dx, CFL, Flux, DiffMat, Cdt, epsilon)
  INTEGER                       :: time_scheme
  INTEGER                       :: N, percentage
  REAL(kind = dp)               :: dx, dt, limit, Tend, CFL
  REAL(kind = dp)               :: uu(:), tt
  INTEGER                       :: j
  REAL(kind = dp), external     :: Flux, DiffMat
  procedure(compute_dt)         :: Cdt
  REAL(kind = dp), ALLOCATABLE  :: uold(:), utemp(:), utemp2(:)
  LOGICAL                       :: Extra_Viscosity
  REAL(kind = dp)               :: epsilon
  percentage = 0
  limit = Tend/5
  ALLOCATE(uold(N), utemp(N), utemp2(N))
  uold = 0.0_dp;
  utemp = 0.0_dp; utemp2 = 0.0_dp
  Print *, "Starting computing with Entropy Stable scheme (conservative diffusion)"
  tt = 0.0_dp
  DO WHILE (tt <= Tend)
    uold = uu
    CALL cdt(uold, CFL, dx, dt)
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_EC(uu, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, DiffMat)
    !Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_EC(utemp, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, DiffMat)
      !SECOND STEP
      CALL update_u_EC(utemp2, utemp, N, dx, dt, Extra_Viscosity, epsilon, Flux, DiffMat)
      uu = 0.5*(uold + utemp2)
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + Tend/5
      print *, percentage, "% completed"
    END IF
    tt = tt + dt
  END DO
  print *, "completed..."
  DEALLOCATE(uold, utemp, utemp2)
end subroutine Entropy_Conservative

subroutine update_u_EC(uu, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, DiffMat)
    REAL(kind = dp)               :: epsilon, dx, dt
    REAL(kind = dp)               :: uu(:), uold(:)
    REAL(kind = dp), ALLOCATABLE  :: KK(:)
    INTEGER                       :: N, j
    REAL(kind = dp)               :: uleft, uright
    REAL(kind = dp)               :: Kleft, Kright
    LOGICAL                       :: Extra_Viscosity
    REAL(kind = dp), external     :: Flux, DiffMat
    !Allocate vectors
    ALLOCATE(KK(N))
    KK = 0.0_dp
    ! Update ghost cells
    uleft = uold(1); uright = uold(N)
    Kleft = DiffMat(uleft); Kright = DiffMat(uright)

    ! Compute diffusion flux
    DO j = 1,N
      KK(j) = DiffMat(uold(j))
    END DO

    j = 1
    uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uleft, uold(j))) +&
    dt/dx**2*(KK(j+1) - 2*KK(j) + Kleft) +&
    epsilon*dt/dx**2*merge(uold(j+1)-2*uold(j)+uleft,0.0_dp,Extra_Viscosity)
    DO j = 2,(N-1)
      uu(j) = uold(j) - dt/dx * (Flux(uold(j), uold(j+1))-Flux(uold(j-1), uold(j))) +&
      dt/dx**2*(KK(j+1) - 2*KK(j) + KK(j-1))+&
    epsilon*dt/dx**2*merge(uold(j+1)-2*uold(j)+uold(j-1),0.0_dp,Extra_Viscosity)
    END DO
    j = N
    uu(j) = uold(j) - dt/dx * (Flux(uold(j), uright)-Flux(uold(j-1), uold(j))) +&
    dt/dx**2*(Kright - 2*KK(j) + KK(j-1))+&
    epsilon*dt/dx**2*merge(uright-2*uold(j)+uold(j-1),0.0_dp,Extra_Viscosity)

    DEALLOCATE(KK)
end subroutine update_u_EC
!!!!!!!!!!!!!!!!!!!!!!!!!11 Entropy Stable Scheme, non conservative Diffusion
subroutine Entropy_NonConservative(time_scheme, Extra_Viscosity, uu, N, Tend, dx, CFL, Flux, KK, cdt, epsilon)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, percentage
  REAL(kind = dp)           :: dx, dt, limit, CFL, Tend
  REAL(kind = dp)           :: uu(:), tt
  REAL(kind = dp), ALLOCATABLE    :: uold(:), utemp(:), utemp2(:)
  REAL(kind = dp),external        :: Flux, KK
  procedure(compute_dt)           :: Cdt
  LOGICAL                         :: Extra_Viscosity
  REAL(kind = dp)                 :: epsilon
  percentage = 0
  limit = Tend/5
  ALLOCATE(uold(N), utemp(N), utemp2(N))
  uold = 0.0_dp; utemp = 0.0_dp; utemp2 = 0.0_dp
  Print *, "Starting computing with Entropy Stable scheme (non-conservative diffusion)"
  tt = 0.0_dp
  DO WHILE (tt <= Tend)
    uold = uu
    CALL cdt(uold, CFL, dx, dt)
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_NC(uu, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, KK)
      !Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_NC(utemp, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, KK)
      !SECOND STEP
      CALL update_u_NC(utemp2, utemp, N, dx, dt, Extra_Viscosity, epsilon, Flux,KK)
      uu = 0.5*(uold + utemp2)
    END IF
    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 20
      limit = limit + Tend/5
      print *, percentage, "% completed"
    END IF
    tt = tt + dt
  END DO
  print *, "completed..."
  DEALLOCATE(uold, utemp, utemp2)
end subroutine Entropy_NonConservative

subroutine update_u_NC(uu, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, KK)
    REAL(kind = dp)           :: epsilon, dx, dt
    REAL(kind = dp)           :: uu(:), uold(:)
    INTEGER                   :: N, j
    REAL(kind = dp)           :: uleft, uright
    LOGICAL                   :: Extra_Viscosity
    REAL(kind = dp), external :: Flux, KK

    ! Ghost Cells
    uleft = uold(1); uright = uold(N)

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

!Auxiliar function to compute errors
FUNCTION cumpute_errors(ref, M, uu, N) RESULT(error)
  ! Compute errors in L1 norm
  INTEGER, INTENT(IN)          :: N, M
  REAL(kind = dp), INTENT(IN)  :: ref(:), uu(:)
  REAL(kind=dp), ALLOCATABLE   :: uexact(:)
  REAL(kind = dp)              :: error
  INTEGER                      :: i, j, R
  ! Compute exact values
  ALLOCATE(uexact(N))
  uexact = 0.0_dp
  R = NINT(1.0*M/N)
  DO i = 1, N
    uexact(i) = 1.0_dp/R*sum(ref((R*(i-1)+1):(R*i)))
  END DO
  error = 1.0_dp/N*sum(abs(uu - uexact))
  print *, "Error: ", error
  DEALLOCATE(uexact)
END FUNCTION cumpute_errors
END MODULE numeric_schemes
