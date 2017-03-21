MODULE numeric_schemes_nd
USE decimal
USE tipos
IMPLICIT NONE
  abstract interface
    subroutine compute_dt(uu, CFL, dx, dt)
      REAL(kind=8), INTENT(in)  :: uu(:,:), CFL, dx 
      REAL(kind=8)              :: dt
    end subroutine
  end interface
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!11 Entropy Stable Scheme, non conservative Diffusion
subroutine Entropy_NonConservative_nd(time_scheme, Extra_Viscosity, uu, N, eqs, Tend, dx, CFL, Flux, KK, cdt, epsilon)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, percentage,eqs
  REAL(kind = dp)           :: dx, dt, limit, CFL, Tend
  REAL(kind = dp)           :: uu(:,:), tt
  REAL(kind = dp), ALLOCATABLE    :: uold(:,:), utemp(:,:), utemp2(:,:)
  procedure(compute_dt)           :: Cdt
  LOGICAL                         :: Extra_Viscosity
  REAL(kind = dp)                 :: epsilon

    interface
      function Flux(ul, ur) result(f)
          import  :: dp
          real(kind = dp), intent(in) :: ul(:), ur(:)
          real(kind = dp) :: f(SIZE(ul,1))
      end function Flux
    end interface

    interface
      function KK(ul, ur) result(k)
          import :: dp
          real(kind = dp), intent(in) :: ul(:), ur(:)
          real(kind = dp) :: k(SIZE(ul,1), SIZE(ul,1))
      end function KK
    end interface

  percentage = 0
  limit = Tend/5
  ALLOCATE(uold(N,eqs), utemp(N,eqs), utemp2(N,eqs))
  uold = 0.0_dp; utemp = 0.0_dp; utemp2 = 0.0_dp
  Print *, "Starting computing with Entropy Stable scheme (non-conservative diffusion)"
  tt = 0.0_dp
  DO WHILE (tt <= Tend)
    uold = uu
    CALL cdt(uold, CFL, dx, dt)
    !Scheme with forward Euler
    IF (time_scheme == FORWARD_EULER) THEN
      CALL update_u_NC_nd(uu, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, KK,eqs)
      !Scheme with TVD_RK2
    ELSE IF (time_scheme == TVD_RK2) THEN
      !FIRST STEP
      CALL update_u_NC_nd(utemp, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, KK,eqs)
      !SECOND STEP
      CALL update_u_NC_nd(utemp2, utemp, N, dx, dt, Extra_Viscosity, epsilon, Flux,KK,eqs)
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
end subroutine Entropy_NonConservative_nd

subroutine update_u_NC_nd(uu, uold, N, dx, dt, Extra_Viscosity, epsilon, Flux, KK,eqs)
    REAL(kind = dp)               :: epsilon, dx, dt
    REAL(kind = dp)               :: uu(:,:), uold(:,:)
    REAL(kind = dp), ALLOCATABLE  :: kl(:,:), kr(:,:)
    INTEGER                       :: N, j,eqs
    REAL(kind = dp), ALLOCATABLE  :: uleft(:), uright(:)
    LOGICAL                       :: Extra_Viscosity

    interface
      function Flux(ul, ur) result(f)
          import  :: dp
          real(kind = dp), intent(in) :: ul(:), ur(:)
          real(kind = dp) :: f(SIZE(ul,1))
      end function Flux
    end interface

    interface
      function KK(ul, ur) result(k)
          import :: dp
          real(kind = dp), intent(in) :: ul(:), ur(:)
          real(kind = dp) :: k(SIZE(ul,1), SIZE(ul,1))
      end function KK
    end interface

    ! Ghost Cells
    ALLOCATE(uleft(eqs), uright(eqs), kl(eqs,eqs), kr(eqs,eqs))
    uleft = 0.0_dp; uright = 0.0_dp; kl = 0.0_dp; kr = 0.0_dp
    uleft = uold(1,:); uright = uold(N,:)

    j = 1
    kl = KK(uleft,uold(j,:))
    kr = KK(uold(j,:),uold(j+1,:))
    uu(j,:) = uold(j,:) - dt/dx*(Flux(uold(j,:), uold(j+1,:))-Flux(uleft, uold(j,:))) +&
    dt/dx**2*(MATMUL(kr,(uold(j+1,:)-uold(j,:))) - MATMUL(kl,(uold(j,:)-uleft))) +&
    epsilon*dt/dx**2*merge(uold(j+1,:)-2*uold(j,:)+uleft,0.0_dp,Extra_Viscosity)
    DO j = 2,(N-1)
    kl = KK(uold(j-1,:),uold(j,:))
    kr = KK(uold(j,:),uold(j+1,:))
      uu(j,:) = uold(j,:) - dt/dx*(Flux(uold(j,:), uold(j+1,:))-Flux(uold(j-1,:), uold(j,:))) +&
      dt/dx**2*(MATMUL(kr,(uold(j+1,:)-uold(j,:))) - MATMUL(kl,(uold(j,:)-uold(j-1,:))))+&
      epsilon*dt/dx**2*merge(uold(j+1,:)-2*uold(j,:)+uold(j-1,:),0._dp,Extra_Viscosity)
    END DO
    j = N
    kl = KK(uold(j-1,:),uold(j,:))
    kr = KK(uold(j,:),uright)
    uu(j,:) = uold(j,:) - dt/dx*(Flux(uold(j,:), uright)-Flux(uold(j-1,:), uold(j,:))) +&
    dt/dx**2*(MATMUL(kr,(uright-uold(j,:))) - MATMUL(kl,(uold(j,:)-uold(j-1,:))))+&
    epsilon*dt/dx**2*merge(uright-2*uold(j,:)+uold(j-1,:),0.0_dp,Extra_Viscosity)
    DEALLOCATE(uleft, uright, kl, kr)
end subroutine update_u_NC_nd

END MODULE numeric_schemes_nd
