! Test Problem 1
! Burgers Equation
! u_t + (u^2/2)_x=K(u)_xx   (x,t) in [-2,2]x[0,T]
! where
! K(u) = mu*u^2
MODULE example_burger
USE decimal
USE plot
USE tipos
USE numeric_schemes
IMPLICIT NONE
REAL(kind = dp), PARAMETER :: mu = 0.01
CONTAINS
subroutine burger_runexample(initial_condition)
  INTEGER                   :: initial_condition
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL
  REAL(kind = dp)           :: dt
  INTEGER                   :: ntime
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:), uu(:), uold(:)
  REAL(kind = dp), ALLOCATABLE    :: fplus(:), fminus(:), KK(:)
  INTEGER                         :: tt, i, j
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright, fplusleft, fminusright
  REAL(kind = dp), ALLOCATABLE    :: Kleft, Kright
  CHARACTER(LEN=32)               :: name           ! File name to save plot data

  ! Initialize variables
  Tend = 0.5_dp      ! Final Time
  N = 400          ! Number of nodes
  CFL = 0.9_dp
  dx = 4.0_dp/(N-1)
  dt = CFL/(1.0_dp*(1/dx+4*mu/dx**2))
  ntime = floor(tend/dt)  !Number of time steps

    ! Allocate memory
  ALLOCATE(xx(N), uu(N))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx-dx-2,i=1, N)/)      ! Location of grid points

  ! Initial Conditions
  IF (initial_condition == 1) THEN
    DO j = 1, N
      IF (xx(j) > -1 .AND. xx(j) < 1) THEN
        uu(j) = (1 - xx(j)**2)**2
      ELSE
        uu(j) = 0.0_dp
      END IF
    END DO
  elseif (initial_condition == 2) then    
    DO j = 1, N
      IF (xx(j) > -0.5 .AND. xx(j) < 0.5) THEN
        uu(j) = 1
      ELSE
        uu(j) = 0.0_dp
      END IF
    END DO
  else
    uu = 0.0_dp
  END IF
  !Save initial condition
  ALLOCATE(uinit(N))
  uinit = 0.0_dp;   
  uinit = uu

  !Engquist-Osher Scheme with forward Euler
  name = 'output'
  !CALL Engquist_Osher(FORWARD_EULER, uu, N, ntime, dx, dt, flux, DiffMat)  !MS
  !CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, ntime, dx, dt, fluxEC, DiffMat, 0.0) !ESC
  !CALL Entropy_NonConservative(FORWARD_EULER, .FALSE., uu, N, ntime, dx, dt, fluxEC, KKN, 0.0) !ESNC
  !CALL Entropy_Conservative(FORWARD_EULER, .TRUE., uu, N, ntime, dx, dt, fluxEC, DiffMat, dx*0.1)  !ESC-alpha
  CALL Entropy_NonConservative(FORWARD_EULER, .TRUE., uu, N, ntime, dx, dt, fluxEC, KKN, dx*0.1) !ESNC-alpha
  CALL plot_results(uu, uinit, xx, name)
end subroutine burger_runexample

FUNCTION flux(uu) RESULT(ff)
  !Flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: ff
  ff = 0.5*uu**2
END FUNCTION flux

FUNCTION DiffMat(uu) RESULT(kk)
  !Conservative diffusion in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: kk
  kk = mu*uu**2
END FUNCTION DiffMat

FUNCTION fluxEC(ul, ur) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: ul
  REAL(kind = dp), INTENT(IN)  :: ur
  REAL(kind = dp)              :: ff
  ff = 1/6.0*(ur**2+ur*ul+ul**2)
END FUNCTION fluxEC

FUNCTION KKN(ul, ur) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: ul
  REAL(kind = dp), INTENT(IN)  :: ur
  REAL(kind = dp)              :: kk
  If (ul == 0 .AND. ur ==0) THEN
    kk = 0
  ELSE
    kk = mu*4/3*(ul**2+ul*ur+ur**2)/(ul+ur)
  END IF
END FUNCTION KKN

END MODULE example_burger