! Example 1
! Burgers Equation
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
  CALL Engquist_Osher(FORWARD_EULER, uu, N, ntime, dx, dt, flux, DiffMat)
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
END MODULE example_burger