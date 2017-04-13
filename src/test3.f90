! Test Problem 2
! u_t + (u^2/2)_x=K(u)_xx   (x,t) in [-pi/2,pi]x[0,T]
MODULE test3
USE decimal
USE plot
USE tipos
USE numeric_schemes
IMPLICIT NONE
PUBLIC test3_run
PRIVATE
CONTAINS
subroutine test3_run(run_ref, run_err)
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N, M
  REAL(kind = dp)           :: dx, error
  REAL(kind = dp)           :: CFL
  INTEGER                   :: ntests
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  INTEGER, ALLOCATABLE            ::steps(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:), uu(:), reference(:,:)
  INTEGER                         :: i
  CHARACTER(LEN=32)               :: name           ! File name to save plot data
  REAL(kind = dp), ALLOCATABLE    :: results(:,:)
  CHARACTER(LEN=8), ALLOCATABLE  :: names(:)
  LOGICAL                         :: run_ref, run_err
  ! Zero variables
  Tend = 0.0_dp; error = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 1.0_dp      ! Final Time
  CFL = 0.9_dp
  M = 8000

  !Save reference solution (Change binary flag if needed)
  IF (run_ref) THEN
    name = 'test_3_reference'
    N = M          ! Number of nodes
    CALL setup_problem(dx, N, xx, uu, uinit)
    ntests = 1
    ALLOCATE(results(N, ntests+1), names(ntests+1))
    names = ['x    ', 'REF  ']
    results(:,1) = xx
    CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp, ZERO_FLUX) !ESC-0.2
    results(:,2) = uu
    CALL save_matrix(results, names, name, 0)
    DEALLOCATE(results, uu, names, uinit, xx)
  END IF

  !Compute errors (Change binary flag if needed)
  IF (run_err) THEN
    name = 'test_3_reference'
    ALLOCATE(steps(4))
    steps = [200,400,800,1600]
    !Read reference solution
    CALL read_matrix(name , reference, M, 2)
    ! Prepare to save matrix
    name = 'test_3_errors'
    ntests = 3

    ALLOCATE(results(ntests+1, 4), names(ntests+1))
    results = 0.0_dp
    names = ['N       ', 'MS      ', 'ESC     ', 'ESC-0.3 ']
    results(1,:) = steps
    DO i = 1, 4
      N = steps(i)
      print *, "Starting numerical tests with N = ", N
      CALL setup_problem(dx, N, xx, uu, uinit)
      CALL Engquist_Osher(FORWARD_EULER, uu, N, Tend, dx, CFL, flux, DiffMat, Cdt, ZERO_FLUX)  !MS
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(2, i) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp, ZERO_FLUX) !ESC
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(3, i) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_Conservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.3_dp*dx, ZERO_FLUX) !ESNC
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(4, i) = error
      uu = uinit; error = 0.0_dp
      DEALLOCATE(uu, uinit, xx)
    END DO
    CALL save_matrix(results, names, name, 1)
    DEALLOCATE(results, names)
    STOP
  END IF

  !Run numerical schemes
  N = 100          ! Number of nodes
  CALL setup_problem(dx, N, xx, uu, uinit)
  name = 'test_3_100'
  ntests = 4
  ALLOCATE(results(N, ntests+1), names(ntests+1))
  names = ['X       ', 'u0      ','MS      ', 'ESC     ', 'ESC-0.3 ']
  results(:,1) = xx
  results(:,2) = uinit
  CALL Engquist_Osher(FORWARD_EULER, uu, N, Tend, dx, CFL, flux, DiffMat, Cdt, ZERO_FLUX)  !MS
  results(:,3) = uu
  uu = uinit
  CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp, ZERO_FLUX) !ESC-0.2
  results(:,4) = uu
  uu = uinit
  CALL Entropy_Conservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.3_dp*dx, ZERO_FLUX) !ESNC-0.2
  results(:,5) = uu
  CALL save_matrix(results, names, name, 0)
  !CALL plot_results(uu, uinit, xx, name)
  ! Clean memory
  DEALLOCATE(results, uu, names, uinit, xx)
end subroutine test3_run

SUBROUTINE setup_problem(dx, N, xx, uu, uinit)
  INTEGER, INTENT(IN)             :: N
  REAL(kind = dp)                 :: dx
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:), uu(:)
  INTEGER                         :: i, j
  dx = (3.0/2.0_dp*pi)/N

    ! Allocate memory
  ALLOCATE(xx(N), uu(N))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx+dx/2-pi/2.0_dp,i=0, N-1)/)      ! Location of grid points

  ! Initial Conditions
  DO j = 1, N
    uu(j) = sin(xx(j))
  END DO
  !Save initial condition
  ALLOCATE(uinit(N))
  uinit = 0.0_dp;   
  uinit = uu
END SUBROUTINE   

FUNCTION flux(uu) RESULT(ff)
  !Flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: ff
  if (abs(uu) < precision) THEN
    ff = 0.0_dp
  ELSE
    ff = uu**2/2.0_dp
  END IF
END FUNCTION flux

!k(u)
FUNCTION DiffFlux(uu) RESULT(kk)
  !Conservative diffusion flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu(:)
  REAL(kind = dp), ALLOCATABLE :: kk(:)
  INTEGER                      :: N, j
  N=SIZE(uu,1)
  ALLOCATE(kk(N))
  kk = 0.0_dp
  DO j = 1,N
    if (abs(uu(j)) < precision) THEN
      kk(j) = 0.0_dp
    ELSE
      kk(j) = uu(j)
    END iF
  END DO
END FUNCTION DiffFlux

!K(x)
FUNCTION DiffMat(uu) RESULT(kk)
  !Conservative diffusion in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: kk
  kk = 0.5_dp*max(uu, 0.0_dp)**2
END FUNCTION DiffMat

FUNCTION fluxEC(ul, ur) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: ul
  REAL(kind = dp), INTENT(IN)  :: ur
  REAL(kind = dp)              :: ff
  If (abs(ul) <= precision .AND. abs(ur) <= precision) THEN
    ff = 0.0_dp
  ELSE
    ff = 1/6.0_dp*(ur**2+ur*ul+ul**2)
  END IF
END FUNCTION fluxEC

SUBROUTINE Cdt(uold, CFL, dx, dt)
  ! Update dt based on CFL condition
  REAL(kind = dp), INTENT(IN)  :: uold(:)
  REAL(kind = dp), INTENT(IN)  :: CFL
  REAL(kind = dp), INTENT(IN)  :: dx
  REAL(kind = dp)              :: dt
  dt = 0.0_dp
  dt = CFL/(1.0_dp/dx*MAXVAL(ABS(uold))+1.0_dp/dx**2*2*MAXVAL(ABS(DiffFlux(uold))))
END SUBROUTINE
END MODULE test3