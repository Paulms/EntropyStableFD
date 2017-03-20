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
PUBLIC burger_runexample
PRIVATE
REAL(kind = dp), PARAMETER :: mu = 0.01_dp
CONTAINS
subroutine burger_runexample(initial_condition, run_ref, run_err)
  INTEGER                   :: initial_condition
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
  Tend = 0.5_dp      ! Final Time
  CFL = 0.9_dp
  M = 16000

  !Save reference solution (Change binary flag if needed)
  IF (run_ref) THEN
    IF (initial_condition == 1) THEN
      name = 'test_1_reference'
    ELSE
      name = 'test_1b_reference'
    END IF
    N = M          ! Number of nodes
    CALL setup_problem(dx, N, xx, uu, uinit, initial_condition)
    ntests = 1
    ALLOCATE(results(N, ntests+1), names(ntests+1))
    names = ['x       ', 'REF     ']
    results(:,1) = xx
    CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp) !ESC
    results(:,2) = uu
    CALL save_matrix(results, names, name, 0)
    DEALLOCATE(results, uu, names, uinit, xx)
  END IF

  !Compute errors (Change binary flag if needed)
  IF (run_err) THEN
    IF (initial_condition == 1) THEN
      name = 'test_1_reference'
    ELSE
      name = 'test_1b_reference'
    END IF
    ALLOCATE(steps(5))
    steps = [200,400,800,1600,3200]
    !Read reference solution
    CALL read_matrix(name , reference, M, 2)
    ! Prepare to save matrix
    IF (initial_condition == 1) THEN
      name = 'test_1_errors'
    ELSE
      name = 'test_1b_errors'
    END IF
    ntests = 5

    ALLOCATE(results(ntests+1, 5), names(ntests+1))
    results = 0.0_dp
    names = ['N       ', 'MS      ', 'ESC     ', 'ESNC    ', 'ESC2    ','ESNC2   ']
    results(1,:) = steps
    DO i = 1, 5
      N = steps(i)
      print *, "Starting numerical tests with N = ", N
      CALL setup_problem(dx, N, xx, uu, uinit, initial_condition)
      CALL Engquist_Osher(FORWARD_EULER, uu, N, Tend, dx, CFL, flux, DiffMat, Cdt)  !MS
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(2,i) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp) !ESC
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(3,i) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_NonConservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.0_dp) !ESNC
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(4,i) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_Conservative(TVD_RK2, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp) !ESC2
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(5,i) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_NonConservative(TVD_RK2, .FALSE., uu, N, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.0_dp) !ESNC2
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(6,i) = error
      DEALLOCATE(uu, uinit, xx)
    END DO
    CALL save_matrix(results, names, name,1)
    DEALLOCATE(results, names)
    STOP
  END IF

  !Run numerical schemes
  N = 400          ! Number of nodes
  CALL setup_problem(dx, N, xx, uu, uinit, initial_condition)
    IF (initial_condition == 1) THEN
    name = 'test_1_400'
  ELSE
    name = 'test_1b_400'
  END IF
  ntests = 5
  ALLOCATE(results(N, ntests+1), names(ntests+1))
  names = ['X       ', 'MS      ', 'ESC     ', 'ESNC    ', 'ESC2    ','ESNC2   ']
  results(:,1) = xx
  CALL Engquist_Osher(FORWARD_EULER, uu, N, Tend, dx, CFL, flux, DiffMat, Cdt)  !MS
  results(:,2) = uu
  uu = uinit
  CALL Entropy_Conservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp) !ESC
  results(:,3) = uu
  uu = uinit
  CALL Entropy_NonConservative(FORWARD_EULER, .FALSE., uu, N, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.0_dp) !ESNC
  results(:,4) = uu
  uu = uinit
  CALL Entropy_Conservative(TVD_RK2, .FALSE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.0_dp) !ESC2
  results(:,5) = uu
  uu = uinit
  CALL Entropy_NonConservative(TVD_RK2, .FALSE., uu, N, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.0_dp) !ESNC2
  results(:,6) = uu
  CALL save_matrix(results, names, name,0)
  !CALL plot_results(uu, uinit, xx, name)
  ! Clean memory
  DEALLOCATE(results, uu, names, uinit, xx)
end subroutine burger_runexample

SUBROUTINE setup_problem(dx, N, xx, uu, uinit, initial_condition)
  INTEGER, INTENT(IN)             :: initial_condition
  INTEGER, INTENT(IN)             :: N
  REAL(kind = dp)                 :: dx
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:), uu(:)
  INTEGER                         :: i, j
  dx = 4.0_dp/N

    ! Allocate memory
  ALLOCATE(xx(N), uu(N))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx+dx/2-2,i=0, N-1)/)      ! Location of grid points

  ! Initial Conditions
  IF (initial_condition == 1) THEN
    DO j = 1, N
      IF (xx(j) > -1 .AND. xx(j) < 1) THEN
        uu(j) = (1.0 - xx(j)**2)**2
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
END SUBROUTINE   

FUNCTION flux(uu) RESULT(ff)
  !Flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: ff
  if (abs(uu**2) < precision) THEN
    ff = 0.0_dp
  ELSE
    ff = 0.5*uu**2
  END IF
END FUNCTION flux

FUNCTION DiffMat(uu) RESULT(kk)
  !Conservative diffusion in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: kk
  if (abs(uu**2) < precision) THEN
    kk = 0.0_dp
  ELSE
    kk = mu*uu**2
  END iF
END FUNCTION DiffMat

FUNCTION fluxEC(ul, ur) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: ul
  REAL(kind = dp), INTENT(IN)  :: ur
  REAL(kind = dp)              :: ff
  If (abs(ul) <= precision .AND. abs(ur) <= precision) THEN
    ff = 0.0_dp
  ELSE
    ff = 1/6.0*(ur**2+ur*ul+ul**2)
  END IF
END FUNCTION fluxEC

FUNCTION KKN(ul, ur) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: ul
  REAL(kind = dp), INTENT(IN)  :: ur
  REAL(kind = dp)              :: kk
  If (abs(ul) <= precision .AND. abs(ur) <= precision) THEN
    kk = 0.0_dp
  ELSE
    kk = mu*4/3*(ul**2+ul*ur+ur**2)/(ul+ur)
  END IF
END FUNCTION KKN

SUBROUTINE Cdt(uold, CFL, dx, dt)
  ! Update dt based on CFL condition
  REAL(kind = dp), INTENT(IN)  :: uold(:)
  REAL(kind = dp), INTENT(IN)  :: CFL
  REAL(kind = dp), INTENT(IN)  :: dx
  REAL(kind = dp)              :: dt
  dt = 0.0_dp
  dt = CFL/(1.0_dp/dx*MAXVAL(ABS(uold))+1.0_dp/dx**2*2*MAXVAL(ABS(2*mu*uold)))
END SUBROUTINE


END MODULE example_burger