! Test Problem 2
! u_t + (u^2)_x=(k(u)_x)_x   (x,t) in [-1,1]x[0,T]
MODULE test2
USE decimal
USE plot
USE tipos
USE numeric_schemes
IMPLICIT NONE
PUBLIC test2_run
PRIVATE
CONTAINS
subroutine test2_run()
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
  ! Zero variables
  Tend = 0.0_dp; error = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 0.15_dp      ! Final Time
  CFL = 0.9_dp
  M = 8000

  !Save reference solution (Change binary flag if needed)
  IF (.FALSE.) THEN
    name = 'test_2_reference'
    N = M          ! Number of nodes
    CALL setup_problem(dx, N, xx, uu, uinit)
    ntests = 1
    ALLOCATE(results(N, ntests+1), names(ntests+1))
    names = ['x    ', 'REF  ']
    results(:,1) = xx
    CALL Entropy_Conservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.2_dp*dx) !ESC-0.2
    results(:,2) = uu
    CALL save_matrix(results, names, name)
    DEALLOCATE(results, uu, names, uinit, xx)
  END IF

  !Compute errors (Change binary flag if needed)
  IF (.TRUE.) THEN
    name = 'test_2_reference'
    ALLOCATE(steps(5))
    steps = [200,400,800,1600]
    !Read reference solution
    CALL read_matrix(name , reference, M, 2)
    ! Prepare to save matrix
    name = 'test_2_errors'
    ntests = 3

    ALLOCATE(results(4, ntests+1), names(ntests+1))
    results = 0.0_dp
    names = ['N       ', 'MS      ', 'ESC-0.2 ', 'ESNC-0.2']
    results(:,1) = steps
    DO i = 1, 4
      N = steps(i)
      print *, "Starting numerical tests with N = ", N
      CALL setup_problem(dx, N, xx, uu, uinit)
      CALL Engquist_Osher(FORWARD_EULER, uu, N, Tend, dx, CFL, flux, DiffMat, Cdt)  !MS
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(i,2) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_Conservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.2_dp) !ESC
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(i,3) = error
      uu = uinit; error = 0.0_dp
      CALL Entropy_NonConservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.2_dp) !ESNC
      error = cumpute_errors(reference(:,2), M, uu, N)
      results(i,4) = error
      uu = uinit; error = 0.0_dp
      DEALLOCATE(uu, uinit, xx)
    END DO
    CALL save_matrix(results, names, name)
    DEALLOCATE(results, names)
    STOP
  END IF

  !Run numerical schemes
  N = 400          ! Number of nodes
  CALL setup_problem(dx, N, xx, uu, uinit)
  name = 'test_2_400'
  ntests = 3
  ALLOCATE(results(N, ntests+1), names(ntests+1))
  names = ['N       ', 'MS      ', 'ESC-0.2 ', 'ESNC-0.2']
  results(:,1) = xx
  CALL Engquist_Osher(FORWARD_EULER, uu, N, Tend, dx, CFL, flux, DiffMat, Cdt)  !MS
  results(:,2) = uu
  uu = uinit
  CALL Entropy_Conservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, DiffMat, Cdt, 0.2_dp) !ESC-0.2
  results(:,3) = uu
  uu = uinit
  CALL Entropy_NonConservative(FORWARD_EULER, .TRUE., uu, N, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.2_dp) !ESNC-0.2
  results(:,4) = uu
  CALL save_matrix(results, names, name)
  !CALL plot_results(uu, uinit, xx, name)
  ! Clean memory
  DEALLOCATE(results, uu, names, uinit, xx)
end subroutine test2_run

SUBROUTINE setup_problem(dx, N, xx, uu, uinit)
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
  DO j = 1, N
    IF (xx(j) <= -0.5) THEN
      uu(j) = 0.0_dp
    ELSEIF (-0.5 < xx(j) .AND. xx(j) < -0.3) THEN
      uu(j) = 5*(xx(j)+0.5)
    ELSEIF (-0.3 < xx(j) .AND. xx(j) < 0.3) THEN
      uu(j) = 1.0_dp
    ELSEIF (0.3 < xx(j) .AND. xx(j) < 0.5) THEN
      uu(j) = 5*(0.5 - xx(j))
    ELSE
      uu(j) = 0
    END IF
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
    ff = uu**2
  END IF
END FUNCTION flux

FUNCTION DiffFlux(uu) RESULT(kk)
  !Conservative diffusion in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu(:)
  REAL(kind = dp), ALLOCATABLE :: kk(:)
  INTEGER                      :: N, j
  N=SIZE(uu,1)
  ALLOCATE(kk(N))
  DO j = 1,N
    if (uu(j) < 0.5) THEN
      kk(j) = 0.0_dp
    ELSEIF (0.5 < uu(j) .AND. uu(j) < 0.6) THEN
      kk(j) = 2.5*uu(j)-1.25
    ELSE
      kk(j) = 0.25
    END iF
  END DO
END FUNCTION DiffFlux

FUNCTION DiffMat(uu) RESULT(kk)
  !Conservative diffusion in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: kk
  if (uu < 0.5) THEN
    kk = 0.0_dp
  ELSEIF (0.5 < uu .AND. uu < 0.6) THEN
    kk = 1.25*uu**2-1.25*uu+5.0/16.0
  ELSE
    kk = 0.25*uu-11/80
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
    ff = 1/3.0*(ur**2+ur*ul+ul**2)
  END IF
END FUNCTION fluxEC

FUNCTION KKN(ul, ur) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: ul
  REAL(kind = dp), INTENT(IN)  :: ur
  REAL(kind = dp)              :: kk
  If (abs(ul) <= precision .AND. abs(ur) <= precision) THEN
    kk = 0.0_dp
  ELSEIF (0.5 <= ul .AND. ul <= 0.6 .AND. 0.5 <= ur .AND. ur <= 0.6) THEN
    kk = 2.0/(ul+ur)*(5.0/6*(ul**2+ul*ur+ur**2)-1/8*(ul+ur))
  elseif (ul >= 0.6 .AND. ur >= 0.6) THEN
    kk = 2.0/(ul+ur)*(1.0/8*(ul+ur))
  else
    kk = 2/(ul+ur)*(r(ur)-r(ul))/(ur-ul)
  end if
END FUNCTION KKN

function r(u) RESULT (rr)
  REAL(kind = dp), INTENT(IN)  :: u
  REAL(kind = dp)              :: rr
  if (u <= 0.5) THEN
    rr = 0.0_dp
  elseif (0.5<u .AND. u <0.6) THEN
    rr = 5.0/6*u**3-5.0/8*u**2+5/96
  else
    rr = u**2/8-91/2400
  end if
end

SUBROUTINE Cdt(uold, CFL, dx, dt)
  ! Update dt based on CFL condition
  REAL(kind = dp), INTENT(IN)  :: uold(:)
  REAL(kind = dp), INTENT(IN)  :: CFL
  REAL(kind = dp), INTENT(IN)  :: dx
  REAL(kind = dp)              :: dt
  dt = 0.0_dp
  dt = CFL/(1.0_dp/dx*MAXVAL(ABS(2*uold))+1.0_dp/dx**2*2*MAXVAL(ABS(DiffFlux(uold))))
END SUBROUTINE
END MODULE test2