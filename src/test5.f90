! Test Problem 5
! Parabolic System
MODULE test5
USE decimal
USE plot
USE tipos
USE numeric_schemes_nd
IMPLICIT NONE
PUBLIC test5_run
PRIVATE
REAL(kind = dp), PARAMETER :: mu = 0.001_dp
REAL(kind = dp), PARAMETER :: gr = 9.8_dp

CONTAINS
subroutine test5_run(run_ref)
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N, M
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:), uu(:,:)
  INTEGER                         :: i
  CHARACTER(LEN=32)               :: name           ! File name to save plot data
  REAL(kind = dp), ALLOCATABLE    :: results(:,:)
  CHARACTER(LEN=8), ALLOCATABLE  :: names(:)
  LOGICAL                         :: run_ref
  ! Zero variables
  Tend = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 0.2_dp      ! Final Time
  CFL = 0.9_dp
  M = 1000

  !Save reference solution (Change binary flag if needed)
  IF (run_ref) THEN
    name = 'test_5_reference'
    N = M          ! Number of nodes
    CALL setup_problem(dx, N, xx, uu, uinit)
    ALLOCATE(results(N, 3), names(3))
    names = ['X       ', 'REFh    ', 'REFq    ']
    results(:,1) = xx
    CALL Entropy_NonConservative_nd(FORWARD_EULER, .FALSE., uu, N, 2, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.0_dp, ZERO_FLUX) !ESC-0.2
    results(:,2) = uu(:,1)
    results(:,3) = uu(:,2)
    CALL save_matrix(results, names, name, 0)
    DEALLOCATE(results, uu, names, uinit, xx)
  END IF

  !Run numerical schemes
  N = 500          ! Number of nodes
  CALL setup_problem(dx, N, xx, uu, uinit)
  name = 'test_5_500'
  ALLOCATE(results(N, 5), names(5))
  names = ['X       ', 'h0      ','q0      ', 'ESCNh   ', 'ESCNq   ']
  results(:,1) = xx
  results(:,2) = uinit(:,1)
  results(:,3) = uinit(:,2)
  CALL Entropy_NonConservative_nd(FORWARD_EULER, .FALSE., uu, N, 2, Tend, dx, CFL, fluxEC, KKN, Cdt, 0.0_dp, ZERO_FLUX) !ESNC2
  results(:,4) = uu(:,1)
  results(:,5) = uu(:,2)
  CALL save_matrix(results, names, name, 0)
  !CALL plot_results(uu, uinit, xx, name)
  ! Clean memory
  DEALLOCATE(results, uu, names, uinit, xx)
end subroutine test5_run

SUBROUTINE setup_problem(dx, N, xx, uu, uinit)
  INTEGER, INTENT(IN)             :: N
  REAL(kind = dp)                 :: dx
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:), uu(:,:)
  INTEGER                         :: i, j
  dx = 10.0/N

    ! Allocate memory
  ALLOCATE(xx(N), uu(N,2))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx+dx/2-5.0_dp,i=0, N-1)/)      ! Location of grid points

  ! Initial Conditions
  DO i = 1,N
    if (xx(i) < 0.0) THEN
      uu(i,1) = 2.0_dp
    else
     uu(i,1) = 1.0_dp
   end if
  end do

  !Save initial condition
  ALLOCATE(uinit(N,2))
  uinit = 0.0_dp;   
  uinit = uu
END SUBROUTINE

FUNCTION vv(u)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: u(:)
  REAL(kind = dp)              :: vv(SIZE(u))
  vv = 0.0_dp
  vv = [gr*u(1)-0.5*(u(2)/u(1))**2, (u(2)/u(1))]
END FUNCTION vv

!k(u) numerico
FUNCTION KKN(ul, ur) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: ul(:)
  REAL(kind = dp), INTENT(IN)  :: ur(:)
  REAL(kind = dp)              :: kk(SIZE(ul,1),SIZE(ul,1))
  kk = 0.0_dp
  kk(1,1) = 0.5*mu*(sum(vv(ul)**2)+sum(vv(ur)**2))
  kk(2,2) = kk(1,1)
END FUNCTION KKN

!k(u)
FUNCTION DiffFlux(uu) RESULT(kk)
  !Conservative diffusion flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu(:)
  REAL(kind = dp), ALLOCATABLE :: kk(:,:)
  INTEGER                      :: N, j
  ALLOCATE(kk(2,2))
  kk = 0.0_dp
  kk(1,1) = mu*sum(vv(uu)**2)
  kk(2,2) = kk(1,1)
END FUNCTION DiffFlux

FUNCTION fluxEC(ul, ur) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: ul(:)
  REAL(kind = dp), INTENT(IN)  :: ur(:)
  REAL(kind = dp)              :: ff(SIZE(ul,1))
  ff = 0.0_dp
  ff(1) = 0.25*(ur(1)+ul(1))*(ur(2)/ur(1)+ul(2)/ul(1))
  ff(2) = (0.5*gr*0.25*(ur(1)+ul(1))**2+0.25*0.5*(ur(1)+ul(1))*(ur(2)/ur(1)+ul(2)/ul(1))**2)
END FUNCTION fluxEC

SUBROUTINE Cdt(uold, CFL, dx, dt)
  ! Update dt based on CFL condition
  REAL(kind = dp), INTENT(IN)  :: uold(:,:)
  REAL(kind = dp), INTENT(IN)  :: CFL
  REAL(kind = dp), INTENT(IN)  :: dx
  REAL(kind = dp)              :: dt
  REAL(kind = dp), ALLOCATABLE :: uu(:), fu(:), h(:), q(:)
  INTEGER                      :: i, N
  dt = 0.0_dp
  N = SIZE(uold,1)
  ALLOCATE(uu(N), fu(N),h(N),q(N))
  uu = 0.0_dp; fu = 0.0_dp; h = 0.0_dp; q = 0.0_dp
  DO i = 1,N
    uu(i) = SQRT(SUM(DiffFlux(uold(i,:))**2))
  END DO
  h = uold(:,1); q = uold(:,2)
  fu = sqrt(1+(gr*h-q**2/h**2)**2+4*q**2/h**2)
  dt = CFL/(1.0_dp/dx*MAXVAL(fu)+1.0_dp/dx**2*2*MAXVAL(uu))
  DEALLOCATE(uu, fu, h, q)
END SUBROUTINE
END MODULE test5