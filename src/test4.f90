! Test Problem 4
! Parabolic System
MODULE test4
USE decimal
USE plot
USE tipos
USE numeric_schemes_nd
IMPLICIT NONE
PUBLIC test4_run
PRIVATE
REAL(kind = dp), PARAMETER :: mu = 0.1_dp
CONTAINS
subroutine test4_run()
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:), uu(:,:)
  INTEGER                         :: i
  CHARACTER(LEN=32)               :: name           ! File name to save plot data
  REAL(kind = dp), ALLOCATABLE    :: results(:,:)
  CHARACTER(LEN=8), ALLOCATABLE  :: names(:)
  ! Zero variables
  Tend = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 3.0_dp      ! Final Time
  CFL = 0.9_dp

  !Run numerical schemes
  N = 1000          ! Number of nodes
  CALL setup_problem(dx, N, xx, uu, uinit)
  name = 'test_4_1000'
  ALLOCATE(results(N, 5), names(5))
  names = ['X       ', 'u1      ','u2      ', 'ESCNu1  ', 'ESCNu2  ']
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
end subroutine test4_run

SUBROUTINE setup_problem(dx, N, xx, uu, uinit)
  INTEGER, INTENT(IN)             :: N
  REAL(kind = dp)                 :: dx
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:), uu(:,:)
  INTEGER                         :: i, j
  dx = 5.0/N

    ! Allocate memory
  ALLOCATE(xx(N), uu(N,2))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx+dx/2-2.5_dp,i=0, N-1)/)      ! Location of grid points

  ! Initial Conditions
  DO i = 1,N
    if (-1.5 < xx(i) .and. xx(i) < -1.3) THEN
      uu(i,1) = 1.0_dp
      uu(i,2) = 0.0_dp
    else if (-0.5 < xx(i) .and. xx(i) < 0.5) THEN
      uu(i,1) = 1.0_dp
      uu(i,2) = 1.0_dp
    else if (1.3 < xx(i) .and. xx(i) < 1.5) THEN
      uu(i,1) = 0.0_dp
      uu(i,2) = 1.0_dp
    else
     uu(i,1) = 0.0_dp
     uu(i,2) = 0.0_dp
   end if
  end do

  !Save initial condition
  ALLOCATE(uinit(N,2))
  uinit = 0.0_dp;   
  uinit = uu
END SUBROUTINE   

!k(u) numerico
FUNCTION KKN(ul, ur) RESULT(kk)
  !Numerical viscosity matrix
  REAL(kind = dp), INTENT(IN)  :: ul(:)
  REAL(kind = dp), INTENT(IN)  :: ur(:)
  REAL(kind = dp)              :: kk(SIZE(ul,1),SIZE(ul,1))
  kk = 0.0_dp
  kk(1,1) = 0.5_dp*mu*(sum(ul**2)+sum(ur**2))
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
  kk(1,1) = mu*(uu(1)**2+uu(2)**2)
  kk(2,2) = kk(1,1)
END FUNCTION DiffFlux

FUNCTION fluxEC(ul, ur) RESULT(ff)
  !Entropy conservative flux
  REAL(kind = dp), INTENT(IN)  :: ul(:)
  REAL(kind = dp), INTENT(IN)  :: ur(:)
  REAL(kind = dp)              :: ff(SIZE(ul,1))
  ff = 0.0_dp
  ff = (ul**2+ul*ur+ur**2)/6.0
END FUNCTION fluxEC

SUBROUTINE Cdt(uold, CFL, dx, dt)
  ! Update dt based on CFL condition
  REAL(kind = dp), INTENT(IN)  :: uold(:,:)
  REAL(kind = dp), INTENT(IN)  :: CFL
  REAL(kind = dp), INTENT(IN)  :: dx
  REAL(kind = dp)              :: dt
  REAL(kind = dp), ALLOCATABLE :: uu(:)
  INTEGER                      :: i, N
  dt = 0.0_dp
  N = SIZE(uold,1)
  ALLOCATE(uu(N))
  uu = 0.0_dp
  DO i = 1,N
    uu(i) = SQRT(SUM(DiffFlux(uold(i,:))**2))
  END DO
  dt = CFL/(1.0_dp/dx*MAXVAL(sqrt(uold(:,1)**2+uold(:,2)**2))+1.0_dp/dx**2*2*MAXVAL(uu))
  DEALLOCATE(uu)
END SUBROUTINE
END MODULE test4