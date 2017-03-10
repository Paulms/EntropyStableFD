PROGRAM EntropyStableFD
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Entropy Stable Schemes for Degenerate Equations      !!
  !!                                                      !!
  !!                                                      !!
  !! Author: Paul Mendez Silva                            !!
  !! e-mail: paul.mendez@udec.cl                          !!
  !! Date: 10/Marzo/2016                                  !!
  !!                                                      !!
  !! Version: 0.1                                         !!  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE decimal               ! Precision
  USE plot                  ! Save data
  IMPLICIT NONE
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL
  REAL(kind = dp)           :: dt
  REAL(kind = dp)           :: mu
  INTEGER                   :: ntime
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:), uu(:), uold(:)
  REAL(kind = dp), ALLOCATABLE    :: fplus(:), fminus(:)
  INTEGER                   :: tt, i, j
  REAL(kind = dp), ALLOCATABLE    :: uleft, uright, fplusleft, fminusright
  REAL(kind = dp), ALLOCATABLE    :: Kleft, Kright
  CHARACTER(LEN=32)               :: name           ! File name to save plot data

  ! Initialize variables
  name = 'ouput.dat'
  Tend = 0.5_dp      ! Final Time
  N = 1600          ! Number of nodes
  CFL = 0.9
  mu = 0.1
  dx = 4.0_dp/(N-1)
  dt = CFL/(1.0_dp*(1/dx+4*mu/dx**2))
  ntime = floor(tend/dt)  !Number of time steps

    ! Allocate memory
  ALLOCATE(xx(N), uu(N))
  xx = 0.0_dp; uu = 0.0_dp
  xx = (/(i*dx-dx-2,i=1, N)/)      ! Location of grid points
  uu = xx                          ! Preallocation of solution

  ! Initial Conditions
  DO j = 1, N
    IF (xx(j) > -1 .AND. xx(j) < 1) THEN
      uu(j) = (1 - xx(j)**2)**2
    ELSE
      uu(j) = 0
    END IF
  END DO
  uleft = 0.0_dp; uright = 0.0_dp

  !Engquist-Osher Scheme with forward Euler
  ALLOCATE(fplus(N), fminus(N), uold(N), uinit(N))
  uinit = 0.0_dp;   fplus = 0.0_dp;   fminus = 0.0_dp
  uold = 0.0_dp
  uinit = uu
  DO tt = 1,ntime
    uold = uu
    DO j = 1,N
      IF (uold(j) > 0) THEN
        fplus(j) = Flux(uold(j)); fminus(j) = 0.0_dp
      ELSE
        fplus(j) = 0.0_dp; fminus(j) = Flux(uold(j))
      END IF
    END DO
    fplusleft = Flux(uleft); fminusright = Flux(uright)
    Kleft = DiffMat(uleft); Kright = DiffMat(uright)
    j = 1
    uu(j) = uold(j) - dt/dx * (fplus(j) + fminus(j+1) - fplusleft-fminus(j)) +&
    dt/dx**2*(DiffMat(uold(j+1)) - 2*DiffMat(uold(j)) + Kleft)
    DO j = 2,(N-1)
      uu(j) = uold(j) - dt/dx * (fplus(j) + fminus(j+1) - fplus(j-1)-fminus(j)) +&
      dt/dx**2*(DiffMat(uold(j+1)) - 2*DiffMat(uold(j)) + DiffMat(uold(j-1)))
    END DO
    j = N
    uu(j) = uold(j) - dt/dx * (fplus(j) + fminusright - fplus(j-1)-fminus(j)) +&
    dt/dx**2*(Kright - 2*DiffMat(uold(j)) + DiffMat(uold(j-1)))
  END DO
   CALL plot_results(uu, uinit, xx, name)
CONTAINS

FUNCTION flux(uu) RESULT(ff)
  !Flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: ff
  ff = 0.5*uu**2
END FUNCTION flux

FUNCTION DiffMat(uu) RESULT(kk)
  !Flux in Burguer's Equation
  REAL(kind = dp), INTENT(IN)  :: uu
  REAL(kind = dp)              :: kk
  kk = mu*uu**2
END FUNCTION DiffMat



END PROGRAM EntropyStableFD