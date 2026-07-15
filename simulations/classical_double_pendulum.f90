!===============================================================================
! classical_double_pendulum.f90
!
! Classical double pendulum simulation via Lagrangian mechanics.
! State: theta1, theta2 (angles from vertical, rad), omega1, omega2 (rad/s)
! Integrated with 4th-order Runge-Kutta (RK4), fixed time step.
!
! Equations of motion derived from the 2x2 mass (inertia) matrix M(theta):
!   M11 = (m1+m2) L1^2         M12 = m2 L1 L2 cos(theta1-theta2)
!   M22 = m2 L2^2
!   M(theta) * [theta1''; theta2''] = -C(theta,omega) - Grav(theta)
! where C are Coriolis/centrifugal terms and Grav are gravity torques.
! This is solved for [theta1''; theta2''] by inverting the 2x2 system
! directly (closed form) at every RK4 stage.
!
! Compile:  gfortran -O2 -o classical classical_double_pendulum.f90
! Run:      ./classical
!
! Output:
!   classical_trajectory.dat  - t, theta1, theta2, omega1, omega2,
!                                x1,y1,x2,y2 (bob positions), energy
!   Console summary table with an energy-conservation check.
!===============================================================================
program classical_double_pendulum
  implicit none
  integer,  parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846_dp

  ! ---------------- user-editable physical parameters ----------------
  real(dp), parameter :: m1 = 1.0_dp, m2 = 1.0_dp      ! bob masses (kg)
  real(dp), parameter :: L1 = 1.0_dp, L2 = 1.0_dp      ! rod lengths (m)
  real(dp), parameter :: grav = 9.81_dp                ! gravitational accel (m/s^2)

  real(dp), parameter :: theta1_0 = 120.0_dp * pi/180.0_dp  ! initial angle 1 (rad)
  real(dp), parameter :: theta2_0 = -10.0_dp * pi/180.0_dp  ! initial angle 2 (rad)
  real(dp), parameter :: omega1_0 = 0.0_dp
  real(dp), parameter :: omega2_0 = 0.0_dp

  real(dp), parameter :: dt    = 0.0005_dp   ! time step (s)
  real(dp), parameter :: t_end = 20.0_dp     ! total simulated time (s)
  integer,  parameter :: print_every = 20    ! write every Nth step to file
  ! ---------------------------------------------------------------------

  integer  :: nsteps, istep, iunit
  real(dp) :: t, y(4), k1(4), k2(4), k3(4), k4(4)
  real(dp) :: E0, Ef, Emax_dev, x1p, y1p, x2p, y2p, Ecur

  nsteps = nint(t_end/dt)

  open(newunit=iunit, file="classical_trajectory.dat", status="replace", action="write")
  write(iunit,'(A)') "# t  theta1  theta2  omega1  omega2  x1  y1  x2  y2  energy"

  y  = [theta1_0, omega1_0, theta2_0, omega2_0]
  E0 = total_energy(y)
  Emax_dev = 0.0_dp
  t = 0.0_dp

  do istep = 0, nsteps
     if (mod(istep, print_every) == 0) then
        call cartesian(y(1), y(3), x1p, y1p, x2p, y2p)
        Ecur = total_energy(y)
        write(iunit,'(10ES16.7)') t, y(1), y(3), y(2), y(4), x1p, y1p, x2p, y2p, Ecur
     end if

     Emax_dev = max(Emax_dev, abs(total_energy(y) - E0))

     k1 = deriv(y)
     k2 = deriv(y + 0.5_dp*dt*k1)
     k3 = deriv(y + 0.5_dp*dt*k2)
     k4 = deriv(y + dt*k3)
     y  = y + dt/6.0_dp * (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)

     t = t + dt
  end do

  Ef = total_energy(y)
  close(iunit)

  print '(A)', "===================================================================="
  print '(A)', "  Classical Double Pendulum -- Simulation Summary"
  print '(A)', "===================================================================="
  print '(A,F10.4,A,F10.4,A)', "  masses (m1, m2)        : ", m1, " kg, ", m2, " kg"
  print '(A,F10.4,A,F10.4,A)', "  rod lengths (L1, L2)   : ", L1, " m,  ", L2, " m"
  print '(A,F10.4,A)',         "  gravity g              : ", grav, " m/s^2"
  print '(A,F10.4,A,F10.4,A)', "  initial angles         : theta1=", theta1_0*180.0_dp/pi, &
                                 " deg, theta2=", theta2_0*180.0_dp/pi, " deg"
  print '(A,ES12.4,A,F10.3,A)',"  time step / duration   : dt=", dt, " s, t_end=", t_end, " s"
  print '(A,I10)',             "  total RK4 steps        : ", nsteps
  print '(A)', "--------------------------------------------------------------------"
  print '(A,F14.8,A)', "  initial energy  E0      : ", E0, " J"
  print '(A,F14.8,A)', "  final energy    Ef      : ", Ef, " J"
  print '(A,ES14.4,A)',"  max |E(t)-E0| over run  : ", Emax_dev, " J  (RK4 accuracy check)"
  print '(A)', "--------------------------------------------------------------------"
  print '(A)', "  Trajectory data written to classical_trajectory.dat"
  print '(A)', "===================================================================="

contains

  function deriv(yv) result(dydt)
    real(dp), intent(in) :: yv(4)
    real(dp) :: dydt(4)
    real(dp) :: th1, w1, th2, w2, dth
    real(dp) :: M11, M12, M22, detM
    real(dp) :: C1, C2, Grav1, Grav2, rhs1, rhs2

    th1 = yv(1);  w1 = yv(2)
    th2 = yv(3);  w2 = yv(4)
    dth = th1 - th2

    M11  = (m1+m2)*L1**2
    M12  = m2*L1*L2*cos(dth)
    M22  = m2*L2**2
    detM = M11*M22 - M12**2

    C1 = m2*L1*L2*sin(dth)*w2**2
    C2 = -m2*L1*L2*sin(dth)*w1**2
    Grav1 = (m1+m2)*grav*L1*sin(th1)
    Grav2 = m2*grav*L2*sin(th2)

    rhs1 = -(C1 + Grav1)
    rhs2 = -(C2 + Grav2)

    dydt(1) = w1
    dydt(3) = w2
    dydt(2) = ( M22*rhs1 - M12*rhs2) / detM
    dydt(4) = (-M12*rhs1 + M11*rhs2) / detM
  end function deriv

  function total_energy(yv) result(Etot)
    real(dp), intent(in) :: yv(4)
    real(dp) :: Etot, th1, w1, th2, w2, T, V
    th1 = yv(1); w1 = yv(2); th2 = yv(3); w2 = yv(4)
    T = 0.5_dp*(m1+m2)*L1**2*w1**2 + 0.5_dp*m2*L2**2*w2**2 &
        + m2*L1*L2*w1*w2*cos(th1-th2)
    V = -(m1+m2)*grav*L1*cos(th1) - m2*grav*L2*cos(th2)
    Etot = T + V
  end function total_energy

  subroutine cartesian(th1, th2, x1c, y1c, x2c, y2c)
    real(dp), intent(in)  :: th1, th2
    real(dp), intent(out) :: x1c, y1c, x2c, y2c
    x1c = L1*sin(th1)
    y1c = -L1*cos(th1)
    x2c = x1c + L2*sin(th2)
    y2c = y1c - L2*cos(th2)
  end subroutine cartesian

end program classical_double_pendulum
