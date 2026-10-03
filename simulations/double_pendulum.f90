! double_pendulum.f90
! Planar double pendulum: RK4 (Lagrangian form) and three symplectic
! integrators for the non-separable Hamiltonian H(theta,p):
!   - implicit midpoint rule            (order 2)
!   - generalized Stormer-Verlet        (order 2)
!   - Yoshida triple-jump composition   (order 4)
!
! Build:  gfortran -O2 -std=f2008 -Wall double_pendulum.f90 -o dp
!
! Conventions: theta_i measured from the downward vertical,
!   Lagrangian state  y = (th1, th2, w1, w2),   w = d(theta)/dt
!   canonical state   z = (th1, th2, p1, p2)
module dp_mod
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  private
  public :: wp, par_t, rhs, energy_v, ham, vel_to_mom, mom_to_vel
  public :: rk4_step, midpoint_step, sv_step, yoshida4_step

  type :: par_t
     real(wp) :: m1 = 1.0_wp, m2 = 1.0_wp
     real(wp) :: l1 = 1.0_wp, l2 = 1.0_wp
     real(wp) :: g  = 9.81_wp
  end type par_t

  integer,  parameter :: maxit = 200
  real(wp), parameter :: tol   = 1.0e-13_wp

contains

  !---------------------------------------------------------------
  ! Lagrangian right-hand side: solve A*(th1'',th2'')^T = b
  !---------------------------------------------------------------
  subroutine rhs(p, y, dy)
    type(par_t), intent(in)  :: p
    real(wp),    intent(in)  :: y(4)
    real(wp),    intent(out) :: dy(4)
    real(wp) :: d, c, s, mu, a11, a12, a21, a22, b1, b2, det

    d  = y(1) - y(2);  c = cos(d);  s = sin(d);  mu = p%m1 + p%m2
    a11 = mu*p%l1;  a12 = p%m2*p%l2*c
    a21 = p%l1*c;   a22 = p%l2
    b1  = -p%m2*p%l2*y(4)**2*s - mu*p%g*sin(y(1))
    b2  =  p%l1*y(3)**2*s      - p%g*sin(y(2))
    det = a11*a22 - a12*a21

    dy(1:2) = y(3:4)
    dy(3)   = (b1*a22 - a12*b2)/det
    dy(4)   = (a11*b2 - a21*b1)/det
  end subroutine rhs

  !---------------------------------------------------------------
  ! Total energy from Lagrangian variables (th, w)
  !---------------------------------------------------------------
  function energy_v(p, y) result(e)
    type(par_t), intent(in) :: p
    real(wp),    intent(in) :: y(4)
    real(wp) :: e, mu, t, v
    mu = p%m1 + p%m2
    t = 0.5_wp*mu*p%l1**2*y(3)**2 + 0.5_wp*p%m2*p%l2**2*y(4)**2        &
        + p%m2*p%l1*p%l2*y(3)*y(4)*cos(y(1) - y(2))
    v = -mu*p%g*p%l1*cos(y(1)) - p%m2*p%g*p%l2*cos(y(2))
    e = t + v
  end function energy_v

  !---------------------------------------------------------------
  ! Momenta <-> velocities (inverse mass matrix, eqs. thd1, thd2)
  !---------------------------------------------------------------
  subroutine vel_to_mom(p, th, w, mom)
    type(par_t), intent(in)  :: p
    real(wp),    intent(in)  :: th(2), w(2)
    real(wp),    intent(out) :: mom(2)
    real(wp) :: c, mu
    c = cos(th(1) - th(2));  mu = p%m1 + p%m2
    mom(1) = mu*p%l1**2*w(1) + p%m2*p%l1*p%l2*c*w(2)
    mom(2) = p%m2*p%l2**2*w(2) + p%m2*p%l1*p%l2*c*w(1)
  end subroutine vel_to_mom

  subroutine mom_to_vel(p, th, mom, w)
    type(par_t), intent(in)  :: p
    real(wp),    intent(in)  :: th(2), mom(2)
    real(wp),    intent(out) :: w(2)
    real(wp) :: c, s, mu, den
    c = cos(th(1) - th(2));  s = sin(th(1) - th(2));  mu = p%m1 + p%m2
    den  = p%m1 + p%m2*s*s
    w(1) = (p%l2*mom(1) - p%l1*c*mom(2)) / (p%l1**2*p%l2*den)
    w(2) = (mu*p%l1*mom(2) - p%m2*p%l2*c*mom(1)) / (p%m2*p%l1*p%l2**2*den)
  end subroutine mom_to_vel

  !---------------------------------------------------------------
  ! Hamiltonian, eq. (H), and its partial derivatives
  !---------------------------------------------------------------
  function ham(p, th, mom) result(h)
    type(par_t), intent(in) :: p
    real(wp),    intent(in) :: th(2), mom(2)
    real(wp) :: h, w(2), mu
    mu = p%m1 + p%m2
    call mom_to_vel(p, th, mom, w)
    h = 0.5_wp*(mom(1)*w(1) + mom(2)*w(2))                              &
        - mu*p%g*p%l1*cos(th(1)) - p%m2*p%g*p%l2*cos(th(2))
  end function ham

  subroutine h_p(p, th, mom, hp)             ! dH/dp = M^{-1} p = theta'
    type(par_t), intent(in)  :: p
    real(wp),    intent(in)  :: th(2), mom(2)
    real(wp),    intent(out) :: hp(2)
    call mom_to_vel(p, th, mom, hp)
  end subroutine h_p

  subroutine h_th(p, th, mom, hq)            ! dH/dtheta at fixed p
    type(par_t), intent(in)  :: p
    real(wp),    intent(in)  :: th(2), mom(2)
    real(wp),    intent(out) :: hq(2)
    real(wp) :: w(2), s, mu, k
    mu = p%m1 + p%m2
    call mom_to_vel(p, th, mom, w)
    s = sin(th(1) - th(2))
    k = p%m2*p%l1*p%l2*w(1)*w(2)*s
    hq(1) =  k + mu*p%g*p%l1*sin(th(1))
    hq(2) = -k + p%m2*p%g*p%l2*sin(th(2))
  end subroutine h_th

  subroutine vfield(p, z, f)                 ! f = (dH/dp, -dH/dtheta)
    type(par_t), intent(in)  :: p
    real(wp),    intent(in)  :: z(4)
    real(wp),    intent(out) :: f(4)
    real(wp) :: hq(2)
    call h_p (p, z(1:2), z(3:4), f(1:2))
    call h_th(p, z(1:2), z(3:4), hq)
    f(3:4) = -hq
  end subroutine vfield

  !---------------------------------------------------------------
  ! Classical RK4 on the Lagrangian first-order system
  !---------------------------------------------------------------
  subroutine rk4_step(p, y, h)
    type(par_t), intent(in)    :: p
    real(wp),    intent(inout) :: y(4)
    real(wp),    intent(in)    :: h
    real(wp) :: k1(4), k2(4), k3(4), k4(4)
    call rhs(p, y,              k1)
    call rhs(p, y + 0.5_wp*h*k1, k2)
    call rhs(p, y + 0.5_wp*h*k2, k3)
    call rhs(p, y +        h*k3, k4)
    y = y + h*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)/6.0_wp
  end subroutine rk4_step

  !---------------------------------------------------------------
  ! Implicit midpoint rule  z1 = z0 + h f((z0+z1)/2)
  ! (symplectic and symmetric for any H); fixed-point iteration,
  ! convergent for h*|Hessian H| < 1.  Use Newton for larger h.
  !---------------------------------------------------------------
  subroutine midpoint_step(p, z, h, ok)
    type(par_t), intent(in)    :: p
    real(wp),    intent(inout) :: z(4)
    real(wp),    intent(in)    :: h
    logical,     intent(out)   :: ok
    real(wp) :: f(4), zn(4), zt(4), err
    integer  :: it

    call vfield(p, z, f)
    zn = z + h*f                              ! explicit Euler predictor
    ok = .false.
    do it = 1, maxit
       call vfield(p, 0.5_wp*(z + zn), f)
       zt  = z + h*f
       err = maxval(abs(zt - zn))
       zn  = zt
       if (err < tol*(1.0_wp + maxval(abs(zn)))) then
          ok = .true.;  exit
       end if
    end do
    z = zn
  end subroutine midpoint_step

  !---------------------------------------------------------------
  ! Generalized (implicit) Stormer-Verlet, Hairer-Lubich-Wanner:
  !   p_{1/2} = p_0 - h/2 H_th(q_0, p_{1/2})
  !   q_1     = q_0 + h/2 [H_p(q_0,p_{1/2}) + H_p(q_1,p_{1/2})]
  !   p_1     = p_{1/2} - h/2 H_th(q_1, p_{1/2})
  !---------------------------------------------------------------
  subroutine sv_step(p, z, h, ok)
    type(par_t), intent(in)    :: p
    real(wp),    intent(inout) :: z(4)
    real(wp),    intent(in)    :: h
    logical,     intent(out)   :: ok
    real(wp) :: q0(2), p0(2), ph(2), pt(2), q1(2), qt(2)
    real(wp) :: hq(2), hp0(2), hp1(2), err
    integer  :: it
    logical  :: ok1, ok2

    q0 = z(1:2);  p0 = z(3:4)

    ! (1) implicit half step in p
    ph = p0;  ok1 = .false.
    do it = 1, maxit
       call h_th(p, q0, ph, hq)
       pt  = p0 - 0.5_wp*h*hq
       err = maxval(abs(pt - ph));  ph = pt
       if (err < tol*(1.0_wp + maxval(abs(ph)))) then
          ok1 = .true.;  exit
       end if
    end do

    ! (2) implicit full step in q
    call h_p(p, q0, ph, hp0)
    q1 = q0 + h*hp0;  ok2 = .false.
    do it = 1, maxit
       call h_p(p, q1, ph, hp1)
       qt  = q0 + 0.5_wp*h*(hp0 + hp1)
       err = maxval(abs(qt - q1));  q1 = qt
       if (err < tol*(1.0_wp + maxval(abs(q1)))) then
          ok2 = .true.;  exit
       end if
    end do

    ! (3) explicit half step in p
    call h_th(p, q1, ph, hq)
    z(1:2) = q1
    z(3:4) = ph - 0.5_wp*h*hq
    ok = ok1 .and. ok2
  end subroutine sv_step

  !---------------------------------------------------------------
  ! Order-4 composition (Yoshida 1990 / Suzuki triple jump):
  !   Psi4(h) = Psi2(g1 h) o Psi2(g0 h) o Psi2(g1 h)
  !---------------------------------------------------------------
  subroutine yoshida4_step(p, z, h, ok)
    type(par_t), intent(in)    :: p
    real(wp),    intent(inout) :: z(4)
    real(wp),    intent(in)    :: h
    logical,     intent(out)   :: ok
    real(wp), parameter :: g1 = 1.0_wp/(2.0_wp - 2.0_wp**(1.0_wp/3.0_wp))
    real(wp), parameter :: g0 = 1.0_wp - 2.0_wp*g1
    logical :: o1, o2, o3
    call sv_step(p, z, g1*h, o1)
    call sv_step(p, z, g0*h, o2)
    call sv_step(p, z, g1*h, o3)
    ok = o1 .and. o2 .and. o3
  end subroutine yoshida4_step

end module dp_mod

!=====================================================================
program dp_compare
  use dp_mod
  implicit none
  type(par_t) :: p
  real(wp) :: th0(2), w0(2), mom0(2), h, tfin, e0, emax(4)
  real(wp) :: yr(4), zm(4), zs(4), zy(4)
  integer  :: icase, n, nsteps, nprint
  logical  :: ok

  h      = 1.0e-2_wp
  tfin   = 1000.0_wp
  nsteps = nint(tfin/h)
  nprint = nsteps/10

  do icase = 1, 2
     if (icase == 1) then                      ! low energy: regular motion
        th0 = [0.4_wp, 0.2_wp]
        print '(a)', '# case 1: th0 = (0.4, 0.2), at rest (regular regime)'
     else                                      ! high energy: chaotic
        th0 = [2.0_wp, 2.5_wp]
        print '(a)', '# case 2: th0 = (2.0, 2.5), at rest (chaotic regime)'
     end if
     w0 = 0.0_wp
     call vel_to_mom(p, th0, w0, mom0)

     yr = [th0, w0]
     zm = [th0, mom0];  zs = zm;  zy = zm
     e0 = ham(p, th0, mom0)
     emax = 0.0_wp
     print '(a)', '#    t        RK4          midpoint     Stormer-Verlet  Yoshida4'

     do n = 1, nsteps
        call rk4_step(p, yr, h)
        call midpoint_step(p, zm, h, ok);   if (.not. ok) print *, 'midpoint: no conv.'
        call sv_step(p, zs, h, ok);         if (.not. ok) print *, 'SV: no conv.'
        call yoshida4_step(p, zy, h, ok);   if (.not. ok) print *, 'Y4: no conv.'

        emax(1) = max(emax(1), abs(energy_v(p, yr) - e0))
        emax(2) = max(emax(2), abs(ham(p, zm(1:2), zm(3:4)) - e0))
        emax(3) = max(emax(3), abs(ham(p, zs(1:2), zs(3:4)) - e0))
        emax(4) = max(emax(4), abs(ham(p, zy(1:2), zy(3:4)) - e0))

        if (mod(n, nprint) == 0) then
           print '(f8.1,4es14.3)', n*h,                                 &
              (energy_v(p, yr) - e0)/abs(e0),                           &
              (ham(p, zm(1:2), zm(3:4)) - e0)/abs(e0),                  &
              (ham(p, zs(1:2), zs(3:4)) - e0)/abs(e0),                  &
              (ham(p, zy(1:2), zy(3:4)) - e0)/abs(e0)
        end if
     end do
     print '(a,4es14.3)', '# max |dE/E0| :', emax/abs(e0)
     print *
  end do
end program dp_compare
