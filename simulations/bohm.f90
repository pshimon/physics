program node_count_solver
implicit none
integer, parameter :: dp = kind(1.0d0)
integer :: n_target, nodes
real(dp) :: Emin, Emax, E, tol
real(dp) :: xL, xR
integer :: max_iter, iter

! Physical constants (can set to 1 for simplicity)
real(dp), parameter :: hbar = 1.0_dp
real(dp), parameter :: mass = 1.0_dp

! Domain
xL = -10.0_dp
xR =  10.0_dp

! Target quantum number
n_target = 5

! Energy search interval
Emin = 0.0_dp
Emax = 20.0_dp

tol = 1.0e-10_dp
max_iter = 100

call find_eigenvalue(n_target, Emin, Emax, tol, E)

print *, "----------------------------------------"
print *, "Eigenvalue for n =", n_target, " is:"
print *, E
print *, "----------------------------------------"

contains

!==========================================
function V(x) result(val)
real(dp), intent(in) :: x
real(dp) :: val
! Harmonic oscillator V = (1/2) x^2
val = 0.5_dp * x*x
end function
!==========================================

subroutine schrodinger_rhs(x, y, dydx, E)
real(dp), intent(in) :: x, y(2), E
real(dp), intent(out) :: dydx(2)

! y(1) = phi
! y(2) = phi'
dydx(1) = y(2)
dydx(2) = (2.0_dp*mass/hbar**2)*(V(x)-E)*y(1)

end subroutine
!==========================================

subroutine rk4_step(x, y, h, E)
real(dp), intent(inout) :: x, y(2)
real(dp), intent(in) :: h, E
real(dp) :: k1(2), k2(2), k3(2), k4(2), yt(2)

call schrodinger_rhs(x, y, k1, E)

yt = y + 0.5_dp*h*k1
call schrodinger_rhs(x+0.5_dp*h, yt, k2, E)

yt = y + 0.5_dp*h*k2
call schrodinger_rhs(x+0.5_dp*h, yt, k3, E)

yt = y + h*k3
call schrodinger_rhs(x+h, yt, k4, E)

y = y + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
x = x + h

end subroutine
!==========================================

function count_nodes(E) result(nodes)
real(dp), intent(in) :: E
integer :: nodes
real(dp) :: x, h
real(dp) :: y(2), y_prev
integer :: i, Nsteps

nodes = 0
Nsteps = 20000
h = (xR - xL)/real(Nsteps,dp)

! Initial conditions for Dirichlet boundary
x = xL
y(1) = 0.0_dp
y(2) = 1.0_dp   ! arbitrary slope

y_prev = y(1)

do i = 1, Nsteps
    call rk4_step(x, y, h, E)
    if (y(1)*y_prev < 0.0_dp) nodes = nodes + 1
    y_prev = y(1)
end do

end function
!==========================================

subroutine find_eigenvalue(n_target, Emin, Emax, tol, Eout)
integer, intent(in) :: n_target
real(dp), intent(in) :: Emin, Emax, tol
real(dp), intent(out) :: Eout
real(dp) :: E_low, E_high, E_mid
integer :: n_low, n_mid
integer :: iter, max_iter

max_iter = 100

E_low = Emin
E_high = Emax

n_low = count_nodes(E_low)

do iter = 1, max_iter

    E_mid = 0.5_dp*(E_low + E_high)
    n_mid = count_nodes(E_mid)

    if (n_mid > n_target) then
        E_high = E_mid
    else
        E_low = E_mid
    end if

    if (abs(E_high - E_low) < tol) exit

end do

Eout = 0.5_dp*(E_low + E_high)

end subroutine
!==========================================

end program
