module config
   implicit none

   integer, parameter :: wp = kind(1.0d0)
   real(wp), parameter :: pi = 3.14159265358979323846_wp

   ! Grid parameters
   integer, parameter :: nx = 700, ny = 400
   real(wp), parameter :: xl = -10.0_wp, xr = 25.0_wp
   real(wp), parameter :: yb = -10.0_wp, yt = 10.0_wp
   real(wp), parameter :: Lx = xr - xl
   real(wp), parameter :: Ly = yt - yb
   real(wp), parameter :: dx = Lx / real(nx, wp)
   real(wp), parameter :: dy = Ly / real(ny, wp)

   ! Body parameters
   integer, parameter :: nS = 64
   real(wp), parameter :: Lc = 1.0_wp
   real(wp), parameter :: Xc0 = 0.0_wp, Yc0 = 0.0_wp, Th0 = 0.0_wp
   real(wp), parameter :: Vx_mean = 0.0_wp, Vy_mean = 0.0_wp, Om_mean = 0.0_wp
   real(wp), parameter :: Ax = 0.0_wp, Ay = 0.2_wp, Ath = 0.0_wp
   real(wp), parameter :: frx = 0.0_wp, fry = 0.195_wp, frth = 0.0_wp
   real(wp), parameter :: phx = 0.0_wp, phy = 0.0_wp, phth = 0.0_wp

   ! Flow parameters
   real(wp), parameter :: U_inf = 1.0_wp
   real(wp), parameter :: Re = 185.0_wp
   real(wp), parameter :: nu = U_inf * Lc / Re
   character(len=6), parameter :: wall_type = 'slip'

   ! Time parameters
   real(wp), parameter :: dt = 0.002_wp
   integer, parameter :: nt = 100000
   integer, parameter :: n_field = 100
   integer, parameter :: n_force = 50

   ! Solver parameters
   real(wp), parameter :: tol = 1.0e-5_wp
   integer, parameter :: maxiter = 200000
   real(wp), parameter :: omega = 1.9_wp

   ! Eulerian arrays
   real(wp), allocatable :: u(:,:), v(:,:)
   real(wp), allocatable :: u_star(:,:), v_star(:,:)
   real(wp), allocatable :: p(:,:), rhs(:,:)
   real(wp), allocatable :: fu(:,:), fv(:,:)

   ! Lagrangian arrays
   real(wp), allocatable :: X_local(:), Y_local(:)
   real(wp), allocatable :: X(:), Y(:), dS(:)
   real(wp), allocatable :: Ux(:), Uy(:), Fx(:), Fy(:)
   real(wp), allocatable :: Au_inv(:,:), Av_inv(:,:)

   ! Runtime state
   integer :: step, iter
   real(wp) :: t
   real(wp) :: Xc, Yc, Th
   real(wp) :: Vxc, Vyc, Om
   real(wp) :: CD, CL, CM

end module config
