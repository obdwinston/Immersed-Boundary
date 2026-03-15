module solver
   use config
   implicit none

contains

   !=========================================================================
   ! INITIALISATION
   !=========================================================================

   subroutine allocate_arrays()
      allocate(u(0:nx, 0:ny+1), v(0:nx+1, 0:ny))
      allocate(u_star(0:nx, 0:ny+1), v_star(0:nx+1, 0:ny))
      allocate(p(0:nx+1, 0:ny+1), rhs(1:nx, 1:ny))
      allocate(fu(0:nx, 0:ny+1), fv(0:nx+1, 0:ny))
      allocate(X_local(nS), Y_local(nS))
      allocate(X(nS), Y(nS), dS(nS))
      allocate(Ux(nS), Uy(nS), Fx(nS), Fy(nS))
      allocate(Au_inv(nS,nS), Av_inv(nS,nS))
   end subroutine allocate_arrays

   subroutine set_initial_conditions()
      integer :: j

      u = 0.0_wp
      v = 0.0_wp
      p = 0.0_wp

      if (wall_type == 'slip') then
         u = U_inf
      else
         do j = 1, ny
            u(:,j) = 6.0_wp*U_inf * (real(j,wp)-0.5_wp)*dy * (Ly - (real(j,wp)-0.5_wp)*dy) / (Ly*Ly)
         end do
      end if

      step = 0
      t = 0.0_wp
      CD = 0.0_wp
      CL = 0.0_wp
      CM = 0.0_wp

      Xc = Xc0
      Yc = Yc0
      Th = Th0
      Vxc = Vx_mean + 2.0_wp*pi*frx*Ax*cos(phx)
      Vyc = Vy_mean + 2.0_wp*pi*fry*Ay*cos(phy)
      Om = Om_mean + 2.0_wp*pi*frth*Ath*cos(phth)
   end subroutine set_initial_conditions

   !=========================================================================
   ! UTILITIES
   !=========================================================================

   pure function delta(r, h) result(d)
      real(wp), intent(in) :: r, h
      real(wp) :: d, s

      s = abs(r) / h
      if (s >= 2.0_wp) then
         d = 0.0_wp
      else if (s >= 1.0_wp) then
         d = (1.0_wp/(8.0_wp*h)) * (5.0_wp - 2.0_wp*s - sqrt(-7.0_wp + 12.0_wp*s - 4.0_wp*s*s))
      else
         d = (1.0_wp/(8.0_wp*h)) * (3.0_wp - 2.0_wp*s + sqrt(1.0_wp + 4.0_wp*s - 4.0_wp*s*s))
      end if
   end function delta

   function inverse(A, n) result(A_inv)
      integer, intent(in) :: n
      real(wp), intent(in) :: A(n,n)
      real(wp) :: A_inv(n,n)
      integer :: i, j, k, imax
      real(wp) :: tmp, factor, aug(n,2*n)

      ! Build augmented matrix [A | I]
      aug = 0.0_wp
      do i = 1, n
         aug(i,1:n) = A(i,:)
         aug(i,n+i) = 1.0_wp
      end do

      ! Gauss-Jordan elimination
      do k = 1, n
         ! Partial pivoting
         imax = k
         do i = k+1, n
            if (abs(aug(i,k)) > abs(aug(imax,k))) imax = i
         end do
         if (imax /= k) then
            do j = 1, 2*n
               tmp = aug(k,j)
               aug(k,j) = aug(imax,j)
               aug(imax,j) = tmp
            end do
         end if

         ! Scale pivot row
         tmp = aug(k,k)
         aug(k,:) = aug(k,:) / tmp

         ! Eliminate column
         do i = 1, n
            if (i /= k) then
               factor = aug(i,k)
               aug(i,:) = aug(i,:) - factor*aug(k,:)
            end if
         end do
      end do

      ! Extract inverse
      do i = 1, n
         A_inv(i,:) = aug(i,n+1:2*n)
      end do
   end function inverse

   !=========================================================================
   ! BODY
   !=========================================================================

   subroutine read_body(filename)
      character(len=*), intent(in) :: filename
      integer :: k, k_next
      real(wp) :: Xc_local, Yc_local, chord, arc_total, seg_len, dSdx

      ! Read coordinates from file
      open(unit=30, file=filename, status='old', action='read')
      do k = 1, nS
         read(30,*) X_local(k), Y_local(k)
      end do
      close(30)

      ! Compute arc lengths
      arc_total = 0.0_wp
      do k = 1, nS
         k_next = mod(k, nS) + 1
         seg_len = sqrt((X_local(k_next) - X_local(k))**2 + (Y_local(k_next) - Y_local(k))**2)
         dS(k) = seg_len
         arc_total = arc_total + seg_len
      end do

      ! Centre at centroid
      Xc_local = 0.0_wp
      Yc_local = 0.0_wp
      do k = 1, nS
         Xc_local = Xc_local + X_local(k) * dS(k)
         Yc_local = Yc_local + Y_local(k) * dS(k)
      end do
      Xc_local = Xc_local / arc_total
      Yc_local = Yc_local / arc_total
      X_local = X_local - Xc_local
      Y_local = Y_local - Yc_local

      ! Normalise by chord length
      chord = 0.0_wp
      do k = 1, nS
         do k_next = k+1, nS
            seg_len = sqrt((X_local(k_next) - X_local(k))**2 + (Y_local(k_next) - Y_local(k))**2)
            if (seg_len > chord) chord = seg_len
         end do
      end do
      X_local = X_local / chord
      Y_local = Y_local / chord
      dS = dS / chord

      write(*,'(A,A)') 'BODY: ', trim(filename)
      dSdx = maxval(dS) * Lc / dx
      if (dSdx < 0.9_wp .or. dSdx > 1.1_wp) then
         write(*,'(A,F0.3,A)') 'WARNING: dS/dx = ', dSdx, ' (recommended: dS/dx ≈ 1.0)'
      end if
   end subroutine read_body

   subroutine update_body()
      integer :: k, l, i, j, il, ir, jb, jt
      real(wp) :: xu, yu, xv, yv, dk, dl
      real(wp), allocatable :: A_u(:,:), A_v(:,:)

      ! Update kinematics (position and velocity)
      Xc = Xc0 + Vx_mean*t + Ax*sin(2.0_wp*pi*frx*t + phx)
      Yc = Yc0 + Vy_mean*t + Ay*sin(2.0_wp*pi*fry*t + phy)
      Th = Th0 + Om_mean*t + Ath*sin(2.0_wp*pi*frth*t + phth)

      Vxc = Vx_mean + 2.0_wp*pi*frx*Ax*cos(2.0_wp*pi*frx*t + phx)
      Vyc = Vy_mean + 2.0_wp*pi*fry*Ay*cos(2.0_wp*pi*fry*t + phy)
      Om = Om_mean + 2.0_wp*pi*frth*Ath*cos(2.0_wp*pi*frth*t + phth)

      ! Transform local to global coordinates
      do k = 1, nS
         X(k) = Xc + Lc*(X_local(k)*cos(Th) - Y_local(k)*sin(Th))
         Y(k) = Yc + Lc*(X_local(k)*sin(Th) + Y_local(k)*cos(Th))
      end do

      ! Build spread-interpolate matrices
      allocate(A_u(nS,nS), A_v(nS,nS))

      ! Build A_u matrix for u-velocity
      A_u = 0.0_wp
      !$OMP PARALLEL DO PRIVATE(l,il,ir,jb,jt,i,j,xu,yu,dk,dl)
      do k = 1, nS
         do l = 1, nS
            il = max(0, max(floor((X(k)-xl)/dx), floor((X(l)-xl)/dx)) - 1)
            ir = min(nx, min(floor((X(k)-xl)/dx), floor((X(l)-xl)/dx)) + 3)
            jb = max(0, max(floor((Y(k)-yb)/dy + 0.5_wp), floor((Y(l)-yb)/dy + 0.5_wp)) - 1)
            jt = min(ny+1, min(floor((Y(k)-yb)/dy + 0.5_wp), floor((Y(l)-yb)/dy + 0.5_wp)) + 3)
            do j = jb, jt
               do i = il, ir
                  xu = xl + real(i,wp)*dx
                  yu = yb + (real(j,wp) - 0.5_wp)*dy
                  dk = delta(X(k) - xu, dx) * delta(Y(k) - yu, dy)
                  dl = delta(X(l) - xu, dx) * delta(Y(l) - yu, dy)
                  A_u(k,l) = A_u(k,l) + dk*dl * dx*dy * dS(l)*Lc
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! Build A_v matrix for v-velocity
      A_v = 0.0_wp
      !$OMP PARALLEL DO PRIVATE(l,il,ir,jb,jt,i,j,xv,yv,dk,dl)
      do k = 1, nS
         do l = 1, nS
            il = max(0, max(floor((X(k)-xl)/dx + 0.5_wp), floor((X(l)-xl)/dx + 0.5_wp)) - 1)
            ir = min(nx+1, min(floor((X(k)-xl)/dx + 0.5_wp), floor((X(l)-xl)/dx + 0.5_wp)) + 3)
            jb = max(0, max(floor((Y(k)-yb)/dy), floor((Y(l)-yb)/dy)) - 1)
            jt = min(ny, min(floor((Y(k)-yb)/dy), floor((Y(l)-yb)/dy)) + 3)
            do j = jb, jt
               do i = il, ir
                  xv = xl + (real(i,wp) - 0.5_wp)*dx
                  yv = yb + real(j,wp)*dy
                  dk = delta(X(k) - xv, dx) * delta(Y(k) - yv, dy)
                  dl = delta(X(l) - xv, dx) * delta(Y(l) - yv, dy)
                  A_v(k,l) = A_v(k,l) + dk*dl * dx*dy * dS(l)*Lc
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      Au_inv = inverse(A_u, nS)
      Av_inv = inverse(A_v, nS)

      deallocate(A_u, A_v)
   end subroutine update_body

   !=========================================================================
   ! BOUNDARY CONDITIONS
   !=========================================================================

   subroutine set_velocity_bc()
      integer :: j

      ! Inlet (left)
      if (wall_type == 'slip') then
         u(0,1:ny) = U_inf
      else
         do j = 1, ny
            u(0,j) = 6.0_wp*U_inf * (real(j,wp)-0.5_wp)*dy * (Ly - (real(j,wp)-0.5_wp)*dy) / (Ly*Ly)
         end do
      end if
      v(0,0:ny) = -v(1,0:ny)

      ! Outlet (right)
      v(nx+1,0:ny) = v(nx,0:ny)

      ! Top and bottom
      if (wall_type == 'slip') then
         u(0:nx,0) = u(0:nx,1)
         u(0:nx,ny+1) = u(0:nx,ny)
      else
         u(0:nx,0) = -u(0:nx,1)
         u(0:nx,ny+1) = -u(0:nx,ny)
      end if
      v(1:nx,0) = 0.0_wp
      v(1:nx,ny) = 0.0_wp
   end subroutine set_velocity_bc

   subroutine set_pressure_bc()
      ! Inlet (left): Neumann
      p(0,1:ny) = p(1,1:ny)
      ! Outlet (right): Dirichlet (p=0)
      p(nx+1,1:ny) = -p(nx,1:ny)
      ! Bottom: Neumann
      p(1:nx,0) = p(1:nx,1)
      ! Top: Neumann
      p(1:nx,ny+1) = p(1:nx,ny)

      ! Corners
      p(0,0) = p(1,1)
      p(0,ny+1) = p(1,ny)
      p(nx+1,0) = p(nx,1)
      p(nx+1,ny+1) = p(nx,ny)
   end subroutine set_pressure_bc

   !=========================================================================
   ! TIME STEPPING
   !=========================================================================

   subroutine predict_velocity()
      integer :: i, j
      real(wp) :: ue, uw, un, us, vn, vs
      real(wp) :: adv_x, adv_y, diff_x, diff_y

      u_star = u
      v_star = v

      !$OMP PARALLEL DO PRIVATE(i,ue,uw,un,us,adv_x,adv_y,diff_x,diff_y)
      do j = 1, ny
         do i = 1, nx-1
            ue = 0.5_wp*(u(i,j) + u(i+1,j))
            uw = 0.5_wp*(u(i-1,j) + u(i,j))
            un = 0.5_wp*(u(i,j) + u(i,j+1)) * 0.5_wp*(v(i,j) + v(i+1,j))
            us = 0.5_wp*(u(i,j-1) + u(i,j)) * 0.5_wp*(v(i,j-1) + v(i+1,j-1))
            adv_x = (ue*ue - uw*uw) / dx
            adv_y = (un - us) / dy
            diff_x = (u(i+1,j) - 2.0_wp*u(i,j) + u(i-1,j)) / (dx*dx)
            diff_y = (u(i,j+1) - 2.0_wp*u(i,j) + u(i,j-1)) / (dy*dy)
            u_star(i,j) = u(i,j) + dt*(-adv_x - adv_y + nu*(diff_x + diff_y))
         end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i,ue,uw,vn,vs,adv_x,adv_y,diff_x,diff_y)
      do j = 1, ny-1
         do i = 1, nx
            ue = 0.5_wp*(u(i,j) + u(i,j+1)) * 0.5_wp*(v(i,j) + v(i+1,j))
            uw = 0.5_wp*(u(i-1,j) + u(i-1,j+1)) * 0.5_wp*(v(i-1,j) + v(i,j))
            vn = 0.5_wp*(v(i,j) + v(i,j+1))
            vs = 0.5_wp*(v(i,j-1) + v(i,j))
            adv_x = (ue - uw) / dx
            adv_y = (vn*vn - vs*vs) / dy
            diff_x = (v(i+1,j) - 2.0_wp*v(i,j) + v(i-1,j)) / (dx*dx)
            diff_y = (v(i,j+1) - 2.0_wp*v(i,j) + v(i,j-1)) / (dy*dy)
            v_star(i,j) = v(i,j) + dt*(-adv_x - adv_y + nu*(diff_x + diff_y))
         end do
      end do
      !$OMP END PARALLEL DO
   end subroutine predict_velocity

   subroutine apply_forcing()
      integer :: k, l, i, j, il, ir, jb, jt
      real(wp) :: xu, yu, xv, yv, w
      real(wp) :: Fx_total, Fy_total, M_total
      real(wp), allocatable :: rhs_u(:), rhs_v(:)

      allocate(rhs_u(nS), rhs_v(nS))

      ! Interpolate velocity to Lagrangian points
      Ux = 0.0_wp
      Uy = 0.0_wp
      !$OMP PARALLEL DO PRIVATE(il,ir,jb,jt,i,j,xu,yu,xv,yv,w)
      do k = 1, nS
         il = max(0, floor((X(k)-xl)/dx) - 1)
         ir = min(nx, il + 3)
         jb = max(0, floor((Y(k)-yb)/dy + 0.5_wp) - 1)
         jt = min(ny+1, jb + 3)
         do j = jb, jt
            do i = il, ir
               xu = xl + real(i,wp)*dx
               yu = yb + (real(j,wp) - 0.5_wp)*dy
               w = delta(X(k) - xu, dx) * delta(Y(k) - yu, dy)
               Ux(k) = Ux(k) + u_star(i,j) * w*dx*dy
            end do
         end do

         il = max(0, floor((X(k)-xl)/dx + 0.5_wp) - 1)
         ir = min(nx+1, il + 3)
         jb = max(0, floor((Y(k)-yb)/dy) - 1)
         jt = min(ny, jb + 3)
         do j = jb, jt
            do i = il, ir
               xv = xl + (real(i,wp) - 0.5_wp)*dx
               yv = yb + real(j,wp)*dy
               w = delta(X(k) - xv, dx) * delta(Y(k) - yv, dy)
               Uy(k) = Uy(k) + v_star(i,j) * w*dx*dy
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! Compute forcing (velocity mismatch)
      do k = 1, nS
         rhs_u(k) = (Vxc - Om*(Y(k) - Yc) - Ux(k)) / dt
         rhs_v(k) = (Vyc + Om*(X(k) - Xc) - Uy(k)) / dt
      end do

      ! Solve for Lagrangian forces
      !$OMP PARALLEL DO PRIVATE(l)
      do k = 1, nS
         Fx(k) = 0.0_wp
         Fy(k) = 0.0_wp
         do l = 1, nS
            Fx(k) = Fx(k) + Au_inv(k,l)*rhs_u(l)
            Fy(k) = Fy(k) + Av_inv(k,l)*rhs_v(l)
         end do
      end do
      !$OMP END PARALLEL DO

      ! Spread forces to Eulerian grid
      fu = 0.0_wp
      fv = 0.0_wp
      !$OMP PARALLEL DO PRIVATE(il,ir,jb,jt,i,j,xu,yu,xv,yv,w)
      do k = 1, nS
         il = max(0, floor((X(k)-xl)/dx) - 1)
         ir = min(nx, il + 3)
         jb = max(0, floor((Y(k)-yb)/dy + 0.5_wp) - 1)
         jt = min(ny+1, jb + 3)
         do j = jb, jt
            do i = il, ir
               xu = xl + real(i,wp)*dx
               yu = yb + (real(j,wp) - 0.5_wp)*dy
               w = delta(X(k) - xu, dx) * delta(Y(k) - yu, dy)
               !$OMP ATOMIC
               fu(i,j) = fu(i,j) + Fx(k) * w * dS(k)*Lc
            end do
         end do

         il = max(0, floor((X(k)-xl)/dx + 0.5_wp) - 1)
         ir = min(nx+1, il + 3)
         jb = max(0, floor((Y(k)-yb)/dy) - 1)
         jt = min(ny, jb + 3)
         do j = jb, jt
            do i = il, ir
               xv = xl + (real(i,wp) - 0.5_wp)*dx
               yv = yb + real(j,wp)*dy
               w = delta(X(k) - xv, dx) * delta(Y(k) - yv, dy)
               !$OMP ATOMIC
               fv(i,j) = fv(i,j) + Fy(k) * w * dS(k)*Lc
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! Update predictor velocity
      u_star = u_star + dt*fu
      v_star = v_star + dt*fv
      u_star(nx,1:ny) = u_star(nx-1,1:ny)  ! outlet BC

      ! Compute force coefficients
      Fx_total = 0.0_wp
      Fy_total = 0.0_wp
      M_total = 0.0_wp
      !$OMP PARALLEL DO REDUCTION(+:Fx_total,Fy_total,M_total)
      do k = 1, nS
         Fx_total = Fx_total + Fx(k) * dS(k)*Lc
         Fy_total = Fy_total + Fy(k) * dS(k)*Lc
         M_total = M_total + ((X(k) - Xc)*Fy(k) - (Y(k) - Yc)*Fx(k)) * dS(k)*Lc
      end do
      !$OMP END PARALLEL DO

      CD = -2.0_wp*Fx_total / (U_inf**2 * Lc)
      CL = -2.0_wp*Fy_total / (U_inf**2 * Lc)
      CM = -2.0_wp*M_total / (U_inf**2 * Lc**2)

      deallocate(rhs_u, rhs_v)
   end subroutine apply_forcing

   subroutine solve_pressure()
      integer :: i, j, k
      real(wp) :: err, p_old

      ! Compute RHS (divergence of predictor velocity)
      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1, ny
         do i = 1, nx
            rhs(i,j) = (1.0_wp/dt) * ((u_star(i,j) - u_star(i-1,j))/dx + &
               (v_star(i,j) - v_star(i,j-1))/dy)
         end do
      end do
      !$OMP END PARALLEL DO

      ! Red-Black SOR iteration
      do k = 1, maxiter
         err = 0.0_wp

         ! Red sweep
         !$OMP PARALLEL DO PRIVATE(i,p_old) REDUCTION(max:err)
         do j = 1, ny
            do i = 1, nx
               if (mod(i+j,2) == 0) then
                  p_old = p(i,j)
                  p(i,j) = (1.0_wp - omega)*p(i,j) + omega/(2.0_wp/(dx*dx) + 2.0_wp/(dy*dy)) &
                     * ((p(i+1,j) + p(i-1,j))/(dx*dx) + (p(i,j+1) + p(i,j-1))/(dy*dy) - rhs(i,j))
                  err = max(err, abs(p(i,j) - p_old))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
         call set_pressure_bc()

         ! Black sweep
         !$OMP PARALLEL DO PRIVATE(i,p_old) REDUCTION(max:err)
         do j = 1, ny
            do i = 1, nx
               if (mod(i+j,2) == 1) then
                  p_old = p(i,j)
                  p(i,j) = (1.0_wp - omega)*p(i,j) + omega/(2.0_wp/(dx*dx) + 2.0_wp/(dy*dy)) &
                     * ((p(i+1,j) + p(i-1,j))/(dx*dx) + (p(i,j+1) + p(i,j-1))/(dy*dy) - rhs(i,j))
                  err = max(err, abs(p(i,j) - p_old))
               end if
            end do
         end do
         !$OMP END PARALLEL DO
         call set_pressure_bc()

         if (mod(k, 5000) == 0) then
            write(*,'(A,I0,A,F0.5,A,I0,A,ES0.2)') &
               'SOR: step = ', step, ', t = ', t, ', iter = ', k, ', err = ', err
         end if

         if (err < tol) exit
      end do

      iter = k
   end subroutine solve_pressure

   subroutine correct_velocity()
      integer :: i, j

      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1, ny
         do i = 1, nx
            u(i,j) = u_star(i,j) - dt*(p(i+1,j) - p(i,j))/dx
         end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i)
      do j = 1, ny-1
         do i = 1, nx
            v(i,j) = v_star(i,j) - dt*(p(i,j+1) - p(i,j))/dy
         end do
      end do
      !$OMP END PARALLEL DO
   end subroutine correct_velocity

end module solver
