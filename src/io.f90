module io
   use config
   implicit none

   integer, parameter :: forces_unit = 20
   integer, parameter :: fields_unit = 10

contains

   !=========================================================================
   ! CONFIGURATION
   !=========================================================================

   subroutine print_config()
      write(*,'(A,I0,A,I0,A,I0)') 'nx = ', nx, ', ny = ', ny, ', nS = ', nS
      write(*,'(A,ES0.2,A,I0)') 'dt = ', dt, ', nt = ', nt
      write(*,'(A,F0.1,A,ES0.2)') 'Re = ', Re, ', nu = ', nu
      write(*,'(A,F0.3,A,F0.3,A,F0.3)') 'Ax = ', Ax, ', Ay = ', Ay, ', Ath = ', Ath
      write(*,'(A,F0.3,A,F0.3,A,F0.3)') 'frx = ', frx, ', fry = ', fry, ', frth = ', frth
      write(*,'(A,A)') 'wall_type = ', wall_type
   end subroutine print_config

   !=========================================================================
   ! FORCES OUTPUT
   !=========================================================================

   subroutine open_forces_file()
      open(unit=forces_unit, file='outputs/data/forces.dat', status='replace')
      write(forces_unit,'(A)') '# t  CD  CL  CM  Xc  Yc  Th'
   end subroutine open_forces_file

   subroutine log_forces()
      write(*,'(A,I0,A,F0.5,A,I0,A,F0.4,A,F0.4,A,F0.4)') &
         'FORCE: step = ', step, ', t = ', t, &
         ', iter = ', iter, ', CD = ', CD, ', CL = ', CL, ', CM = ', CM
      write(forces_unit,'(7ES16.8)') t, CD, CL, CM, Xc, Yc, Th
   end subroutine log_forces

   subroutine close_forces_file()
      close(forces_unit)
   end subroutine close_forces_file

   !=========================================================================
   ! FIELDS OUTPUT
   !=========================================================================

   subroutine log_fields()
      integer :: i, j
      real(wp) :: div
      character(len=256) :: filename

      ! Compute divergence
      div = 0.0_wp
      !$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:div)
      do j = 1, ny
         do i = 1, nx
            div = max(div, abs((u(i,j) - u(i-1,j))/dx + (v(i,j) - v(i,j-1))/dy))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Print to screen
      write(*,'(A,I0,A,F0.5,A,I0,A,ES0.2,A,F0.4,A,F0.4,A,F0.4)') &
         'FIELD: step = ', step, ', t = ', t, ', iter = ', iter, &
         ', div = ', div, ', CD = ', CD, ', CL = ', CL, ', CM = ', CM

      ! Save to file
      write(filename,'(A,I6.6,A)') 'outputs/data/fields_', step, '.dat'
      open(unit=fields_unit, file=trim(filename), status='replace')
      write(fields_unit,'(A)') '# nx ny t Xc Yc Th'
      write(fields_unit,'(I0,1X,I0,1X,ES16.8,1X,ES16.8,1X,ES16.8,1X,ES16.8)') nx, ny, t, Xc, Yc, Th
      write(fields_unit,'(A)') '# x y u v p'
      do j = 1, ny
         do i = 1, nx
            write(fields_unit,'(5ES16.8)') xl + (real(i,wp) - 0.5_wp)*dx, yb + (real(j,wp) - 0.5_wp)*dy, &
               0.5_wp*(u(i-1,j) + u(i,j)), 0.5_wp*(v(i,j-1) + v(i,j)), p(i,j)
         end do
      end do
      close(fields_unit)
   end subroutine log_fields

end module io
