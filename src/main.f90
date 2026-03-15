program main
   use config
   use solver
   use io
   implicit none

   ! Initialise
   call allocate_arrays()
   call read_body('body.dat')
   call set_initial_conditions()
   call update_body()

   call system('mkdir -p outputs/data')
   call open_forces_file()
   call print_config()
   call log_forces()
   call log_fields()

   ! Time loop
   do step = 1, nt
      t = step * dt

      call update_body()
      call set_velocity_bc()
      call predict_velocity()
      call apply_forcing()
      call solve_pressure()
      call correct_velocity()

      if (mod(step, n_force) == 0) call log_forces()
      if (mod(step, n_field) == 0) call log_fields()
   end do

   call close_forces_file()

end program main
