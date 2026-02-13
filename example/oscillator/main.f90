program demo
   use stdlib_kinds, only: dp
   use stdlib_io_npy, only: save_npy
   use LightControl
   implicit none

   integer, parameter :: n = 2, p = 1, q = 1
   real(dp), parameter :: dt = 1.0e-2_dp
   integer, parameter :: nsteps = int(20.0_dp/dt)
   real(dp), parameter :: r = 1.0_dp
   real(dp) :: A(n, n), B(n), C(n)
   type(pid_controller) :: pid

   !> Get state-space system.
   call state_space_system(A, B, C)

   !> Uncontrolled response.
   block
      real(dp) :: x(n), y(0:nsteps), u, e
      integer :: i
      !> Initial condition.
      x = 0.0_dp; y(0) = dot_product(C, x)
      !> Integrate system.
      do i = 1, nsteps
         e = r - y(i - 1)
         u = e
         x = x + dt*matmul(A, x) + dt*B*u
         y(i) = dot_product(C, x)
      end do
      call save_npy("example/oscillator/uncontrolled_response.npy", y)
   end block

   !> Proportional control.
   block
      real(dp) :: x(n), y(0:nsteps), u, e
      integer :: i
      real(dp), parameter :: kp = 5.0_dp, ki = 0.0_dp
      !> Initialize controller.
      pid = pid_controller(kp=kp, ki=ki, Ts=dt)
      call pid%set_point(r)

      !> Initial condition.
      x = 0.0_dp; y(0) = dot_product(C, x)
      do i = 1, nsteps
         u = pid%output(y(i - 1))
         x = x + dt*matmul(A, x) + dt*B*u
         y(i) = dot_product(C, x)
      end do
      call save_npy("example/oscillator/P_response.npy", y)
   end block

   !> Proportional-Integral control.
   block
      real(dp) :: x(n), y(0:nsteps), u, e
      integer :: i
      real(dp), parameter :: kp = 5.0_dp, ki = 2.0_dp
      !> Initialize controller.
      pid = pid_controller(kp=kp, ki=ki, Ts=dt)
      call pid%set_point(r)

      !> Initial condition.
      x = 0.0_dp; y(0) = dot_product(C, x)
      do i = 1, nsteps
         u = pid%output(y(i - 1))
         x = x + dt*matmul(A, x) + dt*B*u
         y(i) = dot_product(C, x)
      end do
      call save_npy("example/oscillator/PI_response.npy", y)
   end block

   !> Proportional-Integral control with actuator saturation.
   block
      real(dp) :: x(n), y(0:nsteps), u, e
      integer :: i
      real(dp), parameter :: kp = 5.0_dp, ki = 2.0_dp
      real(dp), parameter :: umin = 0.0_dp, umax = 1.0_dp
      !> Initialize controller.
      pid = pid_controller(kp=kp, ki=ki, Ts=dt, umin=umin, umax=umax)
      call pid%set_point(r)

      !> Initial condition.
      x = 0.0_dp; y(0) = dot_product(C, x)
      do i = 1, nsteps
         u = pid%output(y(i - 1))
         x = x + dt*matmul(A, x) + dt*B*u
         y(i) = dot_product(C, x)
      end do
      call save_npy("example/oscillator/PI_saturation_response.npy", y)
   end block

contains
   pure subroutine state_space_system(A, B, C)
      real(dp), intent(out) :: A(:, :), B(:), C(:)
      !> Dynamics matrix.
      A(1, 1) = 0.0_dp; A(1, 2) = 1.0_dp
      A(2, 1) = -1.0_dp; A(2, 2) = -1.0_dp
      !> Input-to-state matrix.
      B(1) = 0.0_dp; B(2) = 1.0_dp
      !> Measurement matrix.
      C(1) = 1.0_dp; C(2) = 0.0_dp
   end subroutine state_space_system
end program demo
