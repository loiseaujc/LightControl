module lightcontrol_pid
   use assert_m, only: assert => assert_always
   use stdlib_kinds, only: dp
   use stdlib_math, only: clip
   use stdlib_optval, only: optval
   implicit none
   private

   type, public :: pid_controller
      private
      !> Proprotional gain.
      real(dp) :: kp
      !> Integral gain.
      real(dp) :: ki
      !> Derivative gain.
      real(dp) :: kd
      !> Sampling time.
      real(dp) :: dt
      !> Internal state (for integral term).
      real(dp) :: integral_state
      !> Reference value.
      real(dp) :: r
      !> Upper and lower bounds for the actuation signal.
      real(dp) :: umin, umax
   contains
      private
      procedure, pass(self), public :: set_point
      procedure, pass(self), public :: output
   end type pid_controller

   interface pid_controller
      module function initialize_pid_controller(kp, ki, kd, Ts, umin, umax) result(pid)
         real(dp), optional, intent(in) :: kp
         real(dp), optional, intent(in) :: ki
         real(dp), optional, intent(in) :: kd
         real(dp), optional, intent(in) :: Ts
         real(dp), optional, intent(in) :: umin, umax
         type(pid_controller) :: pid
      end function initialize_pid_controller
   end interface

   interface
      pure module subroutine set_point(self, r)
         class(pid_controller), intent(inout) :: self
         real(dp), intent(in) :: r
      end subroutine set_point

      module function output(self, y) result(u)
         class(pid_controller), intent(inout) :: self
         real(dp), intent(in) :: y
         real(dp) :: u
      end function output
   end interface

contains

   module procedure initialize_pid_controller
   !> Proportional gain.
   pid%kp = optval(kp, 0.0_dp)
   call assert(assertion=abs(pid%kp) < huge(1.0_dp), &
               description="Proportional gain needs to have a finite value.")
   !> Integral gain.
   pid%ki = optval(ki, 0.0_dp)
   call assert(assertion=abs(pid%ki) < huge(1.0_dp), &
               description="Integral gain needs to have a finite value.")
   !> Derivative gain.
   pid%kd = optval(kd, 0.0_dp)
   call assert(assertion=pid%kd == 0.0_dp, &
               description="Derivative path is not yet supported.")
   call assert(assertion=abs(pid%kd) < huge(1.0_dp), &
               description="Derivative gain needs to have a finite value.")
   !> Sampling time.
   pid%dt = optval(Ts, 1.0_dp)
   call assert(assertion=pid%dt > 0.0_dp, &
               description="Sampling time needs to be strictly positive.")
   call assert(assertion=pid%dt < huge(1.0_dp), &
               description="Sampling time needs to have a finite value.")
   !> Internal state.
   pid%integral_state = 0.0_dp
   !> Set point.
   pid%r = 0.0_dp
   !> Upper and lower bounds for the actuation signal.
   pid%umin = optval(umin, -huge(1.0_dp))
   pid%umax = optval(umax, huge(1.0_dp))
   call assert(assertion=pid%umin <= pid%umax, &
               description="Lower bound for actuation signal needs to be larger than its upper bound.")
   end procedure initialize_pid_controller

   module procedure set_point
   call assert(assertion=abs(r) < huge(1.0_dp), &
               description="Set point needs to have a finite value.")
   self%r = r
   end procedure set_point

   module procedure output
   real(dp) :: err
   call assert(assertion=abs(y) < huge(1.0_dp), &
               description="Measurement needs to have a finite value.")
   !> Error signal.
   err = self%r - y
   !> Integral state.
   self%integral_state = self%integral_state + self%dt*err
   self%integral_state = clip(self%integral_state, self%umin, self%umax) ! Avoid integral wind-up.
   !> Output signal.
   u = self%kp*err + self%ki*self%integral_state
   !> Clamp the signal.
   u = clip(u, self%umin, self%umax)
   end procedure output
end module lightcontrol_pid
