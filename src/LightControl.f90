module LightControl
   use lightcontrol_pid, only: pid_controller
   implicit none
   private
   public :: pid_controller

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, LightControl!"
   end subroutine say_hello
end module LightControl
