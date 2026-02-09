module LightControl
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, LightControl!"
  end subroutine say_hello
end module LightControl
