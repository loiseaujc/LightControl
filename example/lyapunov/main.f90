program demo
   use stdlib_kinds, only: dp
   use stdlib_io_npy, only: save_npy
   use stdlib_linalg, only: outer_product
   use LightControl
   implicit none(type, external)
   integer, parameter :: n = 2
   integer :: i, j
   real(dp) :: A(n, n), Q(n, n), b(n), X(n, n)

   !> State-space model.
   A(1, 1) = 0.0_dp; A(1, 2) = 1.0_dp
   A(2, 1) = -1.0_dp; A(2, 2) = -1.0_dp

   b(1) = 0.0_dp
   b(2) = 1.0_dp

   !> Lyapunov equation.
   Q = -outer_product(b, b)
   X = lyap(A, Q)

   do i = 1, n
      print *, (A(i, j), j=1, n)
   end do

   print *

   do i = 1, n
      print *, (Q(i, j), j=1, n)
   end do

   print *

   do i = 1, n
      print *, (X(i, j), j=1, n)
   end do
end program demo
