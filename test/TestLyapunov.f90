module TestLyapunov
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_kinds, only: dp
   use stdlib_math, only: all_close
   use LightControl, only: lyap
   implicit none

   public :: collect_test_lyapunov
contains
   subroutine collect_test_lyapunov(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Continuous-time Lyapunov", test_lyap) &
                  ]
   end subroutine collect_test_lyapunov

   subroutine test_lyap(error)
      type(error_type), allocatable, intent(out) :: error
      !----------------------------------
      !-----     MATLAB EXAMPLE     -----
      !----------------------------------
      block
         real(dp) :: A(2, 2), Q(2, 2), X(2, 2), Xref(2, 2)
         !> Problem's data.
         A(1, 1) = 1.0_dp; A(1, 2) = 2.0_dp
         A(2, 1) = -3.0_dp; A(2, 2) = -4.0_dp

         Q(1, 1) = 3.0_dp; Q(1, 2) = 1.0_dp
         Q(2, 1) = 1.0_dp; Q(2, 2) = 1.0_dp

         Xref(1, 1) = 6.1667_dp; Xref(1, 2) = -3.8333_dp
         Xref(2, 1) = -3.8333_dp; Xref(2, 2) = 3.0_dp

         !> Solve Lyapunov equation.
         X = lyap(A, -Q)

         !> Check error.
         call check(error, all_close(X, Xref, abs_tol=1e-4_dp))
         if (allocated(error)) return
      end block
   end subroutine test_lyap
end module TestLyapunov
