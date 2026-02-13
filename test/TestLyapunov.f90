module TestLyapunov
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_kinds, only: dp
   use stdlib_math, only: all_close
   use stdlib_linalg, only: eye
   use LightControl, only: lyap, dlyap
   implicit none

   public :: collect_test_lyapunov
contains
   subroutine collect_test_lyapunov(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Continuous-time Lyapunov", test_lyap), &
                  new_unittest("Discrete-time Lyapunov", test_dlyap) &
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

      !---------------------------------
      !-----     SCIPY EXAMPLE     -----
      !---------------------------------
      block
         real(dp) :: A(3, 3), Q(3, 3), X(3, 3), Xref(3, 3)
         !> Problem's data.
         A(1, :) = [-3.0_dp, -2.0_dp, 0.0_dp]
         A(2, :) = [-1.0_dp, -1.0_dp, 0.0_dp]
         A(3, :) = [0.0_dp, -5.0_dp, -1.0_dp]

         Q = eye(3)

         Xref(1, :) = [-0.75_dp, 0.875_dp, -3.75_dp]
         Xref(2, :) = [0.875_dp, -1.375_dp, 5.3125_dp]
         Xref(3, :) = [-3.75_dp, 5.3125_dp, -27.0625_dp]

         !> Solve Lyapunov equation.
         X = lyap(A, Q)

         !> Check error.
         call check(error, all_close(X, Xref, abs_tol=1e-4_dp))
         if (allocated(error)) return
      end block
   end subroutine test_lyap

   subroutine test_dlyap(error)
      type(error_type), allocatable, intent(out) :: error
      !---------------------------------
      !-----     SCIPY EXAMPLE     -----
      !---------------------------------
      block
         real(dp) :: A(2, 2), Q(2, 2), X(2, 2), Xref(2, 2)
         !> Problem's data
         A(1, 1) = 0.2_dp; A(1, 2) = 0.5_dp
         A(2, 1) = 0.7_dp; A(2, 2) = -0.9_dp

         Q = eye(2)

         Xref(1, :) = [0.70872893_dp, 1.43518822_dp]
         Xref(2, :) = [1.43518822_dp, -2.4266315_dp]

         !> Solve Lyapunov equation.
         X = dlyap(A, -Q)

         !> Check error.
         call check(error, all_close(X, Xref, abs_tol=1e-8_dp))
         if (allocated(error)) return
      end block

      !----------------------------------
      !-----     SLICOT EXAMPLE     -----
      !----------------------------------
      block
         real(dp) :: A(3, 3), Q(3, 3), X(3, 3), Xref(3, 3)
         !> Problem's data.
         A(1, :) = [3.0_dp, 1.0_dp, 1.0_dp]
         A(2, :) = [1.0_dp, 3.0_dp, 0.0_dp]
         A(3, :) = [0.0_dp, 0.0_dp, 3.0_dp]

         Q(1, :) = [25.0_dp, 24.0_dp, 15.0_dp]
         Q(2, :) = [24.0_dp, 32.0_dp, 8.0_dp]
         Q(3, :) = [15.0_dp, 8.0_dp, 40.0_dp]

         Xref(1, :) = [2.0_dp, 1.0_dp, 1.0_dp]
         Xref(2, :) = [1.0_dp, 3.0_dp, 0.0_dp]
         Xref(3, :) = [1.0_dp, 0.0_dp, 4.0_dp]

         !> Solve Lyapunov equation.
         X = dlyap(transpose(A), Q)

         !> Check error.
         call check(error, all_close(X, Xref, abs_tol=1e-8_dp))
         if (allocated(error)) return
      end block
   end subroutine test_dlyap
end module TestLyapunov
