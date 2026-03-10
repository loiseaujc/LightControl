module TestRiccati
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_kinds, only: dp
   use stdlib_math, only: all_close
   use stdlib_linalg, only: eye
   use LightControl, only: care, dare
   implicit none(type, external)
   private

   public :: collect_test_riccati
contains
   subroutine collect_test_riccati(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Continuous-time Riccati", test_care), &
                  new_unittest("Discrete-time Riccati", test_dare) &
                  ]
   end subroutine collect_test_riccati

   subroutine test_care(error)
      type(error_type), allocatable, intent(out) :: error

      !---------------------------------
      !-----     SCIPY EXAMPLE     -----
      !---------------------------------
      block
         integer, parameter :: n = 2      ! Number of states.
         integer, parameter :: m = 1      ! Number of inputs.
         real(dp) :: A(n, n), b(n)        ! State-space model.
         real(dp) :: Q(n, n), r           ! LQR cost.
         real(dp) :: Xref(n, n)           ! Reference solution.
         real(dp), allocatable :: X(:, :) ! SLICOT solution.

         !> Problem's data.
         A(1, :) = [4.0_dp, 3.0_dp]
         A(2, :) = [-4.5_dp, -3.5_dp]

         b = [1.0_dp, -1.0_dp]

         Q(1, :) = [9.0_dp, 6.0_dp]
         Q(2, :) = [6.0_dp, 4.0_dp]

         r = 1.0_dp

         !> Riccati solver.
         X = care(A, b, Q, r)

         !> Reference solution.
         Xref(1, :) = [21.72792206_dp, 14.48528137_dp]
         Xref(2, :) = [14.48528137_dp, 9.65685425_dp]

         !> Check error.
         call check(error, all_close(X, Xref, abs_tol=1e-8_dp))
         if (allocated(error)) return

      end block
   end subroutine test_care

   subroutine test_dare(error)
      type(error_type), allocatable, intent(out) :: error

      !----------------------------------
      !-----     MATLAB EXAMPLE     -----
      !----------------------------------
      block
         integer, parameter :: n = 2      ! Number of states.
         integer, parameter :: m = 1      ! Number of inputs.
         real(dp) :: A(n, n), b(n)        ! State-space model.
         real(dp) :: Q(n, n), r           ! LQR cost.
         real(dp) :: Xref(n, n)           ! Reference solution.
         real(dp), allocatable :: X(:, :) ! SLICOT solution.

         !> Problem's data.
         A(1, :) = [-0.9_dp, -0.3_dp]
         A(2, :) = [0.7_dp, 0.1_dp]

         b = 1.0_dp

         Q(1, :) = [1.0_dp, 0.0_dp]
         Q(2, :) = [0.0_dp, 3.0_dp]

         r = 0.1_dp

         !> Riccati solver.
         X = dare(A, b, Q, r)

         !> Reference solution.
         Xref(1, :) = [4.7687_dp, 0.9438_dp]
         Xref(2, :) = [0.9438_dp, 3.2369_dp]

         !> Check error.
         call check(error, all_close(X, Xref, abs_tol=1e-4_dp))
         if (allocated(error)) return

      end block
   end subroutine test_dare

end module TestRiccati
