module TestRiccati
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_kinds, only: dp
   use stdlib_math, only: all_close
   use stdlib_linalg, only: eye
   use LightControl, only: solve_riccati, riccati_workspace
   implicit none

   public :: collect_test_riccati
contains
   subroutine collect_test_riccati(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("Continuous-time Riccati", test_care) &
                  ]
   end subroutine collect_test_riccati

   subroutine test_care(error)
      type(error_type), allocatable, intent(out) :: error

      !---------------------------------
      !-----     SCIPY EXAMPLE     -----
      !---------------------------------
      block
         integer, parameter :: m = 1, n = 2
         character(len=1), parameter :: dico = "c", scale = "n"
         real(dp) :: A(n, n), B(n, m), Q(n, n), R(m, m)
         integer :: i, j
         !> Problem's data.
         A(1, :) = [4.0_dp, 3.0_dp]
         A(2, :) = [-4.5_dp, -3.5_dp]

         B(:, 1) = [1.0_dp, -1.0_dp]

         Q(1, :) = [9.0_dp, 6.0_dp]
         Q(2, :) = [6.0_dp, 4.0_dp]

         R = 1.0_dp

         call solve_riccati(A, B, Q, R, dico, scale)

         do i = 1, n
            print *, (Q(i, j), j=1, n)
         end do
      end block
   end subroutine test_care

end module TestRiccati
