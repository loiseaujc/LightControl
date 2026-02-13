program check
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use TestLyapunov, only: collect_test_lyapunov
   use LightControl
   implicit none

   ! Unit-test related.
   integer :: status, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   status = 0
   testsuites = [ &
                new_testsuite("Lyapunov Test Suite", collect_test_lyapunov) &
                ]

   do is = 1, size(testsuites)
      write (*, *) "-------------------------------"
      write (error_unit, fmt) "Testing :", testsuites(is)%name
      write (*, *) "-------------------------------"
      write (*, *)
      call run_testsuite(testsuites(is)%collect, error_unit, status)
      write (*, *)
   end do

   if (status > 0) then
      write (error_unit, '(i0, 1x, a)') status, "test(s) failed!"
      error stop
   else if (status == 0) then
      write (*, *) "All tests successfully passed!"
      write (*, *)
   end if

end program check
