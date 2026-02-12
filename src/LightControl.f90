module LightControl
   use assert_m, only: assert => assert_always
   use stdlib_kinds, only: dp
   use stdlib_optval, only: optval
   use lightcontrol_pid, only: pid_controller
   implicit none(type, external)
   private
   public :: pid_controller
   public :: lyap, dlyap, solve_lyapunov
   public :: ctrb_gramian, obs_gramian

   !--------------------------------------
   !-----     LYAPUNOV EQUATIONS     -----
   !--------------------------------------

   interface
      pure module function lyap(A, Q) result(X)
         real(dp), intent(in)  :: A(:, :)
         real(dp), intent(in)  :: Q(:, :)
         real(dp), allocatable :: X(:, :)
      end function lyap
   end interface

   interface
      pure module function dlyap(A, Q) result(X)
         real(dp), intent(in)  :: A(:, :)
         real(dp), intent(in)  :: Q(:, :)
         real(dp), allocatable :: X(:, :)
      end function dlyap
   end interface

   interface ctrb_gramian
      module function ctrb_gramian_siso(A, b, discrete) result(P)
         real(dp), intent(in)                     :: A(:, :)
         real(dp), intent(in), contiguous, target :: b(:)
         logical, intent(in), optional            :: discrete
         real(dp), allocatable                    :: P(:, :)
      end function ctrb_gramian_siso

      pure module function ctrb_gramian_mimo(A, B, discrete) result(P)
         real(dp), intent(in)          :: A(:, :)
         real(dp), intent(in)          :: B(:, :)
         logical, intent(in), optional :: discrete
         real(dp), allocatable         :: P(:, :)
      end function ctrb_gramian_mimo
   end interface ctrb_gramian

   interface obs_gramian
      module function obs_gramian_siso(A, c, discrete) result(Q)
         real(dp), intent(in)                     :: A(:, :)
         real(dp), intent(in), contiguous, target :: c(:)
         logical, intent(in), optional            :: discrete
         real(dp), allocatable                    :: Q(:, :)
      end function obs_gramian_siso

      pure module function obs_gramian_mimo(A, C, discrete) result(Q)
         real(dp), intent(in)          :: A(:, :)
         real(dp), intent(in)          :: C(:, :)
         logical, intent(in), optional :: discrete
         real(dp), allocatable         :: Q(:, :)
      end function obs_gramian_mimo
   end interface obs_gramian

   interface
      pure module subroutine solve_lyapunov(A, C, U, dico, op, factorized, job, scale, separation, ferr, wr, wi, iwork, dwork)
         real(dp), intent(inout), target         :: A(:, :)
         real(dp), intent(inout), target         :: C(:, :)
         real(dp), intent(inout), target         :: U(:, :)
         character(len=1), intent(in)            :: dico
         character(len=1), intent(in)            :: op
         logical, intent(in)                     :: factorized
         character(len=1), intent(in)            :: job
         real(dp), optional, intent(in)          :: scale
         real(dp), optional, intent(out)         :: separation
         real(dp), optional, intent(out)         :: ferr
         real(dp), optional, intent(out), target :: wr(:), wi(:)
         integer, optional, intent(out), target  :: iwork(:)
         real(dp), optional, intent(out), target :: dwork(:)
      end subroutine solve_lyapunov
   end interface

   interface
      pure module integer function lyapunov_workspace(n, dico, job, fact) result(ldwork)
         integer, intent(in) :: n
         !> Order of the system.
         character(len=1), intent(in) :: dico
         !> Discrete or continuous time Lyapunov equation.
         character(len=1), intent(in) :: job
         !> Computation to be performed:
         !  - job = "x": compute the solution only.
         !  - job = "s": compute the separation only.
         !  - job = "b": compute both the solution and the separation.
         character(len=1), intent(in) :: fact
         !> Whether A is already in real Schur form or not.
         !  - fact = "f": A has already been factorized.
         !  - fact = "n": A has not been factorized yet.
      end function lyapunov_workspace
   end interface
end module LightControl
