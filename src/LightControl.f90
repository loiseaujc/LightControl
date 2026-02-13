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
      !!    ### Description
      !!
      !!    Solve the (continuous-time) Lyapunov equation
      !!
      !!    \[
      !!        AX + XA^{\top} = Q
      !!    \]
      !!
      !!    where both \(A\) and \(Q\) are \(n \times n\) matrices.
      !!    The matrix \(Q\) morever needs to be symmetric (\(C = C^\top\)).
      !!    The solution matrix \( X \) is an \(n \times n\) symmetric positive definite matrix.
      !!
      !!    ### Syntax
      !!
      !!    ```fortran
      !!        X = lyap(A, Q)
      !!    ```
      !!
      !!    ### Arguments
      !!
      !!    - `A`   :   Double precision rank 2 array of size `n x n`.
      !!                It is an `intent(in)` argument.
      !!
      !!    - `Q`   :   Double precision rank 2 array of size `n x n`.
      !!                It is an `intent(in)` argument.
      !!
      !!    - `X`   :   Double precision rank 2 array of size `n x n`.
      !!                Returned by the function.
      pure module function lyap(A, Q) result(X)
         real(dp), intent(in)  :: A(:, :)
         !> Dynamics matrix of dimension `n x n`.
         real(dp), intent(in)  :: Q(:, :)
         !> Symmetric matrix of dimension `n x n` corresponding to the right-hand
         !> side of the Lyapunov equaiton.
         real(dp), allocatable :: X(:, :)
         !> Solution of the Lyapunov equation, dimension `n x n`.
      end function lyap
   end interface

   interface
      !!    ### Description
      !!
      !!    Solve the (discrete-time) Lyapunov equation
      !!
      !!    \[
      !!        AXA^{\top} - X = Q
      !!    \]
      !!
      !!    where both \(A\) and \(Q\) are \(n \times n\) matrices.
      !!    The matrix \(Q\) morever needs to be symmetric (\(C = C^\top\)).
      !!    The solution matrix \( X \) is an \(n \times n\) symmetric positive definite matrix.
      !!
      !!    ### Syntax
      !!
      !!    ```fortran
      !!        X = dlyap(A, Q)
      !!    ```
      !!
      !!    ### Arguments
      !!
      !!    - `A`   :   Double precision rank 2 array of size `n x n`.
      !!                It is an `intent(in)` argument.
      !!
      !!    - `Q`   :   Double precision rank 2 array of size `n x n`.
      !!                It is an `intent(in)` argument.
      !!
      !!    - `X`   :   Double precision rank 2 array of size `n x n`.
      !!                Returned by the function.
      pure module function dlyap(A, Q) result(X)
         real(dp), intent(in)  :: A(:, :)
         real(dp), intent(in)  :: Q(:, :)
         real(dp), allocatable :: X(:, :)
      end function dlyap
   end interface

   interface ctrb_gramian
      !!    ### Description
      !!
      !!    Given an LTI system in state-space form
      !!
      !!    \[
      !!        \begin{aligned}
      !!            \sigma x & = Ax + Bu \\
      !!            y & = Cx + Du
      !!        \end{aligned}
      !!    \]
      !!
      !!    where it is assumed that \(A\) is stable, and \(\sigma\) denotes either
      !!    the time-derivative (continuous-time) or shift operator (discrete-time),
      !!    this function computes the corresponding controllability gramian. This
      !!    gramian is solution to
      !!
      !!    \[
      !!        AP + PA^{\top} = - BB^\top
      !!    \]
      !!
      !!    for a continuous-time system, and
      !!
      !!    \[
      !!        APA^{\top} - P = -BB^\top
      !!    \]
      !!
      !!    for a discrete-time one.
      !!
      !!    ### Syntax
      !!
      !!    ```fortran
      !!        P = ctrb_gramian(A, B)
      !!    ```
      !!
      !!    ### Arguments
      !!
      !!    - `A`   :   Double precision rank 2 array of dimension `n x n`.
      !!                It is an `intent(in)` argument.
      !!
      !!    - `B`   :   Double precision rank 1 array of dimension `n`, or rank 2 array of
      !!                dimension `n x m`. It is an `intent(in)` argument.
      !!
      !!    - `P`   :   Double precision rank 2 array.
      !!                Returned by the function.
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
      !!    ### Description
      !!
      !!    Given an LTI system in state-space form
      !!
      !!    \[
      !!        \begin{aligned}
      !!            \sigma x & = Ax + Bu \\
      !!            y & = Cx + Du
      !!        \end{aligned}
      !!    \]
      !!
      !!    where it is assumed that \(A\) is stable, and \(\sigma\) denotes either
      !!    the time-derivative (continuous-time) or shift operator (discrete-time),
      !!    this function computes the corresponding observability gramian. This
      !!    gramian is solution to
      !!
      !!    \[
      !!        A^{\top}Q + QA = - C^{\top}C
      !!    \]
      !!
      !!    for a continuous-time system, and
      !!
      !!    \[
      !!        A^{\top}QA - Q = -C^{\top}C
      !!    \]
      !!
      !!    for a discrete-time one.
      !!
      !!    ### Syntax
      !!
      !!    ```fortran
      !!        Q = obs_gramian(A, C)
      !!    ```
      !!
      !!    ### Arguments
      !!
      !!    - `A`   :   Double precision rank 2 array of dimension `n x n`.
      !!                It is an `intent(in)` argument.
      !!
      !!    - `C`   :   Double precision rank 1 array of dimension `n`, or rank 2 array
      !!                of dimension `m x n`. It is an `intent(in)` argument.
      !!
      !!    - `P`   :   Double precision rank 2 array.
      !!                Returned by the function.
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
      !!    ### Description
      !!
      !!    Low-level wrapper for the `sb03md` subroutine from `slicot`. It is intended to
      !!    solve the real continuous-time Lyapunov equation
      !!
      !!    \[
      !!        \mathrm{op}(A)^{\top} X + X \mathrm{op}(A) = \mathrm{scale} \cdot C
      !!    \]
      !!
      !!    or the real discrete-time Lyapunov equation
      !!
      !!    \[
      !!        \mathrm{op}(A)^{\top} X \mathrm{op}(A) - X = \mathrm{scale} \cdot C
      !!    \]
      !!
      !!    and/or estimate an associated condition number (called separation), where
      !!    \(\mathrm{op}(A) = A\) or \(A^{\top}\) and \(C\) is symmetric. \(A\), \(C\)
      !!    and \(X\) are \(n \times n\) matrices, while \( 0 < \mathrm{scale} \leq 1\)
      !!    is an output set by the routine to avoid overflow in \( X \).
      !!
      !!    ### Syntax
      !!
      !!    ```fortran
      !!        call solve_lyapunov(A, C, U, dico, op, factorized, job [, scale, separation, ferr, wr, wi, iwork, dwork])
      !!    ```
      !!
      !!    ### Arguments
      !!
      !!    - `A`   :   Double precision rank-2 array of size `n x n`. On entry, it contains
      !!                the matrix `A`. If `factorized = .true.`, then `A` contains an upper
      !!                quasi-triangular matrix in Schur canonical form: the elements below
      !!                the upper Hessenberg part of `A` are not referenced and `A` is unchanged
      !!                on exit. If `factorized = .false.`, `A` is overwritten on exit by its
      !!                Schur canonical form. It is an `intent(inout)` argument.
      !!
      !!    - `C`   :   Double precision rank-2 array of size `n x n`. On entry, if `job = "x"`
      !!                or `"b"`, it contains the symmetric matrix `C` and is overwritten
      !!                by the symmetric solution matrix `X` on exit. If `job = "s"`, `C`
      !!                is not referenced.
      !!
      !!    - `U`   :   Double precision rank-2 array of size `n x n`. If `factorized = .true.`,
      !!                it must contains the orthogonal matrix of the real Schur factorization
      !!                of `A` and is unchanged on exit. If `factorized = .false.`, its input
      !!                value is irrelevant and will be overwritten by the orthogonal matrix
      !!                of the real Schur factorization of `A` on exit. It is an `intent(inout)`
      !!                argument.
      !!
      !!    - `dico`    :   `character(len=1)` input argument specifying whether the continuous
      !!                    (`dico = "c"`) or discrete (`dico = "d"`) Lyapunov equation is being
      !!                    solved. It is an `intent(in)` argument.
      !!
      !!    - `op`  :   `character(len=1)` input specifying the form of \(\mathrm{op}(A)\) to
      !!                be used as follows:
      !!                - `"n"`: \(\mathrm{op}(A) = A\) (no transpose).
      !!                - `"t"`: \(\mathrm{op}(A) = \A^{\top}\) (transpose).
      !!
      !!    - `factorized`  :   Boolean flag determining whether the input matrix `A` has
      !!                        already been factorized or not. It is an `intent(in)` argument.
      !!
      !!    - `job` :   `character(len=1)` input specifying the computation to be performed as
      !!                follows:
      !!                - `"x"`: compute the solution only.
      !!                - `"s"`: compute the separation only.
      !!                - `"b"`: compute both the solution and the separation.
      !!
      !!    - `scale` (optional)    :   double precision floating point number returned by
      !!                                `sb03md`. It is set between 0 and 1 to prevent the
      !!                                solution overflowing. It is an `intent(out)` argument.
      !!
      !!    - `separation` (optional)   :   double precision floating point number returned by
      !!                                    `sb03md` if `job = "s"` or `job = "b"`. It is an
      !!                                    `intent(out)` argument.
      !!
      !!    - `ferr` (optional) :   double precision floating point number returned by `sb03md`.
      !!                            It is an estimated forward error bound for the solution `X`.
      !!                            If `job = "x"` or `job ="s"`, it is not referenced.
      !!
      !!    - `wr`, `wi` (optional) :   double precision rank 1 array of size `n`.
      !!                                If `factorized = .false.`, they contain on exit the
      !!                                real and imaginary parts of the eigenvalues of `A`.
      !!
      !!    - `iwork` (optional)    :   integer array of dimension `n**2`. It is not referenced
      !!                                if `job = "x"`.
      !!
      !!    - `dwork` (optional)    :   double precision array of size `ldwork`, where `ldwork`
      !!                                is computed with the `lyapunov_workspace` function.
      pure module subroutine solve_lyapunov(A, C, U, dico, op, factorized, job, scale, separation, ferr, wr, wi, iwork, dwork)
         real(dp), intent(inout), target         :: A(:, :)
         real(dp), intent(inout), target         :: C(:, :)
         real(dp), intent(inout), target         :: U(:, :)
         character(len=1), intent(in)            :: dico
         character(len=1), intent(in)            :: op
         logical, intent(in)                     :: factorized
         character(len=1), intent(in)            :: job
         real(dp), optional, intent(out)         :: scale
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
