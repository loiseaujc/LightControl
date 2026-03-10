module LightControl
   use assert_m, only: assert => assert_always
   use stdlib_kinds, only: dp
   use stdlib_optval, only: optval
   implicit none(type, external)
   private

   !> Lyapunov equations.
   public :: lyap, dlyap, solve_lyapunov, lyapunov_workspace
   public :: ctrb_gramian, obs_gramian

   !> Riccati equations.
   public :: dare, care, solve_riccati, riccati_workspace

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
      module function lyap(A, Q) result(X)
         implicit none(type, external)
         real(dp), intent(in)  :: A(:, :)
         !> Dynamics matrix of dimension `n x n`.
         real(dp), intent(in)  :: Q(:, :)
         !> Symmetric matrix of dimension `n x n` corresponding to the right-hand
         !> side of the Lyapunov equation.
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
      module function dlyap(A, Q) result(X)
         implicit none(type, external)
         real(dp), intent(in)  :: A(:, :)
         !> Dynamics matrix of dimension `n x n`.
         real(dp), intent(in)  :: Q(:, :)
         !> Symmetrix matrix of dimension `n x n` corresponding to the right-hande
         !> side of the Lyapunov equation.
         real(dp), allocatable :: X(:, :)
         !> Solution of the Lyapunov equation, dimension `n x n`.
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
         implicit none(type, external)
         real(dp), intent(in)                     :: A(:, :)
         !> Dynamics matrix of dimension `n x n`.
         real(dp), intent(in), contiguous, target :: b(:)
         !> Input-to-state matrix for a single-input system, dimension `n`.
         logical, intent(in), optional            :: discrete
         !> Whether the system is discrete-time (`discrete = .true.`) or
         !> continuous-time (`discrete = .false.`) one. Default is `discrete = .false.`
         real(dp), allocatable                    :: P(:, :)
         !> Controllability gramian of the system, dimension `n x n`.
      end function ctrb_gramian_siso

      module function ctrb_gramian_mimo(A, B, discrete) result(P)
         implicit none(type, external)
         real(dp), intent(in)          :: A(:, :)
         !> Dynamics matrix of dimension `n x n`.
         real(dp), intent(in)          :: B(:, :)
         !> Input-to-state matrix for a multiple-input system, dimension `n x k`.
         logical, intent(in), optional :: discrete
         !> Whether the system is discrete-time (`discrete = .true.`) or
         !> continuous-time (`discrete = .false.`) one. Default is `discrete = .false.`
         real(dp), allocatable         :: P(:, :)
         !> Controllability gramian of the system, dimension `n x n`.
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
         implicit none(type, external)
         real(dp), intent(in)                     :: A(:, :)
         !> Dynamics matrix of the system, dimension `n x n`.
         real(dp), intent(in), contiguous, target :: c(:)
         !> Measurement matrix of the system for single-output system, dimension `n`.
         logical, intent(in), optional            :: discrete
         !> Whether the system is discrete-time (`discrete = .true.`) or
         !> continuous-time (`discrete = .false.`) one. Default is `discrete = .false.`
         real(dp), allocatable                    :: Q(:, :)
         !> Observablity gramian of the system, dimension `n x n`.
      end function obs_gramian_siso

      module function obs_gramian_mimo(A, C, discrete) result(Q)
         implicit none(type, external)
         real(dp), intent(in)          :: A(:, :)
         !> Dynamics matrix of the system, dimension `n x n`.
         real(dp), intent(in)          :: C(:, :)
         !> Measurement matrix of the system for multiple-output system, dimension `k x n`.
         logical, intent(in), optional :: discrete
         !> Whether the system is discrete-time (`discrete = .true.`) or
         !> continuous-time (`discrete = .false.`) one. Default is `discrete = .false.`
         real(dp), allocatable         :: Q(:, :)
         !> Observablity gramian of the system, dimension `n x n`.
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
      module subroutine solve_lyapunov(A, C, U, &
                                       dico, op, factorized, job, &
                                       scale, separation, ferr, wr, wi, iwork, dwork)
         implicit none(type, external)
         real(dp), intent(inout)                 :: A(:, :)
         !> Dynamics matrix of the system, dimension `n x n`.
         real(dp), intent(inout)                 :: C(:, :)
         !> On entry: right-hand side symmetric matrix of dimension `n x n`.
         !> On exit: Overwritten by the solution to the Lyapunov equation.
         real(dp), intent(inout)                 :: U(:, :)
         !> On entry: If `factorized = .true.`, it must contains the orthogonal matrix of the
         !>           real Schur form of `A` and is unchanged on exit. If `factorized = .false.`
         !>           its values on entry are irrelevant.
         !> On exit: If `factorized = .true.`, it is unchanged. If `factorized = .false.`, it
         !>          will contain the orthgonal matrix of the real Schur factorization of `A`.
         character(len=1), intent(in)            :: dico
         !> Whether the continuous (`dico = "c"`) or discrete (`dico = "d"`) time Lyapunov
         !> is being solved.
         character(len=1), intent(in)            :: op
         !> Whether `op(A) = A` (`op = "n"`) or `op(A) = transpose(A)` (`op = "t"`).
         logical, intent(in)                     :: factorized
         !> Whether `A` has already been factorized to its real Schur form.
         !> If `factorized = .true.`, `U` must contain the associated orthogonal transformation.
         character(len=1), intent(in)            :: job
         !> Input specifying the computation to be performed as follows:
         !>                - `"x"`: compute the solution only.
         !>                - `"s"`: compute the separation only.
         !>                - `"b"`: compute both the solution and the separation.
         real(dp), optional, intent(out)         :: scale
         !> double precision floating point number returned by `sb03md`. It is set
         !> between 0 and 1 to prevent the solution overflowing.
         real(dp), optional, intent(out)         :: separation
         !> double precision floating point number returned by `sb03md` if `job = "s"` or
         !> `job = "b"`.
         real(dp), optional, intent(out)         :: ferr
         !> double precision floating point number returned by `sb03md`. It is an estimated
         !> forward error bound for the solution `X`. If `job = "x"` or `job ="s"`,
         !> it is not referenced.
         real(dp), optional, intent(out), target :: wr(:), wi(:)
         !> If `factorized = .false.`, they contain on exit the real and imaginary parts
         !> of the eigenvalues of `A`.
         integer, optional, intent(out), target  :: iwork(:)
         !> Integer workspace array.
         real(dp), optional, intent(out), target :: dwork(:)
         !> double precision workspace array.
      end subroutine solve_lyapunov
   end interface

   interface
      module integer function lyapunov_workspace(n, dico, job, fact) result(ldwork)
         implicit none(type, external)
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

   !-------------------------------------
   !-----     RICCATI EQUATIONS     -----
   !-------------------------------------

   interface
      module integer function riccati_workspace(m, n, dico, jobg, jobl, fact) result(ldwork)
         implicit none(type, external)
         integer, intent(in) :: n
         !> Number of states.
         integer, intent(in) :: m
         !> Number of inputs.
         character(len=1), intent(in) :: dico
         character(len=1), intent(in) :: jobg
         character(len=1), intent(in) :: jobl
         character(len=1), intent(in) :: fact
      end function riccati_workspace
   end interface

   interface
      module subroutine solve_riccati(A, B, Q, R, dico, scale, &
                                      rcond, wr, wi, iwork, dwork, bwork, overwrite_a)
         implicit none(type, external)
         real(dp), intent(inout), target         :: A(:, :)
         real(dp), intent(inout)                 :: B(:, :)
         real(dp), intent(inout)                 :: Q(:, :)
         real(dp), intent(inout)                 :: R(:, :)
         character(len=1), intent(in)            :: dico
         character(len=1), intent(in)            :: scale
         real(dp), optional, intent(out)         :: rcond
         real(dp), optional, intent(out), target :: wr(:), wi(:)
         integer, optional, intent(out), target  :: iwork(:)
         real(dp), optional, intent(out), target :: dwork(:)
         logical, optional, intent(out), target  :: bwork(:)
         logical, optional, intent(in)           :: overwrite_a
      end subroutine solve_riccati
   end interface

   interface care
      module function care_siso(A, b, Q, r) result(X)
         implicit none(type, external)
         real(dp), intent(in) :: A(:, :)
         real(dp), intent(in), contiguous, target :: b(:)
         real(dp), intent(in) :: Q(:, :)
         real(dp), intent(in) :: r
         real(dp), allocatable :: X(:, :)
      end function care_siso

      module function care_mimo(A, B, Q, R) result(X)
         implicit none(type, external)
         real(dp), intent(in) :: A(:, :)
         real(dp), intent(in) :: B(:, :)
         real(dp), intent(in) :: Q(:, :)
         real(dp), intent(in) :: R(:, :)
         real(dp), allocatable :: X(:, :)
      end function care_mimo
   end interface care

   interface dare
      module function dare_siso(A, b, Q, r) result(X)
         implicit none(type, external)
         real(dp), intent(in) :: A(:, :)
         real(dp), intent(in), contiguous, target :: b(:)
         real(dp), intent(in) :: Q(:, :)
         real(dp), intent(in) :: r
         real(dp), allocatable :: X(:, :)
      end function dare_siso

      module function dare_mimo(A, B, Q, R) result(X)
         implicit none(type, external)
         real(dp), intent(in) :: A(:, :)
         real(dp), intent(in) :: B(:, :)
         real(dp), intent(in) :: Q(:, :)
         real(dp), intent(in) :: R(:, :)
         real(dp), allocatable :: X(:, :)
      end function dare_mimo
   end interface dare
end module LightControl
