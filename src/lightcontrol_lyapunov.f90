submodule(lightcontrol) lightcontrol_lyapunov
   use slicot, only: sb03md
   implicit none(type, external)

contains

   !------------------------------------
   !-----     LOW-LEVEL DRIVER     -----
   !------------------------------------

   module procedure lyapunov_workspace
   character(len=1), parameter :: trana = "n"
   real(dp)                    :: a(1, 1), u(1, 1), c(1, 1), wr(1), wi(1), dwork(1)
   real(dp)                    :: scale, sep, ferr
   integer                     :: iwork(1), info

   associate (lda => n, &   ! Leading dimension of A.
              ldu => n, &   ! Leading dimension of U (Schur orthogonal transformation).
              ldc => n)     ! Leading dimension C (rhs of Lyapunov equation).

      !> Workspace query.
      ldwork = -1
      call sb03md(dico, job, fact, trana, &     ! Setup.
                  n, a, lda, u, ldu, c, ldc, &  ! Matrices.
                  scale, sep, ferr, wr, wi, &   ! By-product variables.
                  iwork, dwork, ldwork, info)   ! Workspaces.
      !> Optimal workspace size.
      ldwork = int(ceiling(dwork(1), kind=dp))
      call assert(assertion=info == 0, &
                  description="Error during workspace query in sb03md.")
   end associate

   end procedure lyapunov_workspace

   module procedure solve_lyapunov
   character(len=1)  :: fact
   integer           :: ldwork, info
   real(dp)          :: scale_, sep, ferr_
   integer, pointer  :: iwork_(:)
   real(dp), pointer :: dwork_(:), wr_(:), wi_(:)

   associate (lda => size(A, 1), &  ! Leading dimension of A.
              n => size(A, 2), &    ! Number of states.
              ldc => size(C, 1), &  ! Leading dimension of C.
              ldu => size(U, 1))    ! Leading dimension of U.

      ! ----- Input validation -----
      !> Check matrices sizes.
      call assert(assertion=lda == n, &
                  description="Matrix A needs to be square.")
      call assert(assertion=ldc == size(C, 2), &
                  description="Matrix C needs to be square.")
      call assert(assertion=ldu == size(U, 2), &
                  description="Matrix U needs to be square.")

      call assert(assertion=lda == ldc, &
                  description="Dimensions of A and C are inconsistent.")
      call assert(assertion=lda == ldu, &
                  description="Dimensions of A and U are inconsistent.")

      !> Check setup.
      call assert(assertion=any(job == ["x", "s", "b"]), &
                  description="Invalid job. Needs to be 'x', 's', or 'b'.")
      call assert(assertion=any(op == ["n", "t", "c"]), &
                  description="Invalid op. Needs to be 'n', 't', or 'c'.")

      if (.not. factorized) then
         call assert(assertion=present(wr) .and. present(wi), &
                     description="Both wr and wi need to be passed if factorized = .false.")
         call assert(assertion=size(wr) == n, &
                     description="Dimension of wr is inconsistent with A.")
         call assert(assertion=size(wi) == n, &
                     description="Dimension of wi is inconsistent with A.")
         wi_ => wi
         wr_ => wr
      else
         allocate (wr_(1), wi_(1), source=0.0_dp)
      end if

      ! ----- Workspace validation -----
      !> Is A already factorized?
      fact = merge("f", "n", factorized)

      !> Real workspace.
      ldwork = lyapunov_workspace(n, dico, job, fact)   ! Workspace query.
      if (present(dwork)) then
         dwork_ => dwork
      else
         allocate (dwork_(ldwork), source=0.0_dp)
      end if
      call assert(assertion=size(dwork_) >= ldwork, &
                  description="Dimension of dwork is too small.")

      !> Integer workspace.
      if (present(iwork)) then
         iwork_ => iwork
      else
         allocate (iwork_(n**2), source=0)
      end if
      call assert(assertion=size(iwork_) >= n**2, &
                  description="Dimension of iwork is too small.")

      !----- Solve Lyapunov equation -----
      call sb03md(dico, job, fact, op, &            ! Setup.
                  lda, a, n, u, ldu, c, ldc, &      ! Matrices.
                  scale_, sep, ferr_, wr_, wi_, &   ! By-products.
                  iwork_, dwork_, ldwork, info)     ! Workspaces.

      !----- Optional returns -----
      if (present(scale)) scale = scale_
      if (present(separation)) separation = sep
      if (present(ferr)) ferr = ferr_
   end associate

   end procedure solve_lyapunov

   !----------------------------------------
   !-----     HIGH-LEVEL FUNCTIONS     -----
   !----------------------------------------

   module procedure lyap
   logical, parameter          :: factorized = .false.  ! A is not factorized.
   character(len=1), parameter :: dico = "c"            ! Continuous-time problem.
   character(len=1), parameter :: job = "x"             ! Compute only the solution.
   character(len=1), parameter :: op = "t"              ! A.T @ X + X @ A = C
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:) ! Work arrays.
   real(dp)                    :: scale                             ! Scaling for no-overflow.
   associate (n => size(A, 1))
      !> Prepare input for SLICOT.
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      allocate (x, source=q)
      !> Solve Lyapunov equation.
      call solve_lyapunov(amat, x, u, dico, op, factorized, job, wr=wr, wi=wi, scale=scale)
   end associate
   end procedure lyap

   module procedure dlyap
   logical, parameter          :: factorized = .false.  ! A is not factorized.
   character(len=1), parameter :: dico = "d"            ! Discrete-time problem.
   character(len=1), parameter :: job = "x"             ! Compute only the solution.
   character(len=1), parameter :: op = "t"              ! A.T @ X @ A - X = C.
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:) ! Work arrays.
   real(dp)                    :: scale                             ! Scaling for no-overflow.
   associate (n => size(A, 1))
      !> Prepare input for SLICOT.
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      allocate (x, source=q)
      !> Solve Lyapunov equation.
      call solve_lyapunov(amat, x, u, dico, op, factorized, job, wr=wr, wi=wi, scale=scale)
   end associate
   end procedure dlyap

   !----------------------------------------------------------
   !-----     CONTROLLABILITY/OBSERVABILITY GRAMIANS     -----
   !----------------------------------------------------------

   module procedure ctrb_gramian_siso
   real(dp), pointer :: bmat(:, :)
   bmat(1:size(b, 1), 1:1) => b
   p = ctrb_gramian_mimo(a, bmat, discrete)
   end procedure ctrb_gramian_siso

   module procedure ctrb_gramian_mimo
   logical, parameter          :: factorized = .false.  ! A is not factorized.
   character(len=1), parameter :: job = "x"             ! Compute only the solution.
   character(len=1), parameter :: op = "n"              ! A @ X + X @ A.T = -B @ B.T
   character(len=1)            :: dico                  ! Continuous or discrete time?
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:) ! Work arrays.
   real(dp)                    :: scale                             ! Scale for no overflow.

   associate (lda => size(a, 1), &  ! Leading dimension of A.
              n => size(a, 2))      ! Number of states.
      !> Input validation
      call assert(assertion=lda == n, &
                  description="Matrix A needs to be square.")
      call assert(assertion=size(b, 1) == n, &
                  description="Leading dimension of B inconsistent with A.")
      !> Prepare input for SLICOT
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      p = matmul(b, transpose(b))   ! NOTE: Replace with appropriate BLAS function.
      dico = "c"; if (present(discrete)) dico = merge("d", "c", discrete)
      !> Solve Lyapunov equation.
      call solve_lyapunov(amat, p, u, dico, op, factorized, job, wr=wr, wi=wi, scale=scale)
   end associate

   end procedure ctrb_gramian_mimo

   module procedure obs_gramian_siso
   real(dp), pointer :: cmat(:, :)
   cmat(1:1, 1:size(c, 1)) => c
   q = obs_gramian_mimo(a, cmat, discrete)
   end procedure obs_gramian_siso

   module procedure obs_gramian_mimo
   logical, parameter          :: factorized = .false.  ! A is not factorized.
   character(len=1), parameter :: job = "x"             ! Compute only the solution.
   character(len=1), parameter :: op = "t"              ! A.T @ X + X @ A = - C.T @ C
   character(len=1)            :: dico                  ! Continuous or discrete time?
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:) ! Work arrays.
   real(dp)                    :: scale                             ! Scale for no-overflow.

   associate (lda => size(a, 1), &  ! Leading dimension of A.
              n => size(a, 2))      ! Number of states.
      !> Input validation.
      call assert(assertion=lda == n, &
                  description="Matrix A needs to be square.")
      call assert(assertion=size(c, 2) == n, &
                  description="Number of columns of C inconsistent with A.")
      !> Prepare input for SLICOT.
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      q = matmul(transpose(c), c)   ! NOTE: Replace with appropriate BLAS function.
      dico = "c"; if (present(discrete)) dico = merge("d", "c", discrete)
      !> Solve Lyapunov equation.
      call solve_lyapunov(amat, q, u, dico, op, factorized, job, wr=wr, wi=wi, scale=scale)
   end associate

   end procedure obs_gramian_mimo

end submodule lightcontrol_lyapunov
