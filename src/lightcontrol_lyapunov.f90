submodule(lightcontrol) lightcontrol_lyapunov
   implicit none(type, external)

   interface
      ! Lyapunov solver from SLICOT.
      pure subroutine sb03md(dico, job, fact, trana, n, a, lda, u, ldu, c, ldc, &
                             scale, sep, ferr, wr, wi, iwork, dwork, ldwork, info)
         import dp
         character(len=1), intent(in) :: dico, fact, job, trana
         integer, intent(in)          :: n, lda, ldu, ldc, ldwork
         integer, intent(out)         :: iwork(*), info
         real(dp), intent(inout)      :: a(lda, *), u(ldu, *), c(ldc, *)
         real(dp), intent(in)         :: scale
         real(dp), intent(out)        :: sep, ferr, wr(*), wi(*), dwork(*)
      end subroutine sb03md
   end interface

contains

   !------------------------------------
   !-----     LOW-LEVEL DRIVER     -----
   !------------------------------------

   module procedure lyapunov_workspace
   character(len=1), parameter :: trana = "n"
   real(dp) :: a(1, 1), u(1, 1), c(1, 1), wr(1), wi(1), dwork(1)
   real(dp) :: scale, sep, ferr
   integer  :: iwork(1), info

   ! ----- Instructions -----
   associate (lda => n, ldu => n, ldc => n)
      ldwork = -1 ! Workspace query.
      call sb03md(dico, job, fact, trana, n, a, lda, u, ldu, c, ldc, &
                  scale, sep, ferr, wr, wi, iwork, dwork, ldwork, info)
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
   associate (lda => size(A, 1), n => size(A, 2), ldc => size(C, 1), ldu => size(U, 1))
      ! ----- Input validation -----
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

      call assert(assertion=any(job == ["x", "s", "b"]), &
                  description="Invalid job. Needs to be 'x', 's', or 'b'.")
      call assert(assertion=any(op == ["n", "t", "c"]), &
                  description="Invalid op. Needs to be 'n', 't', or 'c'.")

      scale_ = optval(scale, 1.0_dp)
      call assert(assertion=(0.0_dp < scale_) .and. (scale_ <= 1.0_dp), &
                  description="scale needs to be between 0 (excluded) and 1 (included).")

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
      fact = merge("f", "n", factorized)
      !> Workspace query.
      ldwork = lyapunov_workspace(n, dico, job, fact)

      if (present(dwork)) then
         dwork_ => dwork
      else
         allocate (dwork_(ldwork), source=0.0_dp)
      end if
      call assert(assertion=size(dwork_) >= ldwork, &
                  description="Dimension of dwork is too small.")

      if (present(iwork)) then
         iwork_ => iwork
      else
         allocate (iwork_(n**2), source=0)
      end if
      call assert(assertion=size(iwork_) >= n**2, &
                  description="Dimension of iwork is too small.")

      !----- Solve Lyapunov equation -----
      call sb03md(dico, job, fact, op, lda, a, n, u, ldu, c, ldc, &
                  scale_, sep, ferr_, wr_, wi_, iwork_, dwork_, ldwork, info)
      !----- Optional returns -----
      if (present(separation)) separation = sep
      if (present(ferr)) ferr = ferr_
   end associate
   end procedure solve_lyapunov

   !----------------------------------------
   !-----     HIGH-LEVEL FUNCTIONS     -----
   !----------------------------------------

   module procedure lyap
   logical, parameter          :: factorized = .false.
   character(len=1), parameter :: dico = "c", job = "x", op = "t"
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:)
   associate (n => size(A, 1))
      ! ----- Prepare input for slicot -----
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      allocate (x, source=q)
      ! ----- Solve Lyapunov equation -----
      call solve_lyapunov(amat, x, u, dico, op, factorized, job, wr=wr, wi=wi)
   end associate
   end procedure lyap

   module procedure dlyap
   logical, parameter          :: factorized = .false.
   character(len=1), parameter :: dico = "d", job = "x", op = "t"
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:)
   associate (n => size(A, 1))
      ! ----- Prepare input for slicot -----
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      allocate (x, source=q)
      ! ----- Solve Lyapunov equation -----
      call solve_lyapunov(amat, x, u, dico, op, factorized, job, wr=wr, wi=wi)
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
   logical, parameter          :: factorized = .false.
   character(len=1), parameter :: job = "x", op = "n"
   character(len=1)            :: dico
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:)
   associate (lda => size(a, 1), n => size(a, 2))
      !----- Input validation ------
      call assert(assertion=lda == n, &
                  description="Matrix A needs to be square.")
      call assert(assertion=size(b, 1) == n, &
                  description="Leading dimension of B inconsistent with A.")
      ! ----- Prepare input for slicot -----
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      p = matmul(b, transpose(b))
      dico = "c"; if (present(discrete)) dico = merge("d", "c", discrete)
      ! ----- Solve Lyapunov equation -----
      call solve_lyapunov(amat, p, u, dico, op, factorized, job, wr=wr, wi=wi)
   end associate
   end procedure ctrb_gramian_mimo

   module procedure obs_gramian_siso
   real(dp), pointer :: cmat(:, :)
   cmat(1:1, 1:size(c, 1)) => c
   q = obs_gramian_mimo(a, cmat, discrete)
   end procedure obs_gramian_siso

   module procedure obs_gramian_mimo
   logical, parameter          :: factorized = .false.
   character(len=1), parameter :: job = "x", op = "t"
   character(len=1)            :: dico
   real(dp), allocatable       :: amat(:, :), u(:, :), wr(:), wi(:)
   associate (lda => size(a, 1), n => size(a, 2))
      !----- Input validation ------
      call assert(assertion=lda == n, &
                  description="Matrix A needs to be square.")
      call assert(assertion=size(c, 2) == n, &
                  description="Number of columns of C inconsistent with A.")
      ! ----- Prepare input for slicot -----
      allocate (amat, source=a)
      allocate (u(n, n), wr(n), wi(n), source=0.0_dp)
      q = matmul(transpose(c), c)
      dico = "c"; if (present(discrete)) dico = merge("d", "c", discrete)
      ! ----- Solve Lyapunov equation -----
      call solve_lyapunov(amat, q, u, dico, op, factorized, job, wr=wr, wi=wi)
   end associate
   end procedure obs_gramian_mimo

end submodule lightcontrol_lyapunov
