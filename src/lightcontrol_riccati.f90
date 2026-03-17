submodule(lightcontrol) lightcontrol_riccati
   use slicot, only: sb02md, sb02mt
   implicit none(type, external)

contains

   !------------------------------------
   !-----     LOW-LEVEL DRIVER     -----
   !------------------------------------

   module procedure riccati_workspace
   integer                     :: ldwork_sb02mt, iwork(1), info, ipiv(1), oufact
   character(len=1), parameter :: uplo = "u"
   real(dp)                    :: a(1, 1), b(1, 1), &   ! State-space model.
                                  q(1, 1), r(1, 1), &   ! LQR cost.
                                  l(1, 1), g(1, 1), &   ! Returned arrays.
                                  dwork(1)              ! Real workspace.

   associate (lda => n, &                           ! Leading dimension of A.
              ldb => n, &                           ! Leading dimension of B.
              ldq => merge(1, n, jobl == "z"), &    ! Leading dimension of Q.
              ldr => m, &                           ! Leading dimension of R.
              ldl => merge(1, n, jobl == "z"), &    ! Leading dimension of L.
              ldg => merge(1, n, jobg == "n"))      ! Leading dimension of G.

      !> Workspace query.
      ldwork = -1
      call sb02mt(jobg, jobl, fact, uplo, &         ! Setup.
                  n, m, a, lda, b, ldb, &           ! State-space model.
                  q, ldq, r, ldr, &                 ! LQR cost.
                  l, ldl, ipiv, oufact, g, ldg, &   ! Returned arrays.
                  iwork, dwork, ldwork, info)       ! Workspaces.
      !> Optimal workspace size.
      ldwork = int(ceiling(dwork(1), kind=dp))
      call assert(assertion=info == 0, &
                  description="Error during workspace query in sb02mt.")
      ! Maximum over ldwork for sb02mt and sb02md.
      ldwork = max(ldwork, 6*n)
   end associate

   end procedure riccati_workspace

   module procedure solve_riccati
   logical               :: copy_a
   integer               :: ldwork, info
   real(dp)              :: rcond_
   real(dp), allocatable :: G(:, :)
   logical, pointer      :: bwork_(:)
   integer, pointer      :: iwork_(:)
   real(dp), pointer     :: amat(:, :), dwork_(:), wr_(:), wi_(:)

   associate (lda => size(A, 1), &  ! Leading dimension of A.
              n => size(A, 2), &    ! Number of states.
              ldb => size(B, 1), &  ! Leading dimension of B.
              m => size(B, 2), &    ! Number of inputs.
              ldq => size(Q, 1), &  ! Leading dimension of Q.
              ldr => size(R, 1), &  ! Leading dimension of R.
              ldg => size(A, 2), &  ! Leading dimension of G = B @ inv(R) @ B.T
              jobg => "g", &        ! G = B @ inv(R) @ B.T needs to be computed.
              jobl => "z", &        ! L is the zero matrix.
              fact => "n", &        ! R is not factorized on entry.
              uplo => "u")          ! Use upper triangular format for symmetric matrices.

      !----- Input validation -----
      !> Check matrices sizes.
      call assert(assertion=lda == n, &
                  description="Matrix A needs to be square.")
      call assert(assertion=lda == ldb, &
                  description="Leading dimensions of A and B are inconsistent.")
      call assert(assertion=all(shape(a) == shape(q)), &
                  description="A and Q need to have the same dimensions.")
      call assert(assertion=size(R, 1) == size(R, 2), &
                  description="R needs to be a square matrix.")
      call assert(assertion=size(R, 1) == size(B, 2), &
                  description="Dimension of R is inconsistent with the number of inputs in B.")

      !> Eigenvalue array size.
      call assert(assertion=present(wr) .eqv. present(wi), &
                  description="Both wr and wi need to be passed.")
      if (present(wr)) then
         wr_ => wr
         wi_ => wi
      else
         allocate (wr_(2*n), wi_(2*n), source=0.0_dp)
      end if
      call assert(assertion=(size(wr_) == 2*n) .and. (size(wi_) == 2*n), &
                  description="wr and wi need to be of size 2*n")

      !> Setup parameters.
      call assert(assertion=any(dico == ["c", "d"]), &
                  description="dico needs to be equal to 'c' or 'd'.")
      call assert(assertion=any(scale == ["g", "n"]), &
                  description="scale needs to be equal to 'g' or 'n'.")

      !----- Workspace validation -----
      !> Integer workspace.
      if (present(iwork)) then
         iwork_ => iwork
      else
         allocate (iwork_(max(m, 2*n)), source=0)
      end if
      call assert(assertion=size(iwork_) >= max(m, 2*n), &
                  description="iwork needs to be of size 2*n.")

      !> Logical workspace.
      if (present(bwork)) then
         bwork_ => bwork
      else
         allocate (bwork_(2*n), source=.false.)
      end if
      call assert(assertion=size(bwork_) == 2*n, &
                  description="bwork needs to be of size 2*n.")

      !> Real workspace.
      ldwork = riccati_workspace(size(B, 2), n, dico, jobg, jobl, fact) ! Query optimal size.
      if (present(dwork)) then
         dwork_ => dwork
      else
         allocate (dwork_(ldwork), source=0.0_dp)
      end if
      call assert(assertion=size(dwork_) >= ldwork, &
                  description="Dimension of dwork is too small.")

      !> Performances: overwrite A ?
      copy_a = .not. optval(overwrite_a, .false.)
      if (copy_a) then
         allocate (amat, source=a)
      else
         amat => a
      end if

      !----- Preprocessing -----
      allocate (G(n, n), source=0.0_dp)
      preprocessing: block
         integer, parameter :: ldl = 1
         integer            :: ipiv(m), oufact
         real(dp)           :: l(ldl, ldl)
         !> Compute G = B @ inv(R) @ B.T
         call sb02mt(jobg, jobl, fact, uplo, &          ! Setup.
                     n, m, amat, lda, b, ldb, &         ! State-space model.
                     q, ldq, r, ldr, &                  ! LQR cost.
                     l, ldl, ipiv, oufact, g, ldg, &    ! Working arrays.
                     iwork_, dwork_, ldwork, info)      ! Workspaces.
         call assert(assertion=info == 0, &
                     description="Error in sb02mt.")
      end block preprocessing

      ! ----- Solve Riccati equation -----
      solve: block
         character(len=1), parameter :: hinv = "i"  ! Inverse symplectic matrix computation.
         character(len=1), parameter :: sort = "s"  ! Stable eigenvalues first in Schur.
         integer                     :: lds, ldu
         real(dp), allocatable       :: s(:, :), u(:, :)
         !> Allocate arrays.
         lds = 2*n; ldu = 2*n
         allocate (s(2*n, 2*n), u(2*n, 2*n), source=0.0_dp)
         !> Riccati solver.
         call sb02md(dico, hinv, uplo, scale, sort, &       ! Setup.
                     n, amat, lda, g, ldg, q, ldq, &        ! State space model + LQR cost.
                     rcond_, wr_, wi_, s, lds, u, ldu, &    ! Working arrays.
                     iwork_, dwork_, ldwork, bwork_, info)  ! Workspace.
         call assert(assertion=info == 0, &
                     description="Error in sb02md.")
      end block solve

      ! ----- Optional return values -----
      if (present(rcond)) rcond = rcond_
   end associate

   end procedure solve_riccati

   !----------------------------------------
   !-----     HIGH-LEVEL FUNCTIONS     -----
   !----------------------------------------
   module procedure care_siso
   real(dp), pointer :: bmat(:, :)
   real(dp)          :: rmat(1, 1)
   bmat(1:size(b), 1:1) => b
   rmat = r
   X = care_mimo(a, bmat, q, rmat)
   end procedure care_siso

   module procedure care_mimo
   character(len=1), parameter :: dico = "c"    ! Continuous-time.
   character(len=1), parameter :: scal = "n"    ! No scaling.
   real(dp), allocatable       :: amat(:, :), bmat(:, :), rmat(:, :)
   !> Prepare input for slicot.
   allocate (amat, source=a)
   allocate (bmat, source=b)
   allocate (rmat, source=r)
   allocate (x, source=q)
   !> Solve continuous-time algebraic Riccati equation.
   call solve_riccati(amat, bmat, x, rmat, dico, scal, overwrite_a=.true.)
   end procedure care_mimo

   module procedure dare_siso
   real(dp), pointer :: bmat(:, :)
   real(dp)          :: rmat(1, 1)
   bmat(1:size(b), 1:1) => b
   rmat = r
   X = dare_mimo(a, bmat, q, rmat)
   end procedure dare_siso

   module procedure dare_mimo
   character(len=1), parameter :: dico = "d"    ! Discrete-time.
   character(len=1), parameter :: scal = "n"    ! No scaling.
   real(dp), allocatable       :: amat(:, :), bmat(:, :), rmat(:, :)
   !> Prepare input for slicot.
   allocate (amat, source=a)
   allocate (bmat, source=b)
   allocate (rmat, source=r)
   allocate (x, source=q)
   !> Solve continuous-time algebraic Riccati equation.
   call solve_riccati(amat, bmat, x, rmat, dico, scal, overwrite_a=.true.)
   end procedure dare_mimo

end submodule lightcontrol_riccati
