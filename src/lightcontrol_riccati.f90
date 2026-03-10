submodule(lightcontrol) lightcontrol_riccati
   use slicot, only: sb02md, sb02mt
   implicit none(type, external)

contains

   !------------------------------------
   !-----     LOW-LEVEL DRIVER     -----
   !------------------------------------

   module procedure riccati_workspace
   integer :: ldwork_sb02mt, iwork(1), info, ipiv(1), oufact
   character(len=1), parameter :: uplo = "u"
   real(dp) :: a(1, 1), b(1, 1), q(1, 1), r(1, 1), l(1, 1), g(1, 1), dwork(1)
   associate (lda => n, ldb => n, ldq => merge(1, n, jobl == "z"), ldr => m, &
              ldl => merge(1, n, jobl == "z"), ldg => merge(1, n, jobg == "n"))
      ldwork = -1 ! Workspace query.
      call sb02mt(jobg, jobl, fact, uplo, n, m, a, lda, b, ldb, q, ldq, r, ldr, l, ldl, ipiv, &
                  oufact, g, ldg, iwork, dwork, ldwork, info)
      ldwork = int(ceiling(dwork(1), kind=dp))
      call assert(assertion=info == 0, &
                  description="Error during workspace query in sb02mt.")
      ! Maximum over ldwork for sb02mt and sb02md.
      ldwork = max(ldwork, 6*n)
   end associate
   end procedure riccati_workspace

   module procedure solve_riccati
   logical :: copy_a
   integer :: ldwork, info
   real(dp) :: rcond_
   real(dp), allocatable :: G(:, :)
   logical, pointer :: bwork_(:)
   integer, pointer :: iwork_(:)
   real(dp), pointer :: amat(:, :), dwork_(:), wr_(:), wi_(:)
   associate (lda => size(A, 1), n => size(A, 2), ldb => size(B, 1), m => size(B, 2), &
              ldq => size(Q, 1), ldr => size(R, 1), ldg => size(A, 2), &
              jobg => "g", jobl => "z", fact => "n", uplo => "u")
      !----- Input validation -----
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

      call assert(assertion=present(wr) .eqv. present(wi), &
                  description="Both wr and wi need to be passed.")
      if (present(wr)) then
         call assert(assertion=(size(wr) == 2*n) .and. (size(wi) == 2*n), &
                     description="wr and wi need to be of size 2*n")
         wr_ => wr
         wi_ => wi
      else
         allocate (wr_(2*n), wi_(2*n), source=0.0_dp)
      end if

      call assert(assertion=any(dico == ["c", "d"]), &
                  description="dico needs to be equal to 'c' or 'd'.")
      call assert(assertion=any(scale == ["g", "n"]), &
                  description="scale needs to be equal to 'g' or 'n'.")

      !----- Workspace validation -----
      if (present(iwork)) then
         iwork_ => iwork
      else
         allocate (iwork_(max(m, 2*n)), source=0)
      end if
      call assert(assertion=size(iwork_) >= max(m, 2*n), &
                  description="iwork needs to be of size 2*n.")

      if (present(bwork)) then
         bwork_ => bwork
      else
         allocate (bwork_(2*n), source=.false.)
      end if
      call assert(assertion=size(bwork_) == 2*n, &
                  description="bwork needs to be of size 2*n.")

      !> Workspace query.
      ldwork = riccati_workspace(size(B, 2), n, dico, jobg, jobl, fact)

      if (present(dwork)) then
         dwork_ => dwork
      else
         allocate (dwork_(ldwork), source=0.0_dp)
      end if
      call assert(assertion=size(dwork_) >= ldwork, &
                  description="Dimension of dwork is too small.")

      copy_a = .not. optval(overwrite_a, .false.)
      if (copy_a) then
         allocate (amat, source=a)
      else
         amat => a
      end if

      !----- Required allocations -----
      allocate (G(n, n), source=0.0_dp)

      !----- Preprocessing -----
      block
         integer, parameter :: ldl = 1
         integer :: ipiv(m), oufact
         real(dp) :: l(ldl, ldl)
         call sb02mt(jobg, jobl, fact, uplo, n, m, &
                     amat, lda, b, ldb, q, ldq, r, ldr, l, ldl, &
                     ipiv, oufact, g, ldg, iwork_, dwork_, ldwork, info)
         call assert(assertion=info == 0, &
                     description="Error in sb02mt.")
      end block

      ! ----- Solve Riccati equation -----
      block
         character(len=1), parameter :: hinv = "i", sort = "s"
         integer :: lds, ldu, i, j
         real(dp), allocatable :: s(:, :), u(:, :)

         !> Allocate arrays.
         lds = 2*n; ldu = 2*n
         allocate (s(2*n, 2*n), u(2*n, 2*n), source=0.0_dp)

         !> Riccati solver.
         call sb02md(dico, hinv, uplo, scale, sort, &
                     n, amat, lda, g, ldg, q, ldq, &
                     rcond_, wr_, wi_, s, lds, u, ldu, &
                     iwork_, dwork_, ldwork, bwork_, info)
         call assert(assertion=info == 0, &
                     description="Error in sb02md.")
      end block

      ! ----- Optional return values -----
      if (present(rcond)) rcond = rcond_
   end associate
   end procedure solve_riccati

   !----------------------------------------
   !-----     HIGH-LEVEL FUNCTIONS     -----
   !----------------------------------------

   module procedure care
   character(len=1), parameter :: dico = "c"
   character(len=1), parameter :: scale = "n"
   real(dp), allocatable :: amat(:, :), bmat(:, :), rmat(:, :)
   !> Prepare input for slicot.
   allocate (amat, source=a)
   allocate (bmat, source=b)
   allocate (rmat, source=r)
   allocate (x, source=q)
   !> Solve continuous-time algebraic Riccati equation.
   call solve_riccati(amat, bmat, x, rmat, dico, scale)
   end procedure care

   module procedure dare
   character(len=1), parameter :: dico = "d"
   character(len=1), parameter :: scale = "n"
   real(dp), allocatable :: amat(:, :), bmat(:, :), rmat(:, :)
   !> Prepare input for slicot.
   allocate (amat, source=a)
   allocate (bmat, source=b)
   allocate (rmat, source=r)
   allocate (x, source=q)
   !> Solve continuous-time algebraic Riccati equation.
   call solve_riccati(amat, bmat, x, rmat, dico, scale)
   end procedure dare

end submodule lightcontrol_riccati
