submodule(lightcontrol) lightcontrol_riccati
   implicit none(type, external)

   interface
      pure subroutine sb02md(dico, hinv, uplo, scal, sort, n, a, lda, g, ldg, q, ldq, &
                             rcond, wr, wi, s, lds, u, ldu, iwork, dwork, ldwork, bwork, info)
         import dp
         character(len=1), intent(in) :: dico, hinv, uplo, scal, sort
         integer, intent(in)          :: n, lda, ldg, ldq, lds, ldu
         integer, intent(inout)       :: ldwork
         real(dp), intent(inout)      :: a(lda, *), g(ldg, *), q(ldq, *)
         real(dp), intent(out)        :: rcond, wr(*), wi(*), s(lds, *), u(ldu, *), dwork(*)
         integer, intent(out)         :: iwork(*), info
         logical, intent(out)         :: bwork(*)
      end subroutine sb02md

      pure subroutine sb02mt(jobg, jobl, fact, uplo, n, m, a, lda, b, ldb, q, ldq, &
                             r, ldr, l, ldl, ipiv, oufact, g, ldg, iwork, dwork, ldwork, info)
         import dp
         character(len=1), intent(in) :: jobg, jobl, fact, uplo
         integer, intent(in) :: n, m, lda, ldb, ldq, ldr, ldl, ldg
         real(dp), intent(inout) :: a(lda, *), b(ldb, *), q(ldq, *)
         real(dp), intent(inout) :: r(ldr, *), l(ldl, *), g(ldg, *)
         integer, intent(inout) :: ipiv(*), ldwork
         integer, intent(out) :: oufact, info, iwork(*)
         real(dp), intent(out) :: dwork(*)
      end subroutine sb02mt
   end interface

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

end submodule lightcontrol_riccati
