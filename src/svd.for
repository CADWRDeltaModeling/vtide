


c***************************************************************************
      
      subroutine svd(matq,matu,matv,cov,w,p,b,ndef,sig,ic,m,n,mm,nn,toler,jc
     1              ,ssq,res,residuals)
c-----------------------------------------------------------------------
c svd() uses singular-value-decomposition to calculate the least-squares
c solution p to an overdetermined system of linear equations with 
c coefficient matrix q, which includes multiple right hand side vectors packed
c side by side into the matrix b.
c
c Much of the api for this routine was based on one by j. cherniawsky, august 1997,
c but some of it (particularly the contents of matq) are now quite different.
c The code basically wraps calls to lapack, which needs to be linked. Prior code was
c based on Numerical Recipes which has licensing issues (and isn't as good anyway).
c
c there are two ways to use svd:
c 1 given an overdetermined system, svd will orthogonalize
c   a and b and produce the least-squares solution.
c 2 given an orthogonalized a (i.e. output from 1),
c   svd will orthogonalize b with respect to a and produce
c   the least-squares solution. this allows the use of
c   multiple r.h.s. without reorthogonalizing a.
c
c No testing has been performed yet for the case ic = 2, since the prior code didn't really
c use it and multiple right hand sides can be handled all at once with a larger rhs matrix b of
c side-by-side right hand sides.
c
c This incarnation due to Eli Ateljevich, California Dept Water Resources
c-----------------------------------------------------------------------
      implicit none
      real*8,intent(in)    :: matq(m,n)   !< design matrix of the overdetermined problem
      real*8,intent(inout) :: matu(m,n)   !< u matrix from svd decomposition
      real*8,intent(inout) :: matv(n,n)   !< v matrix from sdv decomposition
      
      real*8,intent(out) :: cov(n,n)    !< covariance matrix from design matrix
      real*8,intent(out) :: w(nn)      !< an nn-diagonal vector(matrix) of n eigenvalues of q (u)
      real*8,intent(out) :: p(nn,ndef) !< an output array least-squares coefficients
      real*8,intent(in)  :: b(mm,ndef) !< matrix of ndef right hand sides of size mm     
      integer            :: ndef       !< number of right hand sides
      real*8,intent(in) ::  sig(mm)   !< measurement error (standard deviation) for the rhs from the calling program (can be set to 1.).
      integer,intent(in)::  ic         !< input code: 1 for fresh orthogonalization or 2 to re-use prior orthogonalization of q part
      integer,intent(in) :: m         !< number of equations in least squares system
      integer,intent(in) :: n         !< number of columns in least squares system
      integer,intent(in) :: mm        !< number of columns of q
      integer,intent(in) :: nn        !< number of rows of q
      real*8,intent(in)  :: toler     !< an input tolerance, for acceptable matrix condition number
      integer,intent(out):: jc       !< output code which is set to the index of 1st dependent column, if such a column is detected. !todo: still true?
      
      real*8 :: ssq(ndef)      !< sum squared residuals for each rhs     
      real*8 :: res(ndef)      !< maximum residual in magnitude for each rhs
      real*8 :: residuals(m,ndef) !< residuals for each rhs
      
      
c ======== Locals      
      real*8, parameter :: eps = 1.d-10  !< threshold for labeling an singular value "very small" and zeroing it

      real*8 :: resi
      real*8 :: sum
      real*8 :: thresh
      !integer, parameter :: nwt=302 ! make nwt > expected value of nn   !todo: this was only used for wti. why was it needed? 
      real*8 :: wmax
      real*8 :: wti(nn) 
      integer info,lda, ldu, ldvt, lwork
      integer i,j,k,idef
      real*8  y(n)
      real*8  vt(n,n), work(max(1,3*min(m,n)+max(m,n),5*min(m,n)) )
      
c =====================================      
      
      ldu=m
      ldvt=n	
	lda=mm
	lwork=max(1,3*min(m,n)+max(m,n),5*min(m,n))
	
      jc=0

c no need to solve if only rhs has changed
      if(ic.eq.2) go to 10

c compute svd decomposition of matu(=a), with a being replaced by its upper
c matrix u, viz a=u*w*transpose(v), and vector w is output of a diagonal 
c matrix of singular values w(i), i=1,n.
c      call dsvdcmp(matu,m,n,mm,nn,w,matv)

      matu = matq(1:m,1:n)
      call dgesvd( 'O', 'A', M, N, matu, mm, w, matu, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )
      do i = 1,n
        do j = 1,n
           matv(i,j) = VT(j,i)
        end do
      end do

c check for small singular values
      wmax=0.
      do j=1,n
        if(w(j).gt.wmax) wmax=w(j)
      enddo
      thresh=toler*wmax
      do j=1,n
        if(w(j).lt.thresh) then
          w(j)=0.d0
          if(jc.lt.1) jc=j
        endif
      enddo

c compute summation weights (wti, used below)
  10  do j=1,n
        wti(j)=0.d0
        if(w(j).gt.eps) then
c         wti(j)=sig(j)*sig(j)/(w(j)*w(j))
          wti(j)=1.d0/(w(j)*w(j))
        endif
      enddo
c use back-substitution to compute the solution p(i), i=1,n
c     ! need to filter columns where w < thresh
      do idef=1,ndef
        call dgemv('T',m,n,1.d0,matu,m,b(:,idef),1,0.d0,y,1)
        y = y/w
        call dgemv('T',n,n,1.d0,vt,n,y,1,0.d0,p(1,idef),1)
      end do
      
c compute chisq (=ssq) and the largest residual (res)
      ssq=0.d0
      res=0.d0
      do idef=1,ndef
        do i=1,m
          sum=0.d0
          do j=1,n
            sum=sum+p(j,idef)*matq(i,j)
          enddo
          resi=b(i,idef)-sum
	    residuals(i,idef)=resi
          res(idef)=max(res(idef),abs(resi))
          ssq(idef)=ssq(idef)+resi**2
      end do
      end do
c compute variances, covariances, these may need to be given dimension
c of b(i), e.g., using sig(i), but this is better done after return to main
      do i=1,n
        do j=1,i
          sum=0.d0
          do k=1,n
            sum=sum+matv(i,k)*matv(j,k)*wti(k)
          enddo
          cov(i,j)=sum
          cov(j,i)=sum
        enddo
      enddo
      return
      end subroutine