C Static fortran function(s) for loading in R
C Copyright Metasys, INSA/INRA UMR 792, Toulouse, France.
C Author: Serguei SOKOL
C Date: 2010-10-17

      subroutine f2b(bfpr, n, fwrv, incu, fb, b)
C Calculate rhs term b (dense vector) for cumomer system Ax=b
C and the sum of fluxes fb to be included in the diagonal
C term of the matrix A
C bfpr is 4-row index matrix: irow, ifl, ix1, ix2
C n column number in bfpr
C fwrv is a flux vector
C incu is a composite cumomer vector c(1,xinput,xcumo)
      IMPLICIT NONE
      integer i, ir, bfpr(4, n), n
      double precision f, fwrv(1), incu(1), fb(1), b(1), term
      do 1 i=1,n
         ir=bfpr(1,i)
         f=fwrv(bfpr(2,i))
         fb(ir)=fb(ir)+f
         term=f*incu(bfpr(3,i))*incu(bfpr(4,i))
         b(ir)=b(ir)-term
C         write(0,*) "i,ir,f,fb(ir),b(ir),incu1,incu2,term=",
C     & i, ir, f, fb(ir), b(ir),incu(bfpr(3,i)),incu(bfpr(4,i)),term
    1 continue
      return
      end subroutine f2b

      subroutine f2a(fi, islot, n, fwrv, x, nx)
C Calculate off-diagonal terms for a cumomer system Ax=b
C fi is 2-row index matrix: rep, ifl
C n column number in fi
C fwrv is a flux vector
C x is the resulting slot in A@x
      IMPLICIT NONE
      integer i, rep, nrep, ii, bai, fi(2, n), islot(1), n, nx
      double precision fwrv(1), x(nx)
C     counter in i-slot
      ii=0
      i=1
C      write(0,*) "fi,islot,n,nx,fwrv,x=",fi,"; ",islot,n,nx,fwrv,x
      do while (i .le. n)
C        run through repetitions in islot
         ii=ii+1
         nrep=fi(1,i)
         do rep=1,nrep
            x(ii)=x(ii)+fwrv(fi(2,i))
C            write(0,*) "i,ii,x(ii),fi(2,i),fwrv(fi(2,i))=",
C     &i,ii,x(ii),fi(2,i),fwrv(fi(2,i))
            i=i+1
         enddo
      enddo
      if (ii .ne. nx) stop "f2a: ii!=nx at the end of loop"
      return
      end subroutine f2a

      subroutine x2a_fx(xi, n, islot, x0, x, nx)
C Calculate dA_df*x (@x slot)
C xi is 2-row index matrix: ix, idiag
C n column number in xi
C islot is a_fx@i
C x0 is c(0, incu)
C x is the resulting slot in a_fx@x
      IMPLICIT NONE
      integer i, ii, rep, nrep, xi(3, n), n, islot(1), nx
      double precision x0(1), x(nx)
C      write(0,*) "xi,n,islot,x0,x=",xi,"; ",n,islot,x0,x
      i=1
      ii=0
      do while (i .le. n)
         ii=ii+1
         nrep=xi(1,i)
         do rep=1,nrep
            x(ii)=x(ii)+x0(xi(2,i))-x0(xi(3,i))
            i=i+1
         enddo
C         write(0,*) "xi(1,i),xi(2,i),x0(i)=",xi(1,i),xi(2,i),x0(i)
C         write(0,*) "i,ifl,ii,x0(xi(1,i)),x0(xi(2,i)),x(ii)=",
C     &i,ifl,ii,x0(xi(1,i)),x0(xi(1,i)),x(ii)
      enddo
      if (ii .ne. nx) stop "f2a: ii!=nx at the end of loop"
      return
      end subroutine x2a_fx

      subroutine x2b_f(xi, n, islot, incu, x, nx)
C Calculate db_df (@x slot)
C xi is 2-row index matrix: ix1, ix2
C n column number in xi
C islot is b_f@i
C incu is c(1, xinp, xcumo)
C x is the resulting slot in b_f@x
      IMPLICIT NONE
      integer i, ifl, xi(2, n), n, islot(1), nx
      double precision incu(1), x(nx)
      if (n .ne. nx) stop "x2b_f: n!=nx"
      do 1 i=1,n
         x(i)=x(i)-incu(xi(1,i))*incu(xi(2,i))
    1 continue
      return
      end subroutine x2b_f

      subroutine fx2b_x(fxi, n, islot, fwrv, incu, x, nx)
C Calculate db_df (@x slot)
C fxi is 3-row index matrix: rep, ifl, icof
C n column number in fxi
C islot is b_x@i
C incu is c(1, xinp, xcumo)
C x is the resulting slot in b_x@x
      IMPLICIT NONE
      integer i, ii, rep, nrep, fxi(3, n), n, islot(1), nx
      double precision fwrv(1), incu(1), x(nx)
      i=1
      ii=0
      do while (i .le. n)
         ii=ii+1
         nrep=fxi(1,i)
         do rep=1, nrep
            x(ii)=x(ii)-fwrv(fxi(2,i))*incu(fxi(3,i))
            i=i+1
         enddo
      enddo
      if (ii .ne. nx) stop "fx2b_x: ii!=nx at the end of loop"
      return
      end subroutine fx2b_x

