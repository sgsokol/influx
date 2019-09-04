C These subrutines are used in influx_i to solve a series of problems
C 
C Copyright 2014, INRA, France.
C Author: Serguei SOKOL
C Date: 2014-07-22
      subroutine mult_bxt(bx_x, bx_i, bx_p, nr_bx, ntico, nc_bx, c,
     * ldc, nc_c, a)
C Calculate a_t=bx_t%*%c_t where a, bx, c are matrices, t is index in time
C and %*% is a dot product.
C bx is a sparse matrix of size (nr_bx*ntico, nc_bx) given by its fields
C _x, _i, and _p describing column wise storage. The length of _x and _i is nb_x (virtual parameter)
C the size of c which is dense matrix is (ldc*ntico, nc_c) and those
C of a (also dense) is (nr_bx, nc_c, ntico).
C The parameter ldc may be different from ncol(bx) (if it is, then ldc > ncol(bx)
C a is supposed to be initialized to zero
C The result is added to a so it must be initialized to 0 before call
C if a pure multiplication result is needed.
      IMPLICIT NONE
      integer bx_i(1), bx_p(nc_bx+1), nb_x, nr_bx, ntico, nc_bx,
     * ldc, nc_c
      integer icc, icb, ii, irb, it, irb1
      double precision bx_x(1), c(ldc*ntico,nc_c), a(nr_bx, nc_c, ntico)
      double precision tmp
C      write(0,*) "bx_x, bx_i, bx_p, nr_bx, ntico, nc_bx, c,
C     * ldc, nc_c, a=", bx_x, bx_i, bx_p, nr_bx, ntico, nc_bx, c,
C     * ldc, nc_c, a
      do 3 icc=1,nc_c
         do 2 icb=1,nc_bx
            do 1 ii=bx_p(icb)+1,bx_p(icb+1)
               irb=bx_i(ii)
               it=irb/nr_bx
               irb1=mod(irb, nr_bx)+1
               tmp=c(icb+it*ldc, icc)
               it=it+1
               a(irb1, icc, it)=a(irb1, icc, it)+bx_x(ii)*tmp
    1       continue
    2    continue
    3 continue
      return
      end subroutine mult_bxt

      subroutine solve_lut(lua, nr_a, nlua, pivot, b, nc_b, ntico, ilua)
C call lapack dgters() for solving a series of linea systems a_i%*%(b_{i-1}+b_i)=b_i
C The result is stored inplace in b.
C The sizes:
C  - LU matrices of a : (nr_a, nr_a, nlua)
C  - pivot : (nr_a, nlua)
C  - b : (nr_a, nc_c, ntico)
C  - ilua : (ntico). It is a index vector. ilua[i] indicate which lua corresponds to the i-th time point
      IMPLICIT NONE
      integer nr_a, nlua, pivot(nr_a,nlua), nc_b, ntico, ilua(ntico)
      integer it, i, j, INFO
      double precision lua(nr_a, nr_a, nlua), b(nr_a,nc_b,ntico)
      
C     run in time
      do 1 it=1,ntico
         if (it > 1) then
C           add previous term to b_it
            do 2 j=1,nc_b
               do 3 i=1,nr_a
                  b(i,j,it)=b(i,j,it)+b(i,j,it-1)
    3          continue
    2       continue
         endif
C        solve the it-th system
         call DGETRS('N', nr_a, nc_b, lua(1,1,ilua(it)), nr_a,
     *      pivot(1,ilua(it)), b(1,1,it), nr_a, INFO)
    1 continue
      return
      end subroutine solve_lut
