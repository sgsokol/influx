C************************************************************************
C ALMOST TRIDIAGONAL MATRIX QR SOLVERS
C************************************************************************
 

      SUBROUTINE tridqr(N,l,d,u,dd,uu,qrd)

      IMPLICIT NONE

C-------------------------------------------------------------------------
C QR decomposition of a tridiagonal system of linear equations A*X=b
C where the non-zero diagonal elements are stored in l (below diagonal), 
C d (on the diagonal), u (above diagonal)
C Housholder vectors for Q matrix calculations are stored in l and d (in
C diagonal storage) and R matrix is stored in dd, u, uu
C Once the matrix is decomposed, the system A*x=b can be solved by
C calling tridqrsolv(...).
C-------------------------------------------------------------------------

      INTEGER            N
      DOUBLE PRECISION   l(N)
      DOUBLE PRECISION   d(N)
      DOUBLE PRECISION   u(N)
      DOUBLE PRECISION   dd(N)
      DOUBLE PRECISION   uu(N)
      INTEGER            qrd
C
      INTEGER            I
      DOUBLE PRECISION   temp
C      write(0,*) N, qrd
C      write(0,*) l
C      write(0,*) d
C      write(0,*) u
C-------------------------------------------------------------------------
      DO 10 I=1,N-1
C
C        Q vector
C        alpha
         dd(I)=-sign(sqrt(d(I)*d(I)+l(I+1)*l(I+1)),d(I))
         d(I)=d(I)-dd(I)
         temp=1./sqrt(d(I)*d(I)+l(I+1)*l(I+1))
         d(I)=d(I)*temp
         l(I+1)=l(I+1)*temp
C
C        R vectors
C        next
         temp=2*(d(I)*u(I)+l(I+1)*d(I+1))
         u(I)=u(I)-d(I)*temp
         d(I+1)=d(I+1)-l(I+1)*temp
C        after next
         IF (I .LT. N-1) THEN
            temp=2*l(i+1)*u(i+1)
            uu(I)=-d(I)*temp
            u(I+1)=u(I+1)-l(I+1)*temp
         ENDIF
10    CONTINUE
C     Last column
      dd(N)=-d(N)
      d(N)=-1.
      qrd=1;
      
C      write(0,*) d
C      write(0,*) dd
      RETURN
      END SUBROUTINE tridqr


      SUBROUTINE tridqrsolv(N,l,d,u,dd,uu,b,M)

      IMPLICIT NONE

C-------------------------------------------------------------------------
C Solving a tridiagonal system of linear equations Q*R*x=b
C b can be multiple vectors (in usual for fortran column storage)
C Their number is M.
C x is stored in b so b is lost. If you need it, copy it before calling.
C Housholder vectors for Q matrix calculations are stored in l and d (in
C diagonal storage) and R matrix is stored in dd, u, uu
C To decompose a tridiagonal matrix A in Q*R, call tridqr(...).
C-------------------------------------------------------------------------

      INTEGER            N
      INTEGER            M
      DOUBLE PRECISION   l(N)
      DOUBLE PRECISION   d(N)
      DOUBLE PRECISION   u(N)
      DOUBLE PRECISION   dd(N)
      DOUBLE PRECISION   uu(N)
      DOUBLE PRECISION   b(N,M)
C
      INTEGER            I
      INTEGER            J
      DOUBLE PRECISION   temp
C-------------------------------------------------------------------------
      if (N .EQ. 0 .OR. M .EQ. 0) THEN
         write(*,*) "tridqrsolv: N or M is zero", N, M
         RETURN
      ENDIF
C      write(*,*) "rhs:", N, M
C      write(*,*) b
      DO 10 J=1,M
         DO 110 I=1,N-1
C
C           b=Qt*b
            temp=2*(d(I)*b(I,J)+l(I+1)*b(I+1,J))
            b(I,J)=b(I,J)-d(I)*temp
            b(I+1,J)=b(I+1,J)-l(I+1)*temp
110       CONTINUE
C        Last element
         b(N,J)=-b(N,J)
C         write(*,*) "forward pass:", J
C         write(*,*) b
C
C        R**(-1)*b
C        Last unknown
         b(N,J)=b(N,J)/dd(N)
C        Before last
         if (N .GT. 1) THEN
            b(N-1,J)=(b(N-1,J) - b(N,J)*u(N-1))/dd(N-1)
         ENDIF
C        the rest of unknown vector
         DO 120 I=N-2,1,-1
            b(I,J)=(b(I,J) - b(I+1,J)*u(I) - b(I+2,J)*uu(I))/dd(I)
120      CONTINUE
10    CONTINUE
C      write(*,*) "solution:"
C      write(*,*) b
      RETURN
      END SUBROUTINE tridqrsolv


      SUBROUTINE dense2trid(n,A,l,d,u)

      IMPLICIT NONE

C-------------------------------------------------------------------------
C Extract three main diagonals from a dense matrix A and put them
C in l (lower) d (diagonal) and u (upper) vectors, all of size n.
C The first value in l and the last value in u are not modified
C by this routine.
C Extracted diagonals in A are zeroed
C-------------------------------------------------------------------------

      INTEGER            n
      DOUBLE PRECISION   A(n,n)
      DOUBLE PRECISION   l(n)
      DOUBLE PRECISION   d(n)
      DOUBLE PRECISION   u(n)
C
      INTEGER            i
      INTEGER            j
      
      if (n .LE. 0) then
         RETURN
      endif
      d(1)=a(1,1)
      d(n)=a(n,n)
      a(1,1)=0.
      a(n,n)=0.
      if (n .GT. 1) then
         l(2)=a(2,1)
         u(n-1)=a(n-1,n)
         a(2,1)=0.
         a(n-1,n)=0.
      endif
      DO 10 i=2,n-1
         l(i+1)=a(i+1,i)
         d(i)=a(i,i)
         u(i-1)=a(i-1,i)
         a(i+1,i)=0.
         a(i,i)=0.
         a(i-1,i)=0.
10    CONTINUE

      RETURN
      END SUBROUTINE dense2trid

