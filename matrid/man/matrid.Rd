\name{matrid}
\alias{matrid}
\title{Decomposition and solving of tridiagonal matrix by QR Householder method}
\description{Solves the following linear problem:
 \deqn{A*x=b}
 using tridqr and tridqrsolv FORTRAN subroutines
}
\usage{qr(A)
qr.solve(A,b)}


 \cr Q Householder vectors are stored in l and
 d slots. The diagonals of R are strored in dd,u,uu slots.

\arguments{
 \item{A }{tridiagonal matrix of class matrid containing the coefficients of the system
\eqn{Ax=b}; 
 A has slots l,d,u which are real vectors of lenght N.
 N is size of the system, l, d and u are lower- main- and upper-
 diagonals of the original matrix A.
 \cr The result of QR decomposition is stored in l,d,u and two additional
 vector slots dd and uu.
 \cr A slot qrd is a logical indicating if A is QR
 decomosed (TRUE) or original (FALSE).
 \cr Householder vectors of Q are stored in l and
 d. The diagonals of R are strored in dd,u,uu vectors.
 If the d member of A have a names attribute, it will be used to
 label the output x.
 \cr If qr.solve() is called with
 original (not yet decomposed) matrix, qr() is automatically called
 before solving the system.
 \cr The solution x is stored in b which can be a matrix of multiple
 right hand sides. Each rhs (and its respective solution) takes one column in b.
 }
 \item{b}{one or multiple right hand sides of linear system. After calling qr.solve(),
 it is overwritten with corresponding x values.}
}

\value{
  The QR decomposition of input matrix for qr() and solution x of A*x=b
  for qr.solve().
}
\author{Serguei Sokol <serguei.sokol@insa-toulouse.fr>}
\examples{# simple 3x3 matrix
library(matrid);
A=new("matrid",
c(NA, 1, 2), # lower diagonal (first value is out of matrix => ignored)
c(2, 2, 2), # main diagonal
c(2, 2, NA)) # upper diagonal (last value is out of matrix => ignored)

# show the matrix diagonals
A

# right hand side
b=c(6, 11, 10);

# solve the system
qr.solve(A,b);
b # should be c(1,2,3)
}
\references{Golub, Gene H. & Van Loan, Charles F. (1996),
Matrix Computations (3rd ed.), Johns Hopkins, ISBN 978-0-8018-5414-9.}
\keyword{linear algebra}
\keyword{tridiagonal}
\keyword{QR decomposition}
