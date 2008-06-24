# QR decomposition and solving of a tridiagonal linear system
# Ax=b
#
# 2008-06-03 v1.0 (first release)
# Author: serguei.sokol@insa-toulouse.fr

# example
# a=new("matrid",c(NA, 1, 2),c(2, 2, 2),c(2, 2, NA))
# ap=matrix(c(2,1,0, 2,2,2, 0,2,2), 3, 3)
# b=c(6, 11, 10);
# solution should be (1,2,3)
# qr.solve(a,b);
#
# modified laplace example:
# n=5
# a=new("matrid", rep(1,n), rep(-2,n), rep(1,n));
# s=matrix(rnorm(n),n,1);
# e=c(n);
# A=new("matridm", a, s, e);
# x.exact=1:n;
# b=multmm(A,x.exact)
# qr.solve(A,b);

setClass("matrid", representation(l="numeric", d="numeric", u="numeric",
      dd="numeric", uu="numeric", Dim="numeric", qrd="logical"));
#method.skeleton("qr", "matrid")
#method.skeleton("solve", "matrid")
setMethod("initialize",
          "matrid",
   function(.Object,
         l=numeric(0),
         d=numeric(0),
         u=numeric(0),
         dd=numeric(0),
         uu=numeric(0),
         qrd=logical(0)) {
      if (nargs() > 1) {
         if (length(l) != length(d) || length(d) != length(u))
            stop("specified l, d, u are of different lengths")
         if (mode(l) != "numeric")
            stop("specified l is not of numeric mode")
         if (mode(d) != "numeric")
            stop("specified d is not of numeric mode")
         if (mode(u) != "numeric")
            stop("specified u is not of numeric mode")
         .Object@l <- as.double(l);
         .Object@d <- as.double(d);
         .Object@u <- as.double(u);
         n=length(.Object@d);
         if (length(.Object@dd)) {
            .Object@dd <- as.double(dd);
            .Object@uu <- as.double(uu);
         } else {
            # take additional memory for R
            .Object@dd <- double(n);
            .Object@uu <- double(n);
         }
         if (length(qrd)==0) {
            .Object@qrd <- FALSE
         } else {
            .Object@qrd <- as.logical(qrd)
         }
         .Object@Dim=c(n,n);
      }
      .Object;
   })
   
# Class matridm extends tridiagonal matrices (a) to almost tridiagonal ones (A).
# m dense columns stored in s are added to the matrix: A=a+s*E^t.
# Here E^t is an operator defined by an integer vector e which
# contains numbers of m dense columns in the whole matrix. For example,
# e[i]=j means i-th column in s takes j-th column (both 1-based) in A
setClass("qr");
setClassUnion("matrix_qr", c("matrix", "qr"));
setClass("matridm",
      representation(
      s="matrix",
      e="integer",
      m="integer",
      amm="matrix_qr"),
      
      contains=c("matrid"));
      
#method.skeleton("qr", "matrid")
#method.skeleton("solve", "matrid")
setMethod("initialize",
          "matridm",
   function(.Object,
         a=new("matrid"),
         s=matrix(0),
         e=numeric(0)) {
      if (nargs() > 1) {
         if (mode(s) != "numeric")
            stop("specified s is not of numeric mode")
         if (mode(e) != "numeric")
            stop("specified s is not of numeric mode")
         if (NROW(s) != a@Dim[1])
            stop("specified a and s have different row number")
         if (NCOL(s) != length(e))
            stop("specified s and e have different column number")
         .Object@l <- a@l
         .Object@d <- a@d
         .Object@u <- a@u
         .Object@dd <- a@dd
         .Object@uu <- a@uu
         .Object@qrd <- a@qrd
         .Object@Dim <- a@Dim
         n<-.Object@Dim[1];
         .Object@s<-as.matrix(s,n);
         .Object@e<-as.integer(e);
         .Object@m<-NCOL(.Object@s);
         m<-.Object@m;
         .Object@amm<-matrix(0,m,m);
         if (a@qrd && m) {
            # compleat QR decomposition of low rank part
            # solve m systems with s as rhs
            res<-.Fortran("tridqrsolv",
               N=n,
               l=.Object@l,
               d=.Object@d,
               u=.Object@u,
               dd=.Object@dd,
               uu=.Object@uu,
               s=.Object@s,
               M=.Object@m,
               NAOK=TRUE,
               DUP=FALSE,
               PACKAGE="matrid");
#show(as.matrix(s,n,m));
            # qr decomposition of dense m by m auxiliary matrix
            .Object@amm=qr(diag(1.,.Object@m)+.Object@s[.Object@e,]);
         }
      }
      .Object;
   })

setMethod("qr",
          "matrid",
function(x,...) {
   if (x@qrd) {
      # already decomposed. Nothing to do
      return(x);
   }
   # input consistency
   if (length(x@l) == 0) stop ("x@l is of length 0");
   if (length(x@d) == 0) stop ("x@d is of length 0");
   if (length(x@u) == 0) stop ("x@u is of length 0");
   if (length(x@l) != length(x@d) || length(x@d) != length(x@u))
         stop("specified l, d, u are of different lengths");
   n=x@Dim[1];
   # fortran call
   res<-.Fortran("tridqr",
         N=n,
         l=x@l,
         d=x@d,
         u=x@u,
         dd=x@dd,
         uu=x@uu,
         qrd=x@qrd,
         NAOK=TRUE,
         DUP=FALSE,
         PACKAGE="matrid");
   if (class(x) == "matridm" && x@m) {
      # solve m systems with s as rhs
      res<-.Fortran("tridqrsolv",
         N=n,
         l=x@l,
         d=x@d,
         u=x@u,
         dd=x@dd,
         uu=x@uu,
         s=x@s,
         M=x@m,
         NAOK=TRUE,
         DUP=FALSE,
         PACKAGE="matrid");
      # qr decomposition of dense m by m auxiliary matrix
#cat("DEBUG before amm update\n");
#show(x);
      x@amm=qr(diag(1.,x@m)+x@s[x@e,]);
   }
#cat("DEBUG\n");
#show(x);
   return(x);
})
setMethod("qr.solve",
          "matrid",
function(a,b,tol=1e-07) {
   if (!a@qrd) {
      # not yet decomposed. do it
      a=qr(a);
   }
   # input consistency
   if (length(a@l) == 0) stop ("a@l is of length 0");
   if (length(a@d) == 0) stop ("a@d is of length 0");
   if (length(a@u) == 0) stop ("a@u is of length 0");
   if (length(a@dd) == 0) stop ("a@dd is of length 0");
   if (length(a@uu) == 0) stop ("a@uu is of length 0");
   if (length(a@l) != length(a@d) || length(a@d) != length(a@u) ||
         length(a@d) != length(a@dd) || length(a@d) != length(a@uu))
      stop("specified l, d, u, dd, uu and dd are of different lengths");
   n=length(a@d);
   if (n != NROW(b)) stop ("NROW(b) is diffent from length(a@d)");
   if (is.matrix(b)) {
      m=NCOL(b);
   } else {
      m=1;
   }
   cb=class(b);
   # make local copy of b
   b=as.double(b);
   # fortran call
#tridqrsolv(N,l,d,u,dd,uu,b,M)
   res<-.Fortran("tridqrsolv",
         N=n,
         l=a@l,
         d=a@d,
         u=a@u,
         dd=a@dd,
         uu=a@uu,
         b=b,
         M=m,
         NAOK=TRUE,
         DUP=FALSE,
         PACKAGE="matrid");
#cat("DEBUG\n")
#show(b);
#show(a);
   if (class(a)=="matridm") {
      # solving slightly modified system a*x=(at+s*E^t)x=b
      # by Sherman-Morrison-Woodbury formula.
      # where at is an easly qr invertable n*n matrix (tridiagonal in this class),
      # s is n*m matrix with m << n (very low rank matrix)
      # e is an m-length integer vector with column numbers that s columns
      # take in a. So if e[i]=k than i-th column in s is a k-th column
      # in a. Naturally, it must be 1<=i<=m, 1<=k<=n.
      # (E^t in problem formulation is an operator based on e vector.)
      # b is a rhs n-length vector.
      # For repeated system solving with various rhs,
      # call a=qr(a) before tridm.solve() repeated calls.
      #print("here matridm");
      bm=as.matrix(b,n)[a@e,];
      bm=qr.solve(a@amm,bm);
#show(bm);
      #diag(1.,a@m)+v[e,],b[e]);
      b=b-a@s%*%bm;
   }
   if (cb=="matrix") {
      return(as.matrix(b,n,m));
   } else {
      return(b);
   }
})
setMethod("show",
          "matrid",
function(object) {
   # show tridiagonal matrix just like three column matrix
   # first l and last u values must be ignored
   n=object@Dim[1];
   alm="";
   suf="";
   extended=(class(object)=="matridm");
   if (extended) {
      alm="almost ";
      suf="m";
   }
   cat(paste(alm, "tridiagonal matrix of class matrid",suf,"\n", sep=""));
   if (object@qrd) {
      cat("QR decomposed\n");
      cat("Q Householder vectors in diagonal strorage:\n");
      print(cbind(l=object@l, d=object@d));
      cat("\nR in diagonal strorage:\n");
      print(cbind(dd=object@dd, u=object@u, uu=object@uu));
      if (extended) {
         cat("\ndense auxiliary matrix:\n");
         print(object@amm);
      }
   } else {
      print(cbind(l=object@l, d=object@d, u=object@u));
   }
   if (extended) {
      cat("added column number(s) e:\n");
      print(object@e);
      cat("added column(s) s:\n");
      print(object@s);
   }
})
setMethod("dim",
          "matrid",
function(x) {
   # retrun dim slot (vector c(n,n) where n is the dimension of the matrix)
   x@Dim;
})
setMethod("NROW",
          "matrid",
function(x) {
   # retrun dim slot (vector c(n,n) where n is the dimension of the matrix)
   x@Dim[1];
})
setMethod("NCOL",
          "matrid",
function(x) {
   # retrun dim slot (vector c(n,n) where n is the dimension of the matrix)
   x@Dim[2];
})

multmm<-function(a,b) {
   # multiply (almost) tridiagonal matrix a by a dense matrix b
   if (a@qrd) stop ("Cannot multiply QR decomposed matrix");
   n=a@Dim[1];
   b=as.matrix(b, n);
   m=NCOL(b);
   res=matrix(0,n,m);
   i=1:(n-1);
   for (j in 1:m) {
      res[,j]=a@d*b[,j]+c(0,a@l[i+1]*b[i,j])+c(a@u[i]*b[i+1,j],0);
   }
   if (class(a)=="matridm" && a@m) {
      res=res+a@s%*%b[e,];
   }
   res;
}

# for debugging purposes
qrh<-function(a) {
   # householder qr
   n=NROW(a);
   q=matrix(0, n, n);
   r=matrix(a, n, n);
   for (j in 1:n) {
      ac=r[j:n,j];
      alph=-sign(ac[1])*norm(ac);
      u=ac;
      u[1]=u[1]-alph;
      v=u/norm(u);
      q[j:n,j]=v;
      for (jj in j:n) {
         r[j:n,jj]=house(v,r[j:n,jj]);
      }
   }
   return(list(a=a,q=q,r=r));
}
qrhsolve<-function(l,b) {
   q=l$q;
   r=l$r;
   n=NROW(b);
   # b<-qt*b
   for (i in 1:n) {
      b=house(q[i:n,i],b,st=i);
#print(b);
   }
   # b<-R**(-1) * b
   b[n]=b[n]/r[i,i];
   for (i in seq(n-1,1,-1)) {
      b[i]=(b[i]-(b[(i+1):n]%*%r[i,(i+1):n]))/r[i,i];
   }
   b;
}
