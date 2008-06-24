# R code for reading the linear system from files and
# solving it.

# 2008-04-14 sokol: started
# 2008-04-16 sokol: started all cumomer solving (one iteration for future min pb)
# 2008-04-21 sokol: started trisparse_solv()

trisparse_solv=function(Al, Ac, Au, spA, b, method="dense") {
   # solve A*x=b where A=tridiag(Al,Ac,Au)+spA and b is dense
   # temporary solution by qr()
   n=length(Ac);
   if (method=="dense") {
      # fulfill a matrix
      A=matrix(0., n , n);
      for (i in 1:n) {
         if (i > 1) {
            A[i,i-1]=Al[i];
         }
         A[i,i]=Ac[i];
         if (i < n) {
            A[i,i+1]=Au[i];
         }
      }
      A=A+as.matrix(spA);
      x=solve(A,b);
   } else if (method=="sparse"){
      # sparse
      # fulfill a matrix
      require(Matrix);
      A=Matrix(0., n , n);
      for (i in 1:n) {
         if (i > 1) {
            A[i,i-1]=Al[i];
         }
         A[i,i]=Ac[i];
         if (i < n) {
            A[i,i+1]=Au[i];
         }
      }
      A=A+spA;
      x=solve(A,b);
#q=qr(A);
#print("qr");
#print(q@V);
#print(q@R);
#print(q@p);
   } else if (method=="smw") {
      # Sherman-Morrison-Woodbury for low rank matrix modification
      require(matrid, lib.loc="/home/sokol/R/lib");
      atri=new("matrid", Al, Ac, Au);
      # extract non zero columns from spA
      e=integer(0);
      m=0;
      if (class(spA)=="Matrix") {
         for (j in 1:n) {
            if (spA@p[j+1] > spA@p[j]) {
               e=c(e,j);
               m=m+1;
            }
         }
         s=as.matrix(spA[,e]);
      } else {
         # spA is treated as a dense matrix
         for (j in 1:n) {
            if (any(spA[,j])) {
               e=c(e,j);
               m=m+1;
            }
         }
         s=as.matrix(spA[,e]);
      }
      A=new("matridm", atri, s, e);
      x=qr.solve(A,b);
   } else {
      stop(paste("Unknown method '", method, "'", sep=""));
   }
#print("A");
#print(A);
}
# get cmd line options
opts=commandArgs();
ipref=which(opts=="--pref");
stopifnot(ipref > 0);
pref=opts[ipref+1];

# profile or not profile?
prof=(length(which(opts=="--prof")) > 0);

#pref="example1";

# R profiling
if (prof) {
   Rprof(paste(pref, ".Rprof", sep=""));
}

# read 0-weight system
fsys=paste(pref, "_fl_matrix.txt", sep="");
Af=as.matrix(read.table(fsys, header=TRUE, row.names=1, sep="\t"));
n=NROW(Af);
if (NCOL(Af) != n+1) stop(paste("Wrong row-column number in file", fsys));
A=Af[,1:n];
#print(A);
b=Af[,n+1];
qrA=qr(A);

# make sure that free params choice leads to not singular matrix
if (qrA$rank != n) stop("Bad choice of free fluxes. The matrix is singular.");

ifree=(b==-1);
#print(ifree);
# try to solve fluxes with random positive free fluxes (uniform on [0;1])
nfree=sum(ifree);
b[ifree]=runif(nfree);
#print(b);

fl=solve(qrA, b);
#print(fl);

# make all-weight systems
expr=parse(paste(pref, "_sym.R", sep=""));
for (i in 1:10) {
   eval(expr);
}
# nw is initialized in _sym.R
#str(spA[[1]]) # to see slots of spA[[1]]

# end of R profiling
if (prof) {
   Rprof(NULL);
}
