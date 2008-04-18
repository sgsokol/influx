# R code for reading the linear system from files and
# solving it.

# 2008-04-14 sokol: start
# 2008-04-16 sokol: all cumomer solving (one iteration for future min pb)

# read 0-weight system
fsys="fwd_rev_matrix.txt";
Af=as.matrix(read.table(fsys, header=TRUE, row.names=1, sep="\t"));
n=NROW(Af);
if (NCOL(Af) != n+1) stop(paste("Wrong row-column number in file", fsys));
A=Af[,1:n];
b=Af[,n+1];
qrA=qr(A);

# make sure that free params choice leads to not singular matrix
if (qrA$rank != n) stop("Bad choice of free fluxes. The matrix is singular.");

ifree=(b==-1);
# try to solve fluxes with random positive free fluxes (uniform on [0;1])
nfree=sum(ifree);
b[ifree]=runif(nfree);

fl=solve(qrA, b);

# read x-weight systems
nw=10;
for (iw in 1:1) {
   # read symbolic matrix
   fsys=paste("cumo_w", iw, "_system.txt", sep="");
   sAb=as.matrix(read.table(fsys, header=TRUE, row.names=1, sep="\t", colClasses="character"));
   # make tridiagonal and rest symbolic matrix
   striA=
   # evaluate symbolic expression of fluxes into numeric values
}
