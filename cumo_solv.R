# R code for reading the linear system from files and
# solving it.

# 2008-04-16 sokol: started all cumomer solving (one iteration for future min pb)
# 2008-04-21 sokol: started trisparse_solv()

source("opt_cumo_tools.R");
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
bfl=Af[,n+1];
qrA=qr(A);

# make sure that free params choice leads to not singular matrix
if (qrA$rank != n) stop("Bad choice of free fluxes. The matrix is singular.");

ifree=(bfl==-1);
#print(ifree);
# try to solve fluxes with random positive free fluxes (uniform on [0;1])
nfree=sum(ifree);
bfl[ifree]=runif(nfree);
#print(b);

fl=solve(qrA, bfl);
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
