suppressPackageStartupMessages(library(nnls)); # for non negative least square
#suppressPackageStartupMessages(library(Matrix, warn=F, verbose=F)); # for sparse matrices
#options(Matrix.quiet=TRUE)
suppressPackageStartupMessages(library(slam)); # for quick sparse matrices
suppressPackageStartupMessages(library(parallel))
#use_magma=suppressWarnings(suppressPackageStartupMessages(require(magma, quietly=T)))
#use_magma=F
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(RcppArmadillo))
suppressPackageStartupMessages(library(rmumps))
suppressPackageStartupMessages(library(arrApply)); # for fast apply() on arrays
suppressPackageStartupMessages(library(multbxxc)); # auxiliary C++ routines
#suppressPackageStartupMessages(library(compiler));
compiler::enableJIT(0)
#suppressPackageStartupMessages(library(Rdsm)); # for shared memory on cluster
# get some common tools
source(file.path(dirx, "tools_ssg.R"))
source(file.path(dirx, "nlsic.R"))
source(file.path(dirx, "kvh.R"))
#loadcmp(file.path(dirx, "tools_ssg.Rc"))
#loadcmp(file.path(dirx, "nlsic.Rc"))
#loadcmp(file.path(dirx, "kvh.Rc"))
