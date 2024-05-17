suppressPackageStartupMessages(library(nlsic)); # for non negative least square
#suppressPackageStartupMessages(library(Matrix, warn=F, verbose=F)); # for sparse matrices
#options(Matrix.quiet=TRUE)
suppressPackageStartupMessages(library(slam)); # for quick sparse matrices
suppressPackageStartupMessages(library(parallel))
#use_magma=suppressWarnings(suppressPackageStartupMessages(require(magma, quietly=T)))
#use_magma=F
##suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(RcppArmadillo))
suppressPackageStartupMessages(library(bitops))
suppressPackageStartupMessages(library(arrApply)) # for fast apply() on arrays
suppressMessages(suppressPackageStartupMessages(library(rmumps)))
suppressMessages(suppressPackageStartupMessages(library(multbxxc))) # auxiliary C++ routines
suppressPackageStartupMessages(library(kvh));
compiler::enableJIT(0)
#suppressPackageStartupMessages(library(Rdsm)); # for shared memory on cluster
# get some common tools
source(file.path(dirr, "tools_ssg.R"), echo=FALSE)
source(file.path(dirr, "psoptim_ic.R"), echo=FALSE)
if (case_i)
    source(file.path(dirr, "funlab.R"), echo=FALSE)
