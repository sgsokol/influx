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
