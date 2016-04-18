// Calculate a_t+=bx_t%*%c_t where a, bx, c are matrices, t is index in time
// and %*% is a dot product.
// bx is a sparse matrix of size (nr_bx*ntico, nc_bx) given by its fields
// _x, _i, and _p describing column wise storage. The length of _x and _i is nb_x (virtual parameter)
// the size of c which is dense matrix is (ldc*ntico, nc_c) and those
// of a (also dense) is (nr_bx, nc_c, ntico).
// The parameter ldc may be different from ncol(bx) (if it is, then ldc > ncol(bx)
// a is supposed to be initialized to zero
// The result is added to a so it must be initialized to 0 before call
// if a pure multiplication result is needed.
// 
// Author: Serguei Sokol, INRA, Toulouse, FRANCE
// Copyright 2016, INRA
// v 0.1 2016-02-24

#include <Rcpp.h>
using namespace Rcpp;

void mult_bxxc(NumericalVector a, List bx, (NumericalVector c, int ntico, dirx) {
