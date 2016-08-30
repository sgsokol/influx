/*
Calculate inplace a_t-=b_t%*%c_t where a, b, c are matrices,
t is index in time and %*% is a dot product.
b is a sparse matrix of size (nr_b*ntico, nc_b) given by its fields
_v, _i, and _j describing triplet storage.
the size of c which is dense array is (ldc, nc_c, ntico) and those
of a (also a 3D array) is (nr_b, nc_c, ntico).
The parameter ldc must be >= ncol(b)
The result is subtracted from a, so if a pure multiplication result is needed,
user must initialized a to 0 before call

To compile do in R
Sys.setenv(PKG_LIBS=file.path(system.file(package="rmumps"), "libs", paste("rmumps", .Platform$dynlib.ext, sep="")))
sourceCpp("mult_bxxc.cpp")

Author: Serguei Sokol, INRA, Toulouse, FRANCE
Copyright 2016, INRA
v 0.1 2016-02-24
v 0.2 2016-03-04 added solve_ieu(), margins of c are permuted (ntico is last now)
*/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace std;

// [[Rcpp::depends(rmumps)]]
#include <rmumps.h>

using namespace Rcpp;
using namespace std::placeholders;
int my_mod(int a, int b) { return a%b; }
int my_div(int a, int b) { return a/b; }

// [[Rcpp::export]]
void mult_bxxc(NumericVector a, List b, NumericVector c) {
   if (!b.inherits("simple_triplet_matrix")) {
      print(b.attr("class"));
      stop("Parameter b must be a simple_triplet_matrix");
   }
   IntegerVector bi=as<IntegerVector>(b["i"])-1; // 0-based
   IntegerVector bj=as<IntegerVector>(b["j"])-1; // 0-based
   NumericVector bv=as<NumericVector>(b["v"]);
   IntegerVector dia=a.attr("dim");
   if (dia.size() != 3) {
      stop("Parameter a must be a 3D array");
   }
   IntegerVector dic=c.attr("dim");
   int nr_b=dia[0], nc_c=dia[1], ntico=dia[2], nc_b=b["ncol"], ldc=dic[0];
   int iic, iia;
   if (dic[2] != ntico)
      stop("dim(c)[3] must be equal to dim(a)[3]");
   if (dic[1] != nc_c)
      stop("dim(c)[2] must be equal to dim(a)[2]");
   if (ldc < nc_b)
      stop("dim(c)[1] must be greater or equal to ncol(b)");
   
   IntegerVector it(bi.size()), irb(bi.size());
   std::transform(bi.begin(), bi.end(), irb.begin(), std::bind(my_mod, _1, nr_b));
   std::transform(bi.begin(), bi.end(), it.begin(), std::bind(my_div, _1, nr_b));
   
//Rcout << "nr_b=" << nr_b << "; ntico=" << ntico << std::endl;
//Rcout << "irb=" << std::endl;
//print(irb);
//Rcout << "it=" << std::endl;
//print(it);
   
   for (int icc=0; icc < nc_c; icc++) {
      for (int i=0; i < bv.size(); i++) {
         iic=bj[i]+(icc+it[i]*nc_c)*ldc;
         iia=irb[i]+(icc+it[i]*nc_c)*nr_b;
//Rcout << "irb=" << irb[i] << ";\t jcb=" << bj[i] << ";\tit=" << it[i] << ";\tiic=" << iic << ";\tiia=" << iia << ";\tc=" << c[iic] << std::endl;
         a[iia]-=bv[i]*c[iic];
      }
   }
   //return R_NilValue;
}

// [[Rcpp::export]]
void solve_ieu(vec& invdt, mat& x0, mat& M, ListOf<XPtr<Rmumps>> ali, cube s, ivec& ilua) {
   // solve an ode system by implicite euler scheme
   // The system is defgined as M*dx/dt=a*x+s where M is a diagonal
   // matrix given by its diagonal vector M (which has a form of matrix for
   // term-by-term multiplication with x0)
   // In discrete terms
   // (M/dt_i-a)*x_i=(M/dt_i)*x_(i-1)+s_i
   // The rmumps matrix (M/dt_i-a) is stored in list ali as XPtr<Rmumps>
   // invdt is 1/dt
   // x0 is the startin value at t0
   // The source term is in array s, its last margin is time
   // The ilua[i] gives the list item number in ali for a given dt_i.
   // Calculations are done in-place so s is modified and contains the
   // solution on exit. The others parameters are not modified.
   // s_i can be a matrix or a vector(== 1-column matrix)
   unsigned int nti=invdt.size();
   unsigned int nxrow=x0.n_rows, nxcol=x0.n_cols;
   // sanity control
   if (M.n_rows != nxrow)
      stop("nrow(M) != nrow(x0)");
   if (M.n_cols != nxcol)
      stop("nrow(M) != nrow(x0)");
   if (s.n_rows != nxrow)
      stop("dim(s)[1] != nrow(x0)");
   if (s.n_cols != nxcol)
      stop("dim(s)[2] != ncol(x0)");
   if (s.n_slices != nti)
      stop("dim(s)[3] != length(invdt)");
   if (ilua.size() != nti)
      stop("length(ilua) != length(invdt)");
   // prepare starting s
   for (unsigned int i=0; i < nti; i++) {
//Rcout << "i=" << i << std::endl;
      if (i == 0)
         s.slice(i)=s.slice(i)+(M%x0)*invdt[i];
      else
         s.slice(i)=s.slice(i)+(M%s.slice(i-1))*invdt[i];
//Rcout << "s1 s_i=" << s.begin_slice(i) << std::endl;
//Rcout << "s2, ilua_i=" << ilua[i] << std::endl;
      //RObject e(ali(ilua[i]-1));
//print(as<Environment>(e.attr(".xData"))["copy"]);
//Rcout << "s3" << std::endl;
      //XPtr<Rmumps> ptr(as<XPtr<Rmumps>>(ali(ilua[i]-1)));
//Rcout << "s4" << std::endl;
      ali[ilua[i]-1]->solveptr(s.begin_slice(i), nxrow, nxcol);
//Rcout << "s5" << std::endl;
   }
}

#include <stdint.h>
typedef struct ij32 ij32;
struct ij32 {
   int i;
   int j;
};

typedef union ui64 ui64;
union ui64 {
   ij32     ij;
   uint64_t b;
};

typedef std::pair<size_t, ui64> iui64;

// [[Rcpp::export]]
IntegerVector match_ij(IntegerVector ix, IntegerVector jx, IntegerVector ti, IntegerVector tj) {
   // match ix,jx-couple in ti,tj-table and return their 1-based positions (0 for non matched couples)
   size_t ni=ix.size(), nti=ti.size();
   
   // Method: put i and j in 64 unsigned vectors, sort and match that vectors
   // init ij and tij
   std::vector<iui64> ij(ni), tij(nti);
   for (size_t i=0; i < ni; i++) {
      ij[i].second.ij.i=ix[i];
      ij[i].second.ij.j=jx[i];
      ij[i].first=i;
   }
   for (size_t i=0; i < nti; i++) {
      tij[i].second.ij.i=ti[i];
      tij[i].second.ij.j=tj[i];
      tij[i].first=i;
   }
   // sort and order them
   std::sort(ij.begin(), ij.end(), [](iui64 a, iui64 b) {return a.second.b < b.second.b;});
   std::sort(tij.begin(), tij.end(), [](iui64 a, iui64 b) {return a.second.b < b.second.b;});
   // run through to get matches
   IntegerVector m(ni, 0);
   for (size_t i=0, it=0; i < ni && it < nti;) {
      if (ij[i].second.b == tij[it].second.b) {
         m[ij[i].first]=tij[it].first+1;
         i++;
         it++;
         continue;
      } else if (ij[i].second.b < tij[it].second.b) {
         i++;
         continue;
      } else {
         it++;
      }
   }
   return(m);
}
/* no timing enhancement
// [[Rcpp::export]]
NumericVector crossprod_st(List x, NumericVector y_) {
   // dot product of simple triplet matrix x (m x n) and a dense matrix y (n x k)
   //Rcout << "here" << endl;
   int m=x["nrow"], n=x["ncol"], k;
   //print(x);
   Dimension ydim(2);
   //Rcout << "here 2" << endl;
   if (y_.hasAttribute("dim")) {
      ydim=y_.attr("dim");
   } else {
      ydim=Dimension(y_.size(), 1);
   }
   k=ydim[1];
   if (n != ydim[0])
      stop("ncol(x) != nrow(y)");
   IntegerVector ix_=as<IntegerVector>(x["i"])-1, jx_=as<IntegerVector>(x["j"])-1;
   NumericVector vx_=as<NumericVector>(x["v"]);
   int *ix=ix_.begin(), *jx=jx_.begin();
   double *vx=vx_.begin(), *y=y_.begin();
   NumericVector r_(m*k, 0.);
   r_.attr("dim") = Dimension(m, k);
   double *r=r_.begin();
   //print(ix_);
   //print(jx_);
   //print(vx_);
   for (int jy=0; jy < k; jy++, y+=n, r+=m) {
      for (int iv=0; iv < vx_.size(); iv++) {
         r[ix[iv]]+=vx[iv]*y[jx[iv]];
      }
   }
   return r_;
}
*/
