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
#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace std;
using namespace Rcpp;

#include <string>

// [[Rcpp::depends(rmumps)]]
#include <rmumps.h>

#include <Rinternals.h>
#include <R.h>

using namespace std;
using namespace std::placeholders;
int my_mod(int a, int b) { return a%b; }
int my_div(int a, int b) { return a/b; }

// not exported (auxiliary) function)
constexpr unsigned int s2i(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (s2i(str, h+1) * 33) ^ str[h];
}

#define NMES 512
char mes[NMES]={0};

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
   //#pragma omp parallel for ordered private(iic,iia) schedule(dynamic) num_threads(2) // little or not at all acceleration
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
void solve_ieu(vec& invdt, const SEXP& x0_, mat& M, ListOf<RObject> ali, cube s, ivec& ilua) {
   // solve an ode system by implicite euler scheme
   // The system is defined as M*dx/dt=a*x+s where M is a diagonal
   // matrix given by its diagonal vector M (which has a form of matrix for
   // term-by-term multiplication with x0)
   // In discrete terms
   // (M/dt_i-a)*x_i=(M/dt_i)*x_(i-1)+s_i
   // The rmumps matrix (M/dt_i-a) is stored in list ali as XPtr<Rmumps>
   // or a plain dense inverted matrix.
   // invdt is 1/dt
   // x0 is the starting value at t0 (may be NULL)
   // The source term is in array s, its last margin is time
   // The ilua[i] gives the list item number in ali for a given dt_i.
   // Calculations are done in-place so s is modified and contains the
   // solution on exit. The others parameters are not modified.
   // s_i can be a matrix or a vector(== 1-column matrix)
   bool addx0=!Rf_isNull(x0_);
//Rcout << "addx0=" << addx0 << std::endl;
   mat x0;
   unsigned int nxrow, nxcol;
   if (addx0) {
      x0=as<mat>(x0_);
      nxrow=x0.n_rows;
      nxcol=x0.n_cols;
   } else {
      nxrow=M.n_rows;
      nxcol=M.n_cols;
   }
   unsigned int nti=invdt.size(), iali;
   NumericVector tmp;
   
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
   // prepare pointers to Rmumps, matrix and a switch
   vector<Rmumps*> alip(ali.size());
   vector<mat> alim(ali.size());
   vector<bool> alismu(ali.size());
   
   for (int i=0; i<ali.size(); i++) {
      alismu[i]=ali[i].inherits("Rcpp_Rmumps");
      if (alismu[i]) {
         alip[i]=as<XPtr<Rmumps>>(as<Environment>(ali[i].slot(".xData"))[".pointer"]);
      } else {
         tmp=as<NumericVector>(ali[i]);
         alim[i]=mat(tmp.begin(), nxrow, nxrow, false);
      }
   }
   
   // prepare starting s
   for (unsigned int i=0; i < nti; i++) {
//Rcout << "i=" << i << std::endl;
      if (i == 0 && addx0)
         s.slice(i)=s.slice(i)+(M%x0)*invdt[i];
      else if (i > 0)
         s.slice(i)=s.slice(i)+(M%s.slice(i-1))*invdt[i];
//Rcout << "s1 s_i=" << s.begin_slice(i) << std::endl;
//Rcout << "s2, ilua_i=" << ilua[i] << std::endl;
      //RObject e(ali(ilua[i]-1));
//print(as<Environment>(e.attr(".xData"))["copy"]);
//Rcout << "s3" << std::endl;
      //XPtr<Rmumps> ptr(as<XPtr<Rmumps>>(ali(ilua[i]-1)));
//Rcout << "s4" << std::endl;
      iali=ilua[i]-1;
      if (alismu[iali])
         alip[iali]->solveptr(s.begin_slice(i), nxrow, nxcol);
      else
         s.slice(i)=alim[iali]*s.slice(i);
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
// [[Rcpp::export]]
void bop(NumericVector& dst, const IntegerVector& mv, const string& sop, NumericVector& src) {
   // bop=bloc operation in place
   // src array is added (if sop=="+=") to dst[...]
   // or any other manipulation is made according to sop parameter
   // Both arrays are supposed to be of type 'double'
   // The operation is done 'in place' without new memory allocation for dst
   // src is reshaped and possibly replicated to fit the designated block of dst.
   // mv can be:
   //  - a 1 or 3 component vector describing the block: 1-margin number of dst, 2-offset, 3-length
   //    if only the margin is present than offest is 0 and length is the total length of this margin
   //  - a matrix of indexes. Its column number must be equal to the length(dim(dst)))
   //    each row of this matrix is a multidimensional index in dst array.
   // sop is one off: "=" (copy src to dst[]), "+=", "-=", "*=", "/="
   
   // layaout arrays as cubes
   uvec did, dis, dimv; // array dimensions
   bool indmat=mv.hasAttribute("dim");
   
   if (dst.hasAttribute("dim")) {
      did=as<uvec>(dst.attr("dim"));
   } else {
      did=uvec(1).fill(dst.size());
   }
   if (indmat) {
      // proceed index matrix
      dimv=as<uvec>(mv.attr("dim"));
      if (dimv.size() != 2)
         stop("Parameter mv must be a matrix, i.e. length(dim(mv))==2. We got "+to_string(dimv.size()));
      if (dimv[1] != did.size())
         stop("Column number in index matrix mv ("+to_string(dimv[1])+") must be equal to length(dim(dst)) ("+to_string(did.size())+")");
      // prepare index vector
      umat  mvu=as<umat>(mv);
//print(wrap(iv));
      vec vdst=vec(dst.begin(), dst.size(), false);
      vec vsrc=vec(src.begin(), src.size(), false);
      did.tail(did.size()-1)=cumprod(did.head(did.size()-1));
      did(0)=1;
      mvu=mvu-1;
#define iv sum(mvu.each_row() % did.t(), 1)
      vdst(iv)=vsrc;
#undef iv
      return;
   }
   uvec dfi(3), dla(3), dlen; // destination first and last indexes by margin, and lengths
   uvec cdid(3), cdis(3); // cube dimensions
   uvec b=as<uvec>(mv);
   if (b.size() != 3 && b.size() != 1)
      stop("Block vector b must be of length 1 or 3 (not "+to_string(b.size())+")");
   unsigned int margin=b[0]-1, ioff, len;
   ioff=b.size()==3 ? b[1] : 0;
   len=b.size()==3 ? b[2] : did[margin];
   // if 0-length operation, leave now and do nothing
   if (len == 0 || prod(did) == 0)
      return;
   if (b[0] < 1 || b[0] > did.size()) {
      stop("Margin value ("+to_string(b[0])+") is invalid. Must be in [1, length(dim(dst)]");
   }
   if (ioff >= did[margin]) {
      stop("Block offset ("+to_string(ioff)+") is invalid. Must be in [0, dim(dst)[b[0]]-1]");
   }
   if (len > did[margin]) {
      stop("Block length ("+to_string(len)+") is invalid. Must be in [0, dim(dst)[b[0]]]");
   }
   // prepare cdst dimensions (cdid)
   if (margin > 0) {
      cdid[0]=prod(did.head(margin));
      cdid[1]=did[margin];
      cdid[2]=prod(did.tail(did.size()-b[0]));
      dfi={0, ioff, 0};
      dla={cdid[0]-1, ioff+len-1, cdid[2]-1};
      //Rcout << "dla=" << dla << endl;
      //Rcout << "dfi=" << dfi << endl;
   } else if (margin == 0) {
      cdid[0]=did[0];
      cdid[1]=prod(did.tail(did.size()-1));
      cdid[2]=1;
      dfi={ioff, 0, 0};
      dla={ioff+len-1, cdid[1]-1, cdid[2]-1};
   }
   dlen=dla-dfi+1;
   //Rcout << "dlen=" << dlen << endl;
   cube cdst(dst.begin(), cdid[0], cdid[1], cdid[2], false);
   if (src.size() == 1) {
      switch(s2i(sop.c_str())) {
      case(s2i("=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]).fill(src[0]);
         break;
      case(s2i("+=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) += src[0];
         break;
      case(s2i("-=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) -= src[0];
         break;
      case(s2i("*=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) *= src[0];
         break;
      case(s2i("/=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) /= src[0];
         break;
      default:
         stop("Unknown operation '"+sop+"'");
      }
      return;
   }
   // resize src if needed
   if (src.hasAttribute("dim")) {
      dis=as<uvec>(src.attr("dim"));
   } else {
      dis=uvec(1).fill(src.size());
   }
   mat msrc(src.begin(), src.size(), 1, false);
   unsigned int pdis=prod(dis), plen=prod(dlen);
   if ( pdis < plen) {
      // replicate msrc as many times as needed
      msrc=repmat(msrc, 1, plen/pdis);
   } else if (pdis > plen) {
      stop("Destination is not big enough ("+to_string(plen)+") to accept source ("+to_string(pdis)+")");
   }
   cube csrc(msrc.begin(), dlen[0], dlen[1], dlen(2), false);
      switch(s2i(sop.c_str())) {
      case(s2i("=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2])=csrc;
         break;
      case(s2i("+=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) += csrc;
         break;
      case(s2i("-=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) -= csrc;
         break;
      case(s2i("*=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) %= csrc;
         break;
      case(s2i("/=")):
         cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2]) /= csrc;
         break;
      default:
         stop("Unknown operation '"+sop+"'");
      }
   if (sop == "=")
      cdst.subcube(dfi[0], dfi[1], dfi[2], dla[0], dla[1], dla[2])=csrc;
}
// [[Rcpp::export]]
void redim(NumericVector& x, uvec& di) {
   // write new dimension vector while keeping the old memory
   uvec dix;
   dix=x.hasAttribute("dim") ? as<uvec>(x.attr("dim")) : uvec(1).fill(x.size());
   if (prod(dix) != prod(di))
      stop("Space in x ("+to_string(prod(dix))+") is not equal to new one ("+to_string(prod(di))+")");
   x.attr("dim")=di;
}
// [[Rcpp::export]]
void resize(SEXP& x_, uvec& di) {
   // write new dimension vector while keeping the old memory
   // new memory cannot be greater than the very first allocation
   RObject x=as<RObject>(x_);
   //Rcout << "len=" << LENGTH(x_) << endl;
   //Rcout << "tlen=" << TRUELENGTH(x_) << endl;
#if R_Version(3, 4, 0) <= R_VERSION
   if (!IS_GROWABLE(x_)) {
      SET_GROWABLE_BIT(x_);
      SET_TRUELENGTH(x_, XLENGTH(x_));
   }
#endif
   //Rcout << "n=" << XLENGTH(x_) << endl;
   unsigned int pdi=prod(di);
   if (TRUELENGTH(x_) >= pdi)
      SETLENGTH(x_, pdi);
   else
      stop("Space in x ("+to_string(TRUELENGTH(x_))+") is not sufficient for new one ("+to_string(pdi)+")");
   x.attr("dim")=di;
}
// [[Rcpp::export]]
List ij2ijv_i(IntegerVector& ir, IntegerVector& jc) {
   // transforms a couple of index vectors ir and jc (ij of a sparse matrix)
   // with possibly repeated values into a vector of unique indexes of non zero values.
   // The response can be then used for repeated creation of sparse
   // matrices with the same pattern by calling iv2v()
   // i and j are supposed to be sorted in increasing order, column-wise (i runs first)
   if (ir.size() != jc.size()) {
      snprintf(mes, NMES, "Sizes of ir (%d) and jc (%d) must be equal", (int) ir.size(), (int) jc.size());
      stop(mes);
   }
   size_t n=ir.size(), last=0;
   IntegerVector iv(n);
   uvec iu(n), ju(n);
   if (n == 0) {
      return List::create(_["i"]=iu, _["j"]=ju, _["iv"]=iv);
   }
   iv[0] = 0;
   iu[0]=ir[0];
   ju[0]=jc[0];
   for (size_t ii=1; ii < n; ii++) {
      last += ir[ii] != ir[ii-1] || jc[ii] != jc[ii-1];
      iv[ii] = last;
      iu[last] = ir[ii];
      ju[last] = jc[ii];
   }
   last++;
   iu.resize(last);
   ju.resize(last);
   return List::create(_["i"]=as<vector<double>>(wrap(iu)), _["j"]=as<vector<double>>(wrap(ju)), _["iv"]=iv);
}
// [[Rcpp::export]]
NumericVector iv2v(IntegerVector& iv, NumericVector& v) {
   // sum values in v according to possibly repeated indexes in iv
   if (iv.size() != v.size()) {
      snprintf(mes, NMES, "Sizes of iv (%d) and v (%d) must be equal", (int) iv.size(), (int) v.size());
      stop(mes);
   }
   NumericVector res(iv[iv.size()-1]+1);
   for (auto i=0; i < iv.size(); i++)
      res[iv[i]] += v[i];
   return res;
}
// [[Rcpp::export]]
NumericVector mm_xpf(List x, NumericVector y_, IntegerVector lsel) {
   // Dot product of simple triplet matrix x (m x n) (measurement matrix) and a dense array y (n x k x l).
   // Only slices of y_ from lsel vector are used.
   // Result array has dimensions (m x len(lsel) x k), i.e. it is permuted on the fly.
   //Rcout << "here" << endl;
   int m=x["nrow"], n=x["ncol"], k; //, l;
   //print(x);
   Dimension ydim(3);
   //Rcout << "here 2" << endl;
   if (y_.hasAttribute("dim")) {
      ydim=y_.attr("dim");
   } else {
      stop("parameter y_ must be a 3D array");
   }
   k=ydim[1];
   // l=ydim[2];
   if (n != ydim[0])
      stop("ncol(x) != dim(y)[1]");
   IntegerVector ix_=as<IntegerVector>(x["i"])-1, jx_=as<IntegerVector>(x["j"])-1;
   lsel=lsel-1;
   NumericVector vx_=as<NumericVector>(x["v"]);
   int *ix=ix_.begin(), *jx=jx_.begin(), ly;
   double *vx=vx_.begin(), *y=y_.begin();
   NumericVector r_(m*k*lsel.size(), 0.);
   r_.attr("dim") = Dimension(m, lsel.size(), k);
   double *r=r_.begin();
   //print(ix_);
   //print(jx_);
   //print(vx_);
   for (int jy=0; jy < k; jy++) {
      for (int ll=0; ll < lsel.size(); ll++, r+=m) {
         ly=lsel[ll];
         y=y_.begin()+n*(jy+k*ly);
         for (int iv=0; iv < vx_.size(); iv++) {
            r[ix[iv]]+=vx[iv]*y[jx[iv]];
         }
      }
   }
   return r_;
}
// [[Rcpp::export]]
void jrhs_ff(List jrhs, List ff, NumericVector xpfw) {
   // xpfw[...]-=jrhs%*%ff
   // dim(xpfw)=nb_row, nb_ff+nb_fgr, ntico
   // dim(jrhs%*%ff)=nb_row*ntico, nb_ff+nb_fgr
   IntegerVector ir=as<IntegerVector>(jrhs["i"])-1;
   IntegerVector jr=as<IntegerVector>(jrhs["j"])-1;
   NumericVector vr=as<NumericVector>(jrhs["v"]);
//Rcout << "*ir=" << ir.begin() << std::endl;
   IntegerVector iff=as<IntegerVector>(ff["i"])-1;
   IntegerVector jff=as<IntegerVector>(ff["j"])-1;
   NumericVector vff=as<NumericVector>(ff["v"]);
   int i, j, k;
   int i1, i3;
   Dimension dix=xpfw.attr("dim");
   int nr=dix[0], nf=dix[1]; //, nt=dix[2];
   int ivr=0; 
   for (int ivff=0; ivff < vff.size(); ivff++) {
      k=jff[ivff];
      if (ivff > 0 && jff[ivff-1] != k)
         ivr=0; // column in ff has changed => restart r
      for (; ivr < vr.size(); ivr++) {
         j=jr[ivr];
         if (j < iff[ivff])
            continue; // find right column in r
         else if (j > iff[ivff])
            break; // move to the next sparse term in vff
         i=ir[ivr];
         i1=i%nr;
         i3=i/nr;
         //i2=k;
         xpfw[i1+nr*(k+nf*i3)]-=vr[ivr]*vff[ivff];
      }
   }
}
