// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
param2fl_x=function(param, cjac=TRUE, labargs) {
   # translate free params (fluxes+scales) to fluxes and cumomers
   # or emus
   # local variabl assignments form labargs
   nm_local=ls(labargs)
   for (item in nm_local) {
      assign(item, labargs[[item]])
   }
   
   nb_xi=length(xi)
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   nb_fgr=nb_f$nb_fgr
   nb_meas=nb_f$nb_meas
   nb_sc=nb_f$nb_sc
   fg=nb_f$mu*param[nm$poolf] # the same alphabetical order
   names(fg)=nm$fgr
   pool[nm$poolf]=param[nm$poolf] # inject variable pools to pool vector

   # calculate all fluxes from free fluxes
   lf=param2fl(param, labargs)
   # prepare measurement pooling operations
   pwe[ipwe]=pool[ip2ipwe]
   spwe=tapply(pwe, pool_factor, sum)
   spwe=1./as.numeric(spwe[nm$measmat])
   pwe=pwe*spwe
   # construct the system A*x=b from fluxes
   # and find x for every weight
   # if fj_rhs is not NULL, calculate jacobian x_f
   if (is.null(labargs$incu)) {
      labargs$incu=incu=c(1, xi, double(nbc_x[nb_w+1L]))
   }
   if (cjac) {
      # derivation of fwrv fluxes by free parameters: free fluxes+concentrations
      mdf_dffp=df_dffp(param, jx_f$lf$flnx, nb_f, nm)
      jx_f$df_dffp=mdf_dffp
      if (is.null(labargs$x_f)) {
         labargs$x_f=x_f=matrix(0., nrow=sum(nb_x), ncol=nb_ff+nb_fgr)
      }
      
      if (length(ijpwef)) {
         dpwe=-pwe*spwe
         # dpwe is shortened to non zero entries in dpw_dpf
         dpwe=ifelse(ipf_in_ppw[ijpwef[,1L]]==ijpwef[,2L], (spwe+dpwe)[ijpwef[,1L]], dpwe[ijpwef[,1L]])
      }
   } else {
      x_f=NULL
   }
   
   # simulate labeling weight by weight
   ba_x=0
   for (iw in seq_len(nb_w)) {
      nb_c=spAb[[iw]]$nb_c
      emuw=ifelse(emu, iw, 1L)
      if (nb_c == 0) {
         next
      }
      ixw=(nbc_x[iw]+1L):nbc_x[iw+1]
      incuw=(1L+nb_xi)+ixw
      if (emu) {
         lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], incu, nm$emu[ixw], emu=emu)
      } else {
         lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], incu, nm$rcumo[ixw], emu=emu)
      }
      if (use_mumps) {
         am=Rmumps$new(lAb$A, FALSE)
         xw=try(am$solve(lAb$b))
         if (inherits(xw, "try-error")) {
            # find 0 rows if any
            izc=apply(lAb$A, 1L, function(v)sum(abs(v))<=1.e-10)
            izf=names(which(abs(lf$fwrv)<1.e-7))
            if (sum(izc) && length(izf)) {
               mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
                  "Zero rows in cumomer matrix A at weight ", iw, ":\n",
                  paste(rownames(lAb$A)[izc], collapse="\n"), "\n",
                  "Zero fluxes are:\n",
                  paste(izf, collapse="\n"), "\n",
                  sep="")
            } else {
               mes="Cumomer matrix is singular.\n"
            }
   #browser()
            return(list(x=NULL, fA=lAb$A, err=1L, mes=mes))
         }
      } else {
         wa=options(warn=2) # to cath singular matrix as error and not just a warning
         lua=try(if (use_magma) magma::lu(magma(as.matrix(lAb$A))) else lu(as.matrix(lAb$A), errSing=T), silent=T)
         options(wa) # restore warning situation
         if (inherits(lua, "try-error")) {
            # find 0 rows if any
            izc=apply(lAb$A, 1L, function(v)sum(abs(v))<=1.e-10)
            izf=names(which(abs(lf$fwrv)<1.e-7))
            if (sum(izc) && length(izf)) {
               mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
                  "Zero rows in cumomer matrix A at weight ", iw, ":\n",
                  paste(rownames(lAb$A)[izc], collapse="\n"), "\n",
                  "Zero fluxes are:\n",
                  paste(izf, collapse="\n"), "\n",
                  sep="")
            } else {
               mes="Cumomer matrix is singular.\n"
            }
   #browser()
            return(list(x=NULL, fA=lAb$A, err=1L, mes=mes))
         }
         b=as.matrix(lAb$b); # may have several columns if emu is TRUE
         #solve the system A*x=b
         #lsolv=trisparse_solv(lAb$A, lAb$b, iw, lf$fwrv, method="sparse")
         xw=lusolve(lua, b)
      }
      if (emu) {
         xw=c(xw, 1.-rowSums(xw))
      }
      incu[incuw]=xw
      if (cjac) {
         # calculate jacobian x_f
         # first, calculate right hand side for jacobian calculation
         # j_rhsw, b_x from sparse matrices
         # bind cumomer vector
         j_b_x=fx2jr(jx_f$lf$fwrv, spAb[[iw]], nb_f, incu)
         j_rhsw=j_b_x$j_rhsw%*%mdf_dffp
         b_x=j_b_x$b_x
         if (iw > 1) {
            if (ba_x > 0) {
               j_rhsw=j_rhsw+b_x%*%x_f[1L:ba_x,,drop=F]
            }
         }
         j_rhsw=as.double(j_rhsw)
         dim(j_rhsw)=c(nb_c, emuw*(nb_ff+nb_fgr))
         if (use_mumps) {
            j_rhsw=am$solve(j_rhsw)
         } else {
            j_rhsw=lusolve(lua, j_rhsw)
         }
         if (emu) {
            dim(j_rhsw)=c(nb_c, iw, nb_ff+nb_fgr)
            x_f[ba_x+seq_len(iw*nb_c),]=j_rhsw
            # m+N component
            #x_f[ba_x+iw*nb_c+seq_len(nb_c),]= -apply(j_rhsw, c(1L,3L), sum)
            x_f[ba_x+iw*nb_c+seq_len(nb_c),]=-arrApply(j_rhsw, 2, "sum")
         } else {
            x_f[ba_x+seq_len(nb_c),]=j_rhsw
         }
      }
      ba_x=ba_x+nb_x[iw]
   }
   names(incu)=c("one", nm$inp, nm$x)
   x=tail(incu, -nb_xi-1L)
   jx_f$x=x
   
   # calculate unreduced and unscaled measurements
   if (length(x) == ncol(measmat)) {
      mx=measmat%*%x+memaone
   } else {
      mx=measmat%*%x[nm$rcumo_in_cumo]+memaone
   }
   if (length(ipooled) > 1L) {
      mv=meas2sum%*%(pwe*mx)
   } else {
      mv=mx
   }
   jx_f$usimcumom=as.numeric(mv)
   names(jx_f$usimcumom)=nm$meas
   if (cjac) {
      # unreduced residuals derivated by scale params
      if (is.null(labargs$dur_dsc)) {
         labargs$dur_dsc=dur_dsc=matrix(0., nrow=nb_meas, ncol=nb_sc)
      }
      if (is.null(labargs$dux_dp)) {
         labargs$dux_dp=dux_dp=matrix(0., nb_meas, nb_ff+nb_sc+nb_poolf)
      }
      # measurement vector before pool ponderation
      # scale part of jacobian
      if (nb_sc > 0) {
         is2m=nb_f$is2m
         dur_dsc[is2m]=mv[is2m[,1]]
      }
#browser()
      # free flux part of jacobian (and partially free pools if present in x_f)
      if (nb_ff+nb_fgr > 0) {
         mffg=measmat%*%x_f
         if (length(ipooled) > 1L) {
            mffg=meas2sum%*%(pwe*mffg)
         }
      } else {
         mffg=Matrix(0., nrow=nb_meas, ncol=0L)
      }
      # free pool part of jacobian
      mff=mffg
      mpf=matrix(0., nrow=nb_meas, ncol=nb_f$nb_poolf)
      if (nb_f$nb_poolf > 0L) {
         if (length(ijpwef) > 0L) {
            # derivation of pool weights
            dpw_dpf@x[]=as.double(mx[ijpwef[,1L]]*dpwe)
            mpf[]=as.matrix(meas2sum%*%dpw_dpf)
            # growth flux depending on free pools
            if (nb_fgr > 0L) {
               mpf=mpf+as.matrix(mffg[,nb_ff+seq_len(nb_fgr),drop=F])
               mff=mffg[,seq_len(nb_ff)]
            } else {
               mff=mffg
            }
         }
      }
      mff=as.matrix(mff)
      if (nb_sc > 0) {
         vsc=c(1., param)[ir2isc]
         mff=vsc*mff
         mpf=vsc*mpf
      }
      # store usefull information in global list jx_f
      dux_dp[, seq_len(nb_ff)]=mff
      dux_dp[, nb_ff+seq_len(nb_sc)]=as.matrix(dur_dsc)
      dux_dp[, nb_ff+nb_sc+seq_len(nb_fgr)]=mpf
      jx_f$param=param
      jx_f$x_f=x_f
      jx_f$dux_dp=dux_dp
   }
   return(list(x=x, lf=jx_f$lf))
}

SEXP fwrv2Abr_rcpp(
NumericVector fwrv_in_,
List spAbr,
NumericVector incu_in_,
CharacterVector nm_rcumo,
bool getb,
bool emu) {

   // auxiliary functions
   Environment base_env_r_=Environment::base_env();
   Function rep_r_=base_env_r_[\"rep\"];
   Function c_r_=base_env_r_[\"c\"];
   // External R function declarations
   Function Matrix_head_r_=Environment(\"package:Matrix\")[\"head\"];
   Function Matrix_Matrix_r_=Environment(\"package:Matrix\")[\"Matrix\"];
   Function Matrix_sparseMatrix_r_=Environment(\"package:Matrix\")[\"sparseMatrix\"];

   // Input variable declarations and conversion
   vec fwrv(fwrv_in_.begin(), fwrv_in_.size(), true);
   vec incu(incu_in_.begin(), incu_in_.size(), true);

   // Output and intermediate variable declarations
   sp_mat a;
   sp_mat A;
   sp_mat b;
   List dimnames(A);
   List dimnames(b);
   ivec i;
   imat ind_a;
   SEXP ind_b;
   double nb_c;
   double w;
   vec x;

   // Translated code starts here
   nb_c=as<double>(spAbr[\"nb_c\"]);
   w=as<double>(spAbr[\"w\"]);
   if (nb_c == 0) {
      a=as<sp_mat>(Matrix_Matrix_r_(0, nb_c, nb_c));
      b=as<sp_mat>(Matrix_Matrix_r_(0, nb_c, 1));
return wrap(List::create(_[\"A\"]=a, _[\"b\"]=b))   }
   ind_a=as<imat>(spAbr[\"ind_a\"]);
   x=fwrv.at((ind_a(span(), (\"indf\")-1))-1);
   i=which(ind_a(span(), (\"ir0\")-1) == ind_a(span(), (\"ic0\")-1));
   x.at((i)-1)=-x.at((i)-1);
   A=as<sp_mat>(Matrix_sparseMatrix_r_(_[\"i\"]=ind_a(span(), (\"ir0\")-1) + 1, _[\"j\"]=ind_a(span(), (\"ic0\")-1) + 1, _[\"x\"]=x, _[\"dims\"]=vec({nb_c, nb_c})));
   dimnames(A)=List::create(NumericVector(NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).begin(), NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).end()), NumericVector(NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).begin(), NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).end()));
   if (getb) {
      if (emu) {
         ind_b=[[(spAbr, \"ind_b_emu\");
         x=-fwrv.at((ind_b(span(), (\"indf\")-1))-1) % incu.at((ind_b(span(), (\"indx1\")-1))-1) % incu.at((ind_b(span(), (\"indx2\")-1))-1);
         b=as<sp_mat>(Matrix_sparseMatrix_r_(_[\"i\"]=ind_b(span(), (\"irow\")-1), _[\"j\"]=ind_b(span(), (\"iwe\")-1), _[\"x\"]=x, _[\"dims\"]=vec({nb_c, w})));
         dimnames(b)=List::create(NumericVector(NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).begin(), NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).end()), NumericVector(linspace<ivec>(1, w, abs(w-1)+1).begin(), linspace<ivec>(1, w, abs(w-1)+1).end()));
      } else {
         ind_b=[[(spAbr, \"ind_b\");
         x=-fwrv.at((ind_b(span(), 0))-1) % incu.at((ind_b(span(), 1))-1) % incu.at((ind_b(span(), 2))-1);
         b=as<sp_mat>(Matrix_sparseMatrix_r_(_[\"i\"]=ind_b(span(), (\"irow\")-1), _[\"j\"]=as<vec>(rep_r_(1, ind_b.n_rows)), _[\"x\"]=x, _[\"dims\"]=vec({nb_c, 1})));
         dimnames(b)=List::create(NumericVector(NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).begin(), NA(nm_rcumo.begin(), std::max(0, std::min((int) (nb_c >= 0 ? nb_c : (nm_rcumo).n_elem+(nb_c)), (int) (nm_rcumo).n_elem))).end()), 1);
      }
return wrap(List::create(_[\"A\"]=A, _[\"b\"]=b))   } else {
return wrap(List::create(_[\"A\"]=A, _[\"b\"]=NumericVector(R_NilValue.begin(), R_NilValue.end())))   }

   return List::create(
      _(\"a\")=a,
      _(\"A\")=A,
      _(\"b\")=b,
      _(\"dimnames(A)\")=dimnames(A),
      _(\"dimnames(b)\")=dimnames(b),
      _(\"i\")=IntegerVector(i.begin(), i.end()),
      _(\"ind_a\")=ind_a,
      _(\"ind_b\")=ind_b,
      _(\"nb_c\")=nb_c,
      _(\"w\")=w,
      _(\"x\")=NumericVector(x.begin(), x.end())
   );

}
