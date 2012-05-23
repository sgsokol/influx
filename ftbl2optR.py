#!/usr/bin/env python

"""
Transform an ftbl to R code which will solve an optimization of flux analysis
problem min(S) over \Theta, where S=||Predicted-Observed||^2_\Sigma^2
and \Theta is a vector of free fluxes (net+xch) and scaling parameters.
Predicted vector is obtained from cumomer vector x (calculated from
free fluxes and divided in chunks according to the cumo weight) by
multiplying it by the measurement matrices and scale factor, boths coming
from ftbl file. Observed values vector xo is extracted from ftbl file.
it is composed of flux and cumomer measurements.
\Sigma^2, covariance diagonal matrices sigma[flux|mass|label|peak]
are orginated from ftbl

usage: ./ftbl2optR.py organism
where organism is the ftbl informative part of file name
(before .ftbl), e.g. organism.ftbl
after execution a file organism.R will be created.
If it already exists, it will be silently overwritten.
The system Afl*flnx=bfl is created from ftbl file.

Important python variables:
Collections:
   netan - (dict) ftbl structured content;
   tfallnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
   collection;
   measures - (dict) exp data;
   rAb - (list) reduced linear systems A*x_cumo=b by weight;
   scale - unique scale names;
   nrow - counts scale names;
   o_sc - ordered scale names;
   o_meas - ordered measurement types;
org - (str) prefix of .ftbl  file like "PPP"
File names (str):
   n_ftbl (descriptor f_ftbl);
   n_R (R code) (f);
   n_fort (fortran code) (ff);
Counts:
   nb_fln, nb_flx, nb_fl (dependent fluxes: net, xch, total), nb_ffn, nb_ffx (free fluxes)
Index translators:
   fwrv2i - flux names to index in R:fwrv;
   cumo2i - cumomer names to index in R:x;
   ir2isc - mapping measurement rows indexes on scale index isc[meas]=ir2isc[meas][ir]
Vector names:
   cumos (list) - names of R:x;
   o_mcumos - cumomers involved in measurements;

Important R variables:
Scalars:
   nb_w, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
   nb_ffn, nb_ffx, nb_ff (free fluxes),
   nb_fcn, nb_fcx, nb_fc (constrained fluxes),
   nb_ineq, nb_param, nb_fmn
Name vectors:
   nm_cumo, nm_fwrv, nm_fallnx, nm_fln, nm_flx, nm_fl, nm_par,
   nm_ffn, nm_ffx,
   nm_fcn, nm_fcx,
   nm_mcumo, nm_fmn
Numeric vectors:
   fwrv - all fluxes (fwd+rev);
   x - all cumomers (weight1+weight2+...);
   param - free flux net, free flux xch, scale label, scale mass, scale peak
   fcn, fcx, fc,
   bp - helps to construct the rhs of flux system
   xi -cumomer input vector
   fallnx - complete flux vector (constr+net+xch)
   bc - helps to construct fallnx
   li - inequality vector (mi%*%fallnx>=li)
   ir2isc - measur row to scale vector replicator
   ci - inequalities for param use (ui%*%param-ci>=0)
   measvec,
   measinvvar,
   imeas,
   fmn
Matrices:
   Afl, qrAfl, invAfl,
   p2bfl - helps to construct the rhs of flux system
   mf, md - help to construct fallnx
   mi - inequality matrix (ftbl content)
   ui - inequality matrix (ready for param use)
   measmat - measmat*(x[imeas];1)=vec of simulated not-yet-scaled measurements
Functions:
   param2fl_x - translate param to flux and cumomer vector (initial approximation)
   cumo_cost - cost function (khi2)
   cumo_grad - finite difference gradient
   cumo_gradj - implicit derivative gradient
"""

# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package
# 2010-10-16 sokol: fortran code is no more generated, R Matrix package is used for sparse matrices.

if __name__ == "__main__":
    import sys;
    import os;
    import stat;
    import time;
    import copy;
    import getopt;

    me=os.path.realpath(sys.argv[0])
    sys.path.append(os.path.dirname(me));
    me=os.path.basename(me);

    from tools_ssg import *;
    import C13_ftbl;

    try:
        import psyco;
        psyco.full();
        #psyco.log();
        #psyco.profile();
    except ImportError:
        #print 'Psyco not installed, the program will just run slower'
        pass;

    def usage():
        sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--DEBUG] network_name[.ftbl]\n");

    #<--skip in interactive session
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "fullsys", "DEBUG"]);
    except getopt.GetoptError, err:
        #pass;
        sys.stderr.write(str(err)+"\n");
        usage();
        #sys.exit(1);

    fullsys=False;
    DEBUG=False;
    for o,a in opts:
        if o in ("-h", "--help"):
            usage();
            sys.exit(0);
        elif o=="--cost":
            cost=True;
        elif o=="--fullsys":
            fullsys=True;
        elif o=="--DEBUG":
            DEBUG=True;
        else:
            #assert False, "unhandled option";
            # unknown options can come from shell
            # which passes all options both to python and R scripts
            # so just ignore unknown options
            pass;
    #aff("args", args);##
    if len(args) != 1:
        usage();
        exit(1);
    org=os.path.basename(args[0])
    dirorg=os.path.dirname(args[0]) or '.'

    # cut .ftbl if any
    if org[-5:]==".ftbl":
        org=org[:-5];
    fullorg=os.path.join(dirorg, org)

    #-->
    #DEBUG=True;
    import ftbl2code;
    ftbl2code.DEBUG=DEBUG;
    #org="ex3";
    #org="PPP_exact";
    #DEBUG=True;
    if DEBUG:
        import pdb;


    n_ftbl=fullorg+".ftbl";
    n_R=fullorg+".R";
    #n_fort=fullorg+".f";
    f_ftbl=open(n_ftbl, "r");
    try:
        os.chmod(n_R, stat.S_IWRITE);
    except:
        pass;
    f=open(n_R, "w");

    # parse ftbl
    ftbl=C13_ftbl.ftbl_parse(f_ftbl);
    f_ftbl.close();

    # analyse network
    # reload(C13_ftbl);

    netan=C13_ftbl.ftbl_netan(ftbl);

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, fullsys);

#    f.write(
    """
# output flux repartition
cat("Dependent fluxes:\\n");
if (nb_fln) {
   print(paste(nm_fln,"net",sep="_"));
}
if (nb_flx) {
   print(paste(nm_flx,"xch",sep="_"));
}
cat("Free fluxes:\\n");
if (nb_ffn) {
   print(paste(nm_ffn,"net",sep="_"));
}
if (nb_ffx) {
   print(paste(nm_ffx,"xch",sep="_"));
}
cat("Constrained fluxes:\\n");
if (nb_fcn) {
   print(paste(nm_fcn,"net",sep="_"));
}
if (nb_fcx) {
   print(paste(nm_fcx,"xch",sep="_"));
}
"""
    f.write("""
if (TIMEIT) {
   cat("preopt  : ", date(), "\\n", sep="");
}
#browser()

# check if initial approximation is feasible
ineq=ui%*%param-ci;
param_old=param
if (any(ineq <= -1.e-10)) {
   cat("The following ", sum(ineq<= -1.e-10), " ineqalities are not respected at starting point:\\n", sep="");
   print(ineq[ineq<= -1.e-10,1]);
   # put them inside
   param=put_inside(param, ui, ci)
   if (any(is.na(param))) {
      if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
         # fatal error occurd
         stop(paste("put_inside: ", attr(param, "mes"), collapse=""))
      }
   } else if (!is.null(attr(param, "err")) && attr(param, "err")==0){
      # non fatal problem
      warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
   }
}

# zero crossing strategy
# inequalities to keep sens of net flux on first call to opt_wwrapper()
# if active they are removed on the second call to opt_wrapper()
mi_zc=NULL
li_zc=NULL
if (nb_fn && zerocross) {
   # add lower limits on [df].net >= zc for positive net fluxes
   # and upper limits on [df].net <= -zc for negative net fluxes
   nm_izc=c()
   ipos=names(which(fallnx[grep("[df].n.", nm_fallnx)]>=0.))
   ineg=names(which(fallnx[grep("[df].n.", nm_fallnx)]<0.))
   mi_zc=matrix(0, nrow=length(ipos)+length(ineg), ncol=nb_fallnx);
   colnames(mi_zc)=nm_fallnx
   if (length(ipos)) {
      nm_izc=c(nm_izc, paste("zc ", ipos, ">=", zc, sep=""));
      mi_zc[(1:length(ipos)),ipos]=diag(1., length(ipos));
   }
   if (length(ineg)) {
      nm_izc=c(nm_izc, paste("zc ", ineg, "<=", -zc, sep=""));
      mi_zc[length(ipos)+(1:length(ineg)),ineg]=diag(-1., length(ineg));
   }
   rownames(mi_zc)=nm_izc
   li_zc=rep(zc, length(nm_izc));
   ui_zc=cbind(mi_zc%*%(md%*%invAfl%*%p2bfl+mf),
      matrix(0., nrow=nrow(mi_zc), ncol=nb_sc));
   ci_zc=li_zc-mi_zc%*%mic;
   # remove redundant/contradictory inequalities
   nb_zc=nrow(ui_zc);
   nb_i=nrow(ui);
   ired=c();
   if (nb_zc > 0) {
      for (i in 1:nb_zc) {
         nmqry=nm_izc[i];
         for (j in 1:nb_i) {
            if ((isTRUE(all.equal(ui[j,],ui_zc[i,])) ||
               isTRUE(all.equal(ui[j,],-ui_zc[i,]))) &&
               abs(abs(ci[j])-abs(ci_zc[i]))<=1.e-2) {
#browser()
               # redundancy
               cat("inequality '", nmqry, "' redundant or contradictory with '", nm_i[j], "' is removed.\\n", sep="");
               ired=c(ired, i);
               break;
            }
         }
      }
   }
   if (!is.null(ired)) {
      # remove all ired inequalities
      ui_zc=ui_zc[-ired,,drop=F];
      ci_zc=ci_zc[-ired];
      nm_izc=nm_izc[-ired];
      mi_zc=mi_zc[-ired,,drop=F];
   }
   if (nrow(ui_zc)) {
      # add zc inequalities
      ui=rbind(ui, ui_zc)
      ci=c(ci, ci_zc)
      nm_i=c(nm_i, nm_izc)
      mi=rbind(mi, mi_zc)
   }
#print(ui)
#print(ci)
}

# set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
# for corresponding measurements
vr=param2fl_x(param, cjac=FALSE, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, spAbr);
if (!is.null(vr$err) && vr$err) {
   stop(vr$mes);
}
simvec=(measmat%*%c(vr$x[irmeas],1.));
if (DEBUG) {
   cat("initial simvec:\\n");
   print(simvec);
}
if (nb_ff < length(param)) {
   ms=measvec*simvec*measinvvar;
   ss=simvec*simvec*measinvvar;
   for (i in (nb_ff+1):length(param)) {
      im=(ir2isc==(i+1));
      param[i]=sum(ms[im])/sum(ss[im]);
   }
}

# prepare flux index conversion
ifwrv=1:nb_fwrv;
names(ifwrv)=nm_fwrv;
ifl_in_fw=ifwrv[paste("fwd", substring(c(nm_fln, nm_flx), 4), sep="")];
iff_in_fw=ifwrv[paste("fwd", substring(c(nm_ffn, nm_ffx), 4), sep="")];
# index couples for jacobian df_dfl, df_dffd
cfw_fl=crv_fl=cbind(ifl_in_fw, 1:nb_fl);
cfw_ff=crv_ff=cbind(iff_in_fw, 1:nb_ff);
crv_fl[,1]=(nb_fwrv/2)+crv_fl[,1];
crv_ff[,1]=(nb_fwrv/2)+crv_ff[,1];

# see if there are any active inequalities at starting point
ineq=ui%*%param-ci;
if (any(abs(ineq)<=1.e-10)) {
   cat("The following ", sum(abs(ineq)<=1.e-10), " ineqalitie(s) are active at starting point:\\n",
      paste(names(ineq[abs(ineq)<=1.e-10,1]), collapse="\\n"), "\\n", sep="");
}
""");

    f.write("""
# formated output in kvh file
fkvh=file("%(fullorg)s_res.kvh", "w");
"""%{
    "fullorg": escape(fullorg, "\\"),
});
    # main part: call optimization
    f.write("""
cat("influx\\n", file=fkvh);
cat("\\tversion\\t", vernum, "\\n", file=fkvh, sep="");
# save options of command line
cat("\\tR command line\\n", file=fkvh);
obj2kvh(opts, "opts", fkvh, indent=2);

# resume system sizes
obj2kvh(nb_sys, "system sizes", fkvh);

# save initial flux and cumomer distribution
cat("starting point\\n", file=fkvh);
names(param)=nm_par;
obj2kvh(param, "starting free parameters", fkvh, indent=1);

x=vr$x;
names(x)=nm_rcumo;
obj2kvh(cumo2mass(x), "starting MID vector", fkvh, indent=1);
# replace :i by #x1's
nm_mask=unlist(lapply(strsplit(names(x), ":"), function(v) {
   n=length(v);
   metab=paste(v[-n], collapse=":");
   mask=int2bit(v[n], len=clen[metab])
   c(metab, mask);
}))
nm_mask=t(matrix(nm_mask, nrow=2));
o=order(nm_mask[,1], nm_mask[,2]);
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T);
nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#");
names(nm_mask)=names(x);
names(x)=nm_mask;
obj2kvh(x[o], "starting cumomer vector", fkvh, indent=1);

fwrv=vr$fwrv;
n=length(fwrv);
names(fwrv)=nm_fwrv;
obj2kvh(fwrv, "starting fwd-rev flux vector", fkvh, indent=1);

f=vr$fallnx;
n=length(f);
names(f)=nm_fallnx;
obj2kvh(f, "starting net-xch01 flux vector", fkvh, indent=1);

rres=cumo_resid(param, cjac=TRUE, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr);
names(rres$res)=c(nm_meas, nm_fmn);
obj2kvh(rres$res, "starting residuals (simulated-measured)/sd_exp", fkvh, indent=1);

rcost=cumo_cost(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr);
obj2kvh(rcost, "starting cost value", fkvh, indent=1);

obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1);
btmp=c(p2bfl%*%param[1:nb_f$nb_ff]+bp);
names(btmp)=dimnames(Afl)[[1]];
obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1);
names(measvec)=nm_meas;
obj2kvh(measvec, "measurement vector", fkvh, indent=1);

#cat("mass vector:\\n");
#print_mass(x);

names(param)=nm_par;
""");
    f.write("""
control_ftbl=list(%(ctrl_ftbl)s);
"""%{
    "ctrl_ftbl": join(", ", (k[8:]+"="+str(v) for (k,v) in netan["opt"].iteritems() if k.startswith("optctrl_"))),
})
    f.write("""
opt_wrapper=function(measvec, fmn, trace=1) {
   if (method == "BFGS") {
      control=list(maxit=500, trace=trace);
      control[names(control_ftbl)]=control_ftbl;
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAbr);
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=trace);
      control[names(control_ftbl)]=control_ftbl;
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-07,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAbr);
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=trace);
      control[names(control_ftbl)]=control_ftbl;
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-07,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, measinvvar, ir2isc,
         fmn, invfmnvar, ifmn, spAbr);
   } else if (method == "nlsic") {
      control=list(trace=trace, btfrac=0.25, btdesc=0.75, maxit=50, errx=1.e-5,
         ci=list(report=F));
      control[names(control_ftbl)]=control_ftbl;
      res=nlsic(param, cumo_resid, 
         ui, ci, control, e=NULL, eco=NULL, flsi=lsi_fun,
         nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi,
         irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar,
         spAbr);
      if (res$err || is.null(res$par)) {
         # store res in kvh
         obj2kvh(res, "optimization aborted here", fkvh);
         warning(res$mes);
         res$par=rep(NA, length(param))
         res$cost=NA
      } else if (nchar(res$mes)) {
         warning(res$mes);
      }
   } else if (method == "ipopt") {
      control=list(max_iter=500, print_level=trace*5);
      control[names(control_ftbl)]=control_ftbl;
      tui=c(t(ui))
      eval_g=function(x, nb_f=nb_f, nm=nm_list, nb_w=nb_rw, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, bp=bp, fc=fc, xi=xi,
         imeas=irmeas, measmat=measmat, measvec=measvec, measinvvar=measinvvar,
         ir2isc=ir2isc, fmn=fmn, invfmnvar=invfmnvar, ifmn=ifmn, spAb=spAbr) {
         return(ui%*%x)
      }
      eval_jac_g=function(x, nb_f=nb_f, nm=nm_list, nb_w=nb_rw, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, bp=bp, fc=fc, xi=xi,
         imeas=irmeas, measmat=measmat, measvec=measvec, measinvvar=measinvvar,
         ir2isc=ir2isc, fmn=fmn, invfmnvar=invfmnvar, ifmn=ifmn, spAb=spAbr) {
         return(tui)
      }
      ui_row_spars=rep.int(1, ncol(ui))
      res=ipoptr(param, cumo_cost, cumo_gradj,
         lb=NULL, ub=NULL,
         eval_g=eval_g,
         eval_jac_g=eval_jac_g,
         eval_jac_g_structure=lapply(1:nrow(ui), function(i)ui_row_spars),
         constraint_lb=ci,
         constraint_ub=rep(1.e19, length(ci)),
         eval_h=NULL,
         eval_h_structure=NULL,
         opts=control,
         ipoptr_environment=new.env(),
         nb_f=nb_f, nm=nm_list, nb_w=nb_rw, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, bp=bp, fc=fc, xi=xi,
         imeas=irmeas, measmat=measmat, measvec=measvec, measinvvar=measinvvar,
         ir2isc=ir2isc, fmn=fmn, invfmnvar=invfmnvar, ifmn=ifmn, spAb=spAbr);
      res$par=res$solution
      names(res$par)=nm_par
      if(res$status != 0) {
         # store res in kvh
         obj2kvh(res$par, "optimization aborted here", fkvh);
         warning(res$message)
      }
   } else {
      stop(paste("Unknown minimization method '", method, "'", sep=""));
   }
   return(res);
}
if (optimize) {
   # check if at starting position all fluxes can be resolved
#browser();
   qrj=qr(jx_f$dr_dff);
   d=diag(qrj$qr)
   qrj$rank=sum(abs(d)>abs(d[1])*1.e-14)
   if (qrj$rank) {
      nm_uns=nm_ff[qrj$pivot[-(1:qrj$rank)]]
   } else {
      nm_uns=nm_ff
   }
   if (qrj$rank < nb_ff && !(least_norm || method!="nlsic")) {
      # Too bad. The jacobian of free fluxes is not of full rank.
      dimnames(jx_f$dr_dff)[[2]]=c(nm_ffn, nm_ffx);
      write.matrix(formatC(jx_f$dr_dff, 15), file="dbg_dr_dff_singular.txt", sep="\t");
      stop(paste("Provided measurements (isotopomers and fluxes) are not
   sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:
   ", paste(nm_uns, sep=", ", collapse=", "),
         "\nJacobian dr_dff is dumped in dbg_dr_dff_singular.txt", sep=""));
   }
   if (TIMEIT) {
      cat("optim   : ", date(), "\\n", sep="");
   }
   # pass control to the chosen optimization method
   res=opt_wrapper(measvec, fmn);
   if (any(is.na(res$par))) {
      obj2kvh(res, "optimization process informations", fkvh);
      stop("Optimization failed")
   }
#browser()
   if (zerocross && !is.null(mi_zc)) {
      # inverse active "zc" inequalities
      nm_inv=names(which((ui%*%res$par-ci)[,1]<=1.e-10))
      i=grep("^zc ", nm_inv, v=T)
      if (length(i) > 0) {
         i=str2ind(i, nm_i)
         cat("The following inequalities are active after first stage
of zero crossing strategy and will be inverted:\\n", paste(nm_i[i], collapse="\\n"), "\\n", sep="")
         ipos=grep(">=", nm_i[i], v=T)
         ineg=grep("<=", nm_i[i], v=T)
         ui[i,]=-ui[i,]
         if (length(ipos)) {
            ipzc=str2ind(ipos, nm_izc)
            ipos=str2ind(ipos, nm_i)
            ci[ipos]=zc+mi_zc[ipzc,,drop=F]%*%mic
            nm_i[ipos]=sub(">=", "<=-", nm_i[ipos])
         }
         if (length(ineg)) {
            inzc=str2ind(ineg, nm_izc)
            ineg=str2ind(ineg, nm_i)
            ci[ineg]=zc+mi_zc[inzc,,drop=F]%*%mic
            nm_i[ineg]=sub("<=-", ">=", nm_i[ineg])
         }
         rownames(ui)=nm_i
         names(ci)=nm_i
      }
      # enforce new inequalities
      reopt=TRUE
      param=put_inside(res$par, ui, ci)
      if (any(is.na(param))) {
         if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
            # fatal error occurd, don't reoptimize
            warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
            param=res$param
            reopt=FALSE
         }
      } else if (!is.null(attr(param, "err")) && attr(param, "err")==0){
         # non fatal problem
         warning(paste("put_inside: ", attr(param, "mes"), collapse=""))
      }
      # reoptimize
      if (reopt) {
         res=opt_wrapper(measvec, fmn);
         if (any(is.na(res$par))) {
            obj2kvh(res, "optimization process informations", fkvh);
            stop("Second optimization failed")
         }
      }
   }
   param=res$par;
   names(param)=nm_par;
   obj2kvh(res, "optimization process informations", fkvh);
}
if (TIMEIT) {
   cat("postopt : ", date(), "\\n", sep="");
}
# active constraints
ine=abs(ui%*%param-ci)<1.e-10;
if (any(ine)) {
   obj2kvh(nm_i[ine], "active inequality constraints", fkvh);
}
#browser();
rcost=cumo_cost(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr);
rres=cumo_resid(param, cjac=T, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr);
obj2kvh(rcost, "final cost", fkvh);
names(rres$res)=c(nm_meas, nm_fmn);
obj2kvh(rres$res, "(simulated-measured)/sd_exp", fkvh);

# simulated measurements -> kvh
obj2kvh(cbind(value=jx_f$usimcumom,sd=1./sqrt(measinvvar)), "simulated unscaled labeling measurements", fkvh);
if (nb_sc > 0) {
   obj2kvh(cbind(value=jx_f$usimcumom*c(1.,param)[ir2isc],sd=1./sqrt(measinvvar)), "simulated scaled labeling measurements", fkvh);
}
if (nb_fmn) {
   obj2kvh(cbind(value=jx_f$fallnx[nm_fmn], sd=1./sqrt(invfmnvar)), "simulated flux measurements", fkvh);
}

# gradient -> kvh
gr=2*c(((jx_f$ures*c(measinvvar, invfmnvar))%tmm%jx_f$udr_dp));
names(gr)=nm_par;
obj2kvh(gr, "gradient vector", fkvh);
obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh);

if (fullsys) {
   nm_flist=nm_list
   nm_flist$rcumo=nm_cumo
   v=param2fl_x(param, cjac=F, nb_f, nm_flist, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAbr_f);
} else {
   v=param2fl_x(param, cjac=F, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, spAbr);
}
x=v$x;
if (fullsys) {
   names(x)=nm_cumo;
} else {
   names(x)=nm_rcumo;
}
obj2kvh(cumo2mass(x), "MID vector", fkvh);
# replace :i by :x1's
nm_mask=unlist(lapply(strsplit(names(x), ":"), function(v) {
   n=length(v);
   metab=paste(v[-n], collapse=":");
   mask=int2bit(v[n], len=clen[metab])
   c(metab, mask);
}))
nm_mask=t(matrix(nm_mask, nrow=2));
o=order(nm_mask[,1], nm_mask[,2]);
# replace 0 by x in mask
nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T);
nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#");
names(nm_mask)=names(x);
names(x)=nm_mask;
obj2kvh(x[o], "cumomer vector", fkvh);

fwrv=v$fwrv;
names(fwrv)=nm_fwrv;
#obj2kvh(fwrv, "fwd-rev flux vector", fkvh);

fallnx=v$fallnx;
names(fallnx)=nm_fallnx;

flnx=v$flnx;
names(flnx)=c(nm_fl);

if (sensitive=="mc") {
   if (TIMEIT) {
      cat("monte-ca: ", date(), "\\n", sep="");
   }
   # Monte-Carlo simulation in parallel way
   mc_inst=library(multicore, warn.conflicts=F, verbose=F, logical.return=T);
   invar=c(measinvvar, invfmnvar);
   simcumom=c(1.,param)[ir2isc]*jx_f$usimcumom;
   simfmn=f[nm_fmn];
   mc_sim=function(i) {
      # random measurement generation
      if (nb_meas) {
         meas_mc=rnorm(nb_meas, simcumom, 1./sqrt(measinvvar));
      } else {
         meas_mc=c();
      }
      if (nb_fmn) {
         fmn_mc=rnorm(nb_fmn, simfmn, 1./sqrt(invfmnvar));
      } else {
         fmn_mc=c();
      }
      #cat("imc=", i, "\\n", sep="");
      # minimization
      res=opt_wrapper(meas_mc, fmn_mc, trace=0);
      # return the solution
      return(list(cost=sum(jx_f$res*jx_f$res), par=res$par));
   }
   # parallel execution
   if (mc_inst) {
      mc_res=mclapply(1:nmc, mc_sim);
   } else {
      mc_res=lapply(1:nmc, mc_sim);
   }
   # strip failed mc iterations
   free_mc=sapply(mc_res, function(l) {if (class(l)=="character" || is.na(l$cost)) { ret=rep(NA, nb_param+1) } else { ret=c(l$cost, l$par) }; ret });
   if (length(free_mc)==0) {
      warning("Parallel exectution of Monte-Carlo simulations has failed.")
      free_mc=matrix(NA, nb_param+1, 0)
   }
   cost_mc=free_mc[1,];
   free_mc=free_mc[-1,,drop=F];
   # remove failed m-c iterations
   ifa=which(is.na(free_mc[1,]))
   if (length(ifa)) {
      if (ncol(free_mc)>length(ifa)) {
         warning("Some Monte-Carlo iterations failed.")
      }
      free_mc=free_mc[,-ifa,drop=F]
      cost_mc=cost_mc[-ifa]
      nmc=length(cost_mc)
   }
#browser();
   dimnames(free_mc)[[1]]=nm_par;
   cat("monte-carlo\\n", file=fkvh);
   indent=1;
   obj2kvh(nmc, "sample number", fkvh, indent);
   avaco=multicore:::detectCores()
   obj2kvh(avaco, "detected cores", fkvh, indent);
   avaco=max(1, avaco, na.rm=T)
   obj2kvh(min(avaco, options("cores")$cores, na.rm=T), "used cores", fkvh, indent);
   # cost section in kvh
   cat("\\tcost\\n", file=fkvh);
   indent=2;
   obj2kvh(mean(cost_mc), "mean", fkvh, indent);
   obj2kvh(median(cost_mc), "median", fkvh, indent);
   obj2kvh(sd(cost_mc), "sd", fkvh, indent);
   obj2kvh(sd(cost_mc)*100/mean(cost_mc), "rsd (%)", fkvh, indent);
   obj2kvh(quantile(cost_mc, c(0.025, 0.975)), "ci", fkvh, indent);
   # free parameters section in kvh
   cat("\\tfree parameters\\n", file=fkvh);
   indent=2;
   # param stats
   # mean
   obj2kvh(apply(free_mc, 1, mean), "mean", fkvh, indent);
   # median
   parmed=apply(free_mc, 1, median);
   obj2kvh(parmed, "median", fkvh, indent);
   # covariance matrix
   covmc=cov(t(free_mc));
   obj2kvh(covmc, "covariance", fkvh, indent);
   # sd
   sdmc=sqrt(diag(covmc));
   obj2kvh(sdmc, "sd", fkvh, indent);
   obj2kvh(sdmc*100/abs(param), "rsd (%)", fkvh, indent);
   # confidence intervals
   ci_mc=t(apply(free_mc, 1, quantile, probs=c(0.025, 0.975)));
   ci_mc=cbind(ci_mc, t(diff(t(ci_mc))));
   dimnames(ci_mc)[[2]][3]="length";
   obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent);
   obj2kvh((ci_mc-cbind(param, param, 0))*100/abs(param),
      "relative 95% confidence intervals (%)", fkvh, indent);
   
   # net-xch01 stats
   fallnx_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fallnx);
   fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fallnx;
   if (length(fallnx_mc)) {
      dimnames(fallnx_mc)[[1]]=nm_fallnx;
      # form a matrix output
      fallout=matrix(0, nrow=nrow(fallnx_mc), ncol=0);
      #cat("\\tall net-xch01 fluxes\\n", file=fkvh);
      # mean
#browser()
      fallout=cbind(fallout, mean=apply(fallnx_mc, 1, mean));
      #obj2kvh(apply(fallnx_mc, 1, mean), "mean", fkvh, indent);
      # median
      parmed=apply(fallnx_mc, 1, median);
      fallout=cbind(fallout, median=parmed);
      #obj2kvh(parmed, "median", fkvh, indent);
      # covariance matrix
      covmc=cov(t(fallnx_mc));
      dimnames(covmc)=list(nm_fallnx, nm_fallnx);
      #obj2kvh(covmc, "covariance", fkvh, indent);
      # sd
      sdmc=sqrt(diag(covmc));
      fallout=cbind(fallout, sd=sdmc);
      #obj2kvh(sdmc, "sd", fkvh, indent);
      fallout=cbind(fallout, "rsd (%)"=sdmc*100/abs(fallnx));
      #obj2kvh(sdmc*100/abs(fallnx), "rsd (%)", fkvh, indent);
      # confidence intervals
      ci_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)));
      ci_mc=cbind(ci_mc, t(diff(t(ci_mc))));
      ci_mc=cbind(ci_mc, (ci_mc-cbind(fallnx, fallnx, 0))*100/abs(fallnx));
      dimnames(ci_mc)[[2]]=c("ci 2.5%", "ci 97.5%", "ci 95% length", "rci 2.5% (%)", "rci 97.5% (%)", "rci 95% length (%)");
      #obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent);
      fallout=cbind(fallout, ci_mc);
      #obj2kvh((ci_mc-cbind(fallnx, fallnx, 0))*100/abs(fallnx),
      #   "relative 95% confidence intervals (%)", fkvh, indent);
      o=order(nm_fallnx);
      obj2kvh(fallout[o,,drop=F], "all net-xch01 fluxes", fkvh, indent);
      obj2kvh(covmc[o,o], "covariance of all net-xch01 fluxes", fkvh, indent);
      
      # fwd-rev stats
      fwrv_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, bp, fc)$fwrv);
      dimnames(fwrv_mc)[[1]]=nm_fwrv;
      fallout=matrix(0, nrow=nrow(fwrv_mc), ncol=0);
      #cat("\\tforward-reverse fluxes\\n", file=fkvh);
      # mean
      fallout=cbind(fallout, mean=apply(fwrv_mc, 1, mean))
      #obj2kvh(apply(fwrv_mc, 1, mean), "mean", fkvh, indent);
      # median
      parmed=apply(fwrv_mc, 1, median);
      fallout=cbind(fallout, median=parmed)
      #obj2kvh(parmed, "median", fkvh, indent);
      # covariance matrix
      covmc=cov(t(fwrv_mc));
      dimnames(covmc)=list(nm_fwrv, nm_fwrv);
      #obj2kvh(covmc, "covariance", fkvh, indent);
      # sd
      sdmc=sqrt(diag(covmc));
      fallout=cbind(fallout, sd=sdmc);
      #obj2kvh(sdmc, "sd", fkvh, indent);
      fallout=cbind(fallout, "rsd (%)"=sdmc*100/abs(fwrv));
      #obj2kvh(sdmc*100/abs(fwrv), "rsd (%)", fkvh, indent);
      # confidence intervals
      ci_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)));
      ci_mc=cbind(ci_mc, t(diff(t(ci_mc))));
      ci_mc=cbind(ci_mc, (ci_mc-cbind(fwrv, fwrv, 0))*100/abs(fwrv))
      dimnames(ci_mc)[[2]]=c("ci 2.5%", "ci 97.5%", "ci 95% length", "rci 2.5% (%)", "rci 97.5% (%)", "rci 95% length (%)");
      #obj2kvh(ci_mc, "95% confidence intervals", fkvh, indent);
      #obj2kvh((ci_mc-cbind(fwrv, fwrv, 0))*100/abs(fwrv),
      #   "relative 95% confidence intervals (%)", fkvh, indent);
      fallout=cbind(fallout, ci_mc);
      o=order(nm_fwrv);
      obj2kvh(fallout[o,,drop=F], "forward-reverse fluxes", fkvh, indent);
      obj2kvh(covmc[o,o], "covariance of forward-reverse fluxes", fkvh, indent);
   }
} else if (length(sensitive) && nchar(sensitive)) {
   warning(paste("Unknown sensitivity '", sensitive, "' method chosen.", sep=""));
}

if (TIMEIT) {
   cat("linstats: ", date(), "\\n", sep="");
}
# Linear method based on jacobian x_f
# reset fluxes and jacobians according to param
rres=cumo_resid(param, cjac=T, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr);
#grj=cumo_gradj(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr);
if (DEBUG) {
   library(numDeriv); # for numerical jacobian
   # numerical simulation
   rj=function(v, ...) {
      r=cumo_resid(v, cjac=F, ...);
      c(r$res);
   }
   #dr_dpn=jacobian(rj, param, method="Richardson", method.args=list(),
   #   nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr);
   dx_dfn=num_jacob(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, "fwrv2rAbcumo");
   dr_dp=num_fjacob(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAbr);
   # to compare with jx_f$dr_dp
   gr=grad(cumo_cost, param, method="Richardson", method.args=list(), nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr);
   grn=cumo_grad(param, nb_f, nm_list, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAbr);
#browser();
}

# covariance matrix of free fluxes
invcov=(jx_f$dr_dff)%tmm%(jx_f$dr_dff*c(measinvvar, invfmnvar));
covff=try(solve(invcov));
if (inherits(covff, "try-error")) {
   # matrix seems to be singular
   #svi=svd(invcov);
   #i=svi$d/svi$d[1]<1.e-14;
   #ibad=apply(abs(svi$v[, i, drop=F]), 2, which.max);
   #warning(paste("Inverse of covariance matrix is singular.\\nStatistically undefined fluxe(s) seems to be:\\n",
   #   paste(nm_ff[ibad], collapse="\\n"), "\\nFor more complete list, see '/linear stats/net-xch01 fluxes (sorted by name)' field in the result file.", sep=""));
   # regularize and inverse
   #covff=solve(invcov+diag(1.e-14*svi$d[1], nrow(invcov)));
   #covff=svi$u%*%diag(1./pmax(svi$d, svi$d[1]*1.e-14))%*%t(svi$v)
   svj=svd(jx_f$dr_dff)
   i=svj$d/svj$d[1]<1.e-14;
   if (all(!i)) {
      # we could not find very small d, take just the last
      i[length(i)]=T
   }
   ibad=apply(svj$v[, i, drop=F], 2, which.contrib);
   ibad=unique(unlist(ibad))
   warning(paste("Inverse of covariance matrix is numerically singular.\\nStatistically undefined fluxe(s) seems to be:\\n",
      paste(nm_ff[ibad], collapse="\\n"), "\\nFor more complete list, see the field '/linear stats/net-xch01 fluxes (sorted by name)'\\nin the result file.", sep=""));
   # "square root" of covariance matrix (to preserve numerical positive definitness)
   rtcov=(svj$u/sqrt(c(measinvvar, invfmnvar)))%*%(t(svj$v)/svj$d)
   covff=svj$v%*%(t(svj$u)/svj$d)%*%(svj$u/c(measinvvar, invfmnvar))%*%(t(svj$v)/svj$d)
} else {
   rtcov=chol(covff)
}
# standart deviations of free fluxes
#sdff=sqrt(diag(covff));
cat("linear stats\\n", file=fkvh);
#mtmp=cbind(param[1:nb_ff], sdff, sdff/abs(param[1:nb_ff]));
#dimnames(mtmp)[[2]]=c("val", "sd", "rsd");
#o=order(nm_ff);
#obj2kvh(mtmp[o,,drop=F], "val, sd, rsd free fluxes (sorted by name)", fkvh, indent=1);
#obj2kvh(covff, "covariance free fluxes", fkvh, indent=1);

# sd free+dependent net-xch01 fluxes
fl=c(param[1:nb_ff], flnx);
nm_flfd=c(nm_ff, nm_fl);
covfl=crossprod(rtcov%mmt%(rbind(diag(1., nb_ff), dfl_dff)));
dimnames(covfl)=list(nm_flfd, nm_flfd)
sdfl=sqrt(diag(covfl));
mtmp=cbind("value"=fl, "sd"=sdfl, "rsd"=sdfl/abs(fl));
rownames(mtmp)=nm_flfd;
o=order(nm_flfd);
obj2kvh(mtmp[o,,drop=F], "net-xch01 fluxes (sorted by name)", fkvh, indent=1);
obj2kvh(covfl[o, o], "covariance net-xch01 fluxes", fkvh, indent=1);

# sd of all fwd-rev
covf=crossprod(rtcov%mmt%(jx_f$df_dff));
sdf=sqrt(diag(covf));
mtmp=cbind(fwrv, sdf, sdf/abs(fwrv));
dimnames(mtmp)[[2]]=c("value", "sd", "rsd");
o=order(nm_fwrv);
obj2kvh(mtmp[o,], "fwd-rev fluxes (sorted by name)", fkvh, indent=1);
obj2kvh(covf, "covariance fwd-rev fluxes", fkvh, indent=1);

# select best defined flux combinations
s=svd(covfl);
# lin comb coeff matrix
l=t(s$u);
iord=apply(l,1,function(v){o=order(abs(v), decreasing=T); n=which(cumsum(v[o]**2)>=0.999); o[1:n[1]]})
dimnames(l)=dimnames(covfl);
nm_fl=dimnames(l)[[1]];
nb_fl=nrow(s$u);
comb=rep("", nb_fl);
for (j in 1:nb_fl) {
   cva=l[j,iord[[j]]]; # abs(factor) ordered j-th lin combination
   cva=sign(cva[1])*cva;
   cva=sprintf("%+.3g*", cva);
   cva[cva=="+1*"]="+";
   cva[cva=="-1*"]="-";
   cva[1]=substring(cva[1], 2);
   comb[j]=paste(cva, nm_fl[iord[[j]]], sep="", collapse="");
}
# from best to worst defined sd
sta=sqrt(abs(diag(l%*%covfl%mmt%(l))));
names(sta)=comb;
o=order(sta);
sta=sta[o];
obj2kvh(sta, "sd for uncorrelated net-xch01 linear combinations", fkvh, indent=1);
if (prof) {
   Rprof(NULL);
}
close(fkvh);
""");
    f.write("""
# write edge.netflux property
fedge=file("%(d)s/edge.netflux.%(org)s", "w");
cat("netflux (class=Double)\\n", sep="", file=fedge);
nm_edge=names(edge2fl)
cat(paste(nm_edge, fallnx[edge2fl], sep=" = "), sep="\\n" , file=fedge);

# write edge.xchflux property
fedge=file("%(d)s/edge.xchflux.%(org)s", "w")
flxch=paste(".x", substring(edge2fl, 4), sep="")
ifl=charmatch(flxch, substring(names(fallnx), 2))
cat("xchflux (class=Double)\\n", sep="", file=fedge)
cat(paste(nm_edge, fallnx[ifl], sep=" = "), sep="\\n" , file=fedge)
close(fedge)

if (TIMEIT) {
   cat("rend    : ", date(), "\\n", sep="");
}
"""%{
    "org": escape(org, "\\"),
    "d": escape(dirorg, "\\"),
});

    f.close();
    # make output files just readable to avoid later casual edition
    os.chmod(n_R, stat.S_IREAD);
