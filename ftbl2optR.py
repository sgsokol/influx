#!/usr/bin/env python

r"""
Transform an ftbl to R code which will solve an optimization of flux analysis
problem :math:`\arg \min_{\Theta} S`, where :math:`S=||\mathrm{Predicted}-\mathrm{Observed}||^{2}_{\Sigma}`
and :math:`\Theta` is a vector of parameters to fit: free fluxes (net+xch), scaling parameters and metabolite concentrations
pools.
Predicted vector is obtained from cumomer vector x (calculated from
free fluxes and divided in chunks according to the cumo weight) by
multiplying it by the measurement matrices, weighted by metabolite
pools (in case of pooling) and scale factor, boths coming
from ftbl file. Observed values vector xo is extracted from ftbl file.
it is composed of flux and cumomer measurements.
:math:`\Sigma^2`, covariance diagonal matrices sigma[flux|mass|label|peak|metab.pool]
is orginated from the ftbl file.

usage: ./ftbl2optR.py organism
where organism is the ftbl informative part of file name
(before .ftbl), e.g. organism.ftbl
after execution a file organism.R will be created.
If it already exists, it will be silently overwritten.
The system Afl*flnx=bfl is created from the ftbl file.

Important python variables:

Collections:
   * netan - (dict) ftbl structured content
   * tfallnx - (3-tuple[reac,["d"|"f"|"c"], ["net"|"xch"]] list)- total flux
    collection
   * measures - (dict) exp data
   * rAb - (list) reduced linear systems A*x_cumo=b by weight
   * scale - unique scale names
   * nrow - counts scale names
   * o_sc - ordered scale names
   * o_meas - ordered measurement types
   
File names (str):
   * n_ftbl (descriptor f_ftbl)
   * n_R (R code) (f)
   * n_fort (fortran code) (ff)

Counts:
   * nb_fln, nb_flx, nb_fl (dependent fluxes: net, xch, total), nb_ffn, nb_ffx (free fluxes)

Index translators:
   * fwrv2i - flux names to index in R:fwrv
   * cumo2i - cumomer names to index in R:x
   * ir2isc - mapping measurement rows indexes on scale index isc[meas]=ir2isc[meas][ir]

Vector names:
   * cumos (list) - names of R:x
   * o_mcumos - cumomers involved in measurements

Important R variables:

Scalars:
   * nb_w, nb_cumos, nb_fln, nb_flx, nb_fl (dependent or unknown fluxes),
   * nb_ffn, nb_ffx, nb_ff (free fluxes),
   * nb_fcn, nb_fcx, nb_fc (constrained fluxes),
   * nb_ineq, nb_param, nb_fmn

Name vectors:
   * nm_cumo, nm_fwrv, nm_fallnx, nm_fln, nm_flx, nm_fl, nm_par,
   * nm_ffn, nm_ffx,
   * nm_fcn, nm_fcx,
   * nm_mcumo, nm_fmn

Numeric vectors:
   * fwrv - all fluxes (fwd+rev)
   * x - all cumomers (weight1+weight2+...)
   * param - free flux net, free flux xch, scale label, scale mass, scale peak, metabolite concentrations
   * fcn, fcx, fc - constrained fluxes
   * bp - helps to construct the rhs of flux system
   * xi -cumomer input vector
   * fallnx - complete flux vector (constr+net+xch)
   * bc - helps to construct fallnx
   * li - inequality vector (mi%*%fallnx>=li)
   * ir2isc - measur row to scale vector replicator
   * ci - inequalities for param use (ui%*%param-ci>=0)
   * measvec - measurement vector
   * fmn

Matrices:
   * Afl, qrAfl, invAfl,
   * p2bfl - helps to construct the rhs of flux system
   * mf, md - help to construct fallnx
   * mi - inequality matrix (ftbl content)
   * ui - inequality matrix (ready for param use)
   * measmat - for measmat*x+memaone=vec of simulated not-yet-scaled measurements

Functions:
   * param2fl_x - translate param to flux and cumomer vector (initial approximation)
   * cumo_cost - cost function (khi2)
   * cumo_gradj - implicit derivative gradient
"""

# 2008-07-11 sokol: initial version
# 2009-03-18 sokol: interface homogenization for influx_sim package
# 2010-10-16 sokol: fortran code is no more generated, R Matrix package is used for sparse matrices.

if __name__ == "__main__":
    import sys
    import os
    import stat
    import time
    import copy
    import getopt
    #import pdb

    me=os.path.realpath(sys.argv[0])
    dirx=os.path.dirname(me)
    sys.path.append(os.path.dirname(me))
    me=os.path.basename(me)

    from tools_ssg import *
    import C13_ftbl

    try:
        import psyco
        psyco.full()
        #psyco.log()
        #psyco.profile()
    except ImportError:
        #print 'Psyco not installed, the program will just run slower'
        pass

    def usage():
        sys.stderr.write("usage: "+me+" [-h|--help] [--fullsys] [--emu] [--clownr] [--DEBUG] [--ropts ROPTS] network_name[.ftbl]\n")

    #<--skip in interactive session
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "fullsys", "emu", "clownr", "DEBUG", "ropts="])
    except getopt.GetoptError, err:
        #pass
        sys.stderr.write(str(err)+"\n")
        usage()
        sys.exit(1)

    fullsys=False
    DEBUG=False
    emu=False
    clownr=False
    ropts=['""']
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o=="--fullsys":
            fullsys=True
        elif o=="--emu":
            emu=True
        elif o=="--DEBUG":
            DEBUG=True
        elif o=="--clownr":
            clownr=True
        elif o=="--ropts":
            ropts=a.split("; ") if len(a) else ['""']
        else:
            #assert False, "unhandled option"
            # unknown options can come from shell
            # which passes all options both to python and R scripts
            # so just ignore unknown options
            #pass
            raise
    #aff("args", args);##
    #aff("opts", opts);##
    if len(args) != 1:
        usage()
        exit(1)
    org=os.path.basename(args[0])
    dirorg=os.path.dirname(args[0]) or '.'

    # cut .ftbl if any
    if org[-5:]==".ftbl":
        org=org[:-5]
    fullorg=os.path.join(dirorg, org)

    #-->
    #DEBUG=True
    import ftbl2code
    ftbl2code.DEBUG=DEBUG
    C13_ftbl.clownr=clownr
    #org="ex3"
    #org="PPP_exact"
    #DEBUG=True
    #if DEBUG:
    #    import pdb


    n_ftbl=fullorg+".ftbl"
    n_R=fullorg+".R"
    #n_fort=fullorg+".f"
    f_ftbl=open(n_ftbl, "r")
    try:
        os.chmod(n_R, stat.S_IWRITE)
    except:
        pass
    f=open(n_R, "w")

    # parse ftbl
    ftbl=C13_ftbl.ftbl_parse(f_ftbl)
    f_ftbl.close()

    # analyse network
    # reload(C13_ftbl)

    netan=C13_ftbl.ftbl_netan(ftbl)
    # prepare rcumo system
    if emu:
        rAb=C13_ftbl.rcumo_sys(netan, C13_ftbl.ms_frag_gath(netan))
    else:
        rAb=C13_ftbl.rcumo_sys(netan)

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, fullsys, emu, ropts)

    f.write("""
#browser()
# metabolite pools are : all (poolall) which is divided in free (poolf) and
# constrained (poolc)

# constrained pool
poolc=c(%(poolc)s)
nm_poolc=c(%(nm_poolc)s)
names(poolc)=nm_poolc

# starting values for free pool
poolf=c(%(poolf)s)
nm_poolf=c(%(nm_poolf)s)
names(poolf)=nm_poolf
nb_poolf=length(poolf)
nb_f$nb_poolf=nb_poolf

nm_poolall=c(nm_poolf, nm_poolc)
poolall=as.numeric(c(poolf, poolc))
names(poolall)=nm_poolall
pool=poolall

# extend param vector by free pools
if (nb_poolf > 0) {
#expp   param=c(param, log(poolf))
   param=c(param, poolf)
   nm_par=c(nm_par, nm_poolf)
   nb_param=length(param)
}
nm_list$par=nm_par

nm_list$poolf=nm_poolf
nm_list$poolc=nm_poolc
nm_list$poolall=nm_poolall

#browser()
if (nb_poolf > 0) {
#expp   # extend inequalities ui, ci by log(cupp) >= poolf >= log(clowp)
# extend inequalities ui, ci by cupp>= poolf >= clowp
   nb_row=nrow(ui)
   nb_col=ncol(ui)
   ui=cBind(ui, Matrix (0., nrow=nrow(ui), ncol=nb_poolf))
   ui=rBind(ui, Matrix (0., nrow=2*nb_poolf, ncol=ncol(ui)))
   ui[nb_row+1:nb_poolf,nb_col+1:nb_poolf]=diag(1., nb_poolf)
   ui[nb_row+nb_poolf+1:nb_poolf,nb_col+1:nb_poolf]=diag(-1., nb_poolf)
#expp   ci=c(ci, rep(log(clowp), nb_poolf))
#expp   ci=c(ci, rep(-log(cupp), nb_poolf))
   ci=c(ci, rep(clowp, nb_poolf))
   ci=c(ci, rep(-cupp, nb_poolf))
   nm_i=c(nm_i, paste(nm_poolf, ">=", clowp, sep=""))
   nm_i=c(nm_i, paste(nm_poolf, "<=", cupp, sep=""))
   rownames(ui)=nm_i
   names(ci)=nm_i
   colnames(ui)=nm_par
}
# prepare metabolite measurements
nb_poolm=%(nb_poolm)d
nb_f$nb_poolm=nb_poolm
nm_poolm=c(%(nm_poolm)s)
nm_list$poolm=nm_poolm

# measured values
vecpoolm=c(%(v_poolm)s)
names(vecpoolm)=nm_poolm

# inverse of variance for flux measurements
poolmdev=c(%(poolmdev)s)

# simulated metabolite measurements are calculated as
# measmatpool*poolall=>poolm
measmatpool=Matrix(0., nrow=nb_poolm, ncol=length(poolall))
dimnames(measmatpool)=list(nm_poolm, nm_poolall)
i=matrix(1+c(%(imeasmatpool)s), ncol=2, byrow=T)
measmatpool[i]=1.

# gather all measurement information
measurements=list(
   vec=list(labeled=measvec, flux=fmn, pool=vecpoolm),
   dev=list(labeled=measdev, flux=fmndev, pool=poolmdev),
   mat=list(labeled=measmat, flux=ifmn, pool=measmatpool),
   one=list(labeled=memaone)
)
nm_resid=c(nm_meas, nm_fmn, nm_poolm)
nm_list$resid=nm_resid
"""%{
    "poolf": join(", ", (-netan["met_pools"][m] for m in netan["vpool"]["free"])),
    "nm_poolf": join(", ", netan["vpool"]["free"], '"pf:', '"'),
    "poolc": join(", ", (netan["met_pools"][m] for m in netan["vpool"]["constrained"])),
    "nm_poolc": join(", ", netan["vpool"]["constrained"], '"pc:', '"'),
    "nb_poolm": len(netan["metab_measured"]),
    "nm_poolm": join(", ", netan["metab_measured"].keys(), '"pm:', '"'),
    "v_poolm": join(", ", (item["val"] for item in netan["metab_measured"].values())).replace("nan", "NA"),
    "poolmdev": join(", ", (item["dev"] for item in netan["metab_measured"].values())),
    "imeasmatpool": join(", ", valval((ir,netan["vpool"]["all2i"][m]) for (ir, ml) in enumerate(netan["metab_measured"].keys()) for m in ml.split("+"))),
})
    f.write("""
if (TIMEIT) {
   cat("preopt  : ", date(), "\\n", sep="", file=fclog)
}
#browser()
names(param)=nm_par
# prepare series of starting points
if (nchar(fseries) > 0) {
   pstart=as.matrix(read.table(fseries, header=T, row.n=1, sep="\\t"))
   # skip parameters (rows) who's name is not in nm_par
   i=rownames(pstart) %in% nm_par
   if (!any(i)) {
      cat("Option --fseries is used but no free parameter with known name is found.", "\\n", sep="", file=fcerr)
      stop()
   }
   pstart=pstart[i,,drop=F]
   cat("Using starting values form '", fseries, "' for the following free parameters:\\n", paste(rownames(pstart), collapse="\\n"), "\\n", sep="", file=fclog)
#expp   # take log of free pools
#expp   i=grep("^pf:", rownames(pstart))
#expp   pstart[i,]=log(pstart[i,])
} else {
   pstart=as.matrix(param)
}
nm_pseries=rownames(pstart)

nseries=ncol(pstart)
if (is.null(nseries) || nseries==0) {
   cat("No starting values in the series file '"%s+%fseries%s+%"'.", "\\n", sep="", file=fcerr)
   stop()
}

# prepare series indexes
if (nchar(iseries) > 0) {
   iseries=unique(as.integer(eval(parse(t="c("%s+%iseries%s+%")"))))
   if (initrand) {
      if (nchar(fseries)==0) {
         nseries=max(iseries)
         pstart=matrix(0., nb_param, nseries)
         dimnames(pstart)=list(nm_par, "V" %s+% iseq(nseries))
      }
      pstart[]=runif(nb_param*max(iseries))
      # take log of free pools
      i=grep("^pf:", rownames(pstart))
#expp      pstart[i,]=log(pstart[i,])
   }
   iseries=iseries[iseries<=nseries]
   pstart=pstart[,iseries, drop=F]
   nseries=ncol(pstart)
} else {
   iseries=iseq(nseries)
}

pres=matrix(NA, nb_param, nseries)
rownames(pres)=nm_par
colnames(pres)=colnames(pstart)
costres=rep.int(NA, nseries)
""")
    f.write("""
# formated output in kvh file
fkvh_saved="%s_res.kvh"
"""%escape(fullorg, "\\"))
    f.write("""
for (irun in iseq(nseries)) {
   param[nm_pseries]=pstart[nm_pseries, irun]
   # prepare kvh file name
   if (nseries > 1) {
      runsuf="." %s+% colnames(pstart)[irun]
   } else {
      runsuf=""
   }
   if (length(nseries) > 0) {
      cat("Starting point", runsuf, "\\n", sep="", file=fclog)
   }
   fkvh=file(substring(fkvh_saved, 1, nchar(fkvh_saved)-4) %s+% runsuf %s+% ".kvh", "w");

   # remove zc inequalities from previous interations
   izc=grep("^zc ", nm_i)
   if (length(izc)) {
      ui=ui[-izc,,drop=F]
      ci=ci[-izc]
      nm_i=rownames(ui)
   }
   # check if initial approximation is feasible
   ineq=as.numeric(ui%*%param-ci)
   names(ineq)=rownames(ui)
   param_old=param
   if (any(ineq <= -1.e-10)) {
      cat("The following ", sum(ineq<= -1.e-10), " ineqalities are not respected at starting point", runsuf, ":\\n", sep="", file=fclog)
      i=ineq[ineq<= -1.e-10]
      cat(paste(names(i), i, sep="\\t", collapse="\\n"), "\\n", sep="", file=fclog)
      # put them inside
      capture.output(param <- put_inside(param, ui, ci), file=fclog)
      if (any(is.na(param))) {
         if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
            # fatal error occured
            cat("put_inside", runsuf, ": ", attr(param, "mes"), "\\n",
               file=fcerr, sep="")
            close(fkvh)
            next;
         }
      } else if (!is.null(attr(param, "err")) && attr(param, "err")==0) {
         # non fatal problem
         cat(paste("put_inside: ", attr(param, "mes"), collapse=""), "\\n", file=fcerr)
      }
   }

   # prepare zero crossing strategy
   # inequalities to keep sens of net flux on first call to opt_wrapper()
   # if active they are removed on the second call to opt_wrapper()
   # and finaly all zc constraints are relaxed on the last call to opt_wrapper()
   fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fallnx
   mi_zc=NULL
   li_zc=NULL
   if (zerocross && length(grep("[df].n.", nm_fallnx))>0) {
      # add lower limits on [df].net >= zc for positive net fluxes
      # and upper limits on [df].net <= -zc for negative net fluxes
      nm_izc=c()
      ipos=names(which(fallnx[grep("[df].n.", nm_fallnx)]>=0.))
      ineg=names(which(fallnx[grep("[df].n.", nm_fallnx)]<0.))
      mi_zc=matrix(0, nrow=length(ipos)+length(ineg), ncol=nb_fallnx)
      colnames(mi_zc)=nm_fallnx
      if (length(ipos)) {
         nm_izc=c(nm_izc, paste("zc ", ipos, ">=", zc, sep=""))
         mi_zc[(1:length(ipos)),ipos]=diag(1., length(ipos))
      }
      if (length(ineg)) {
         nm_izc=c(nm_izc, paste("zc ", ineg, "<=", -zc, sep=""))
         mi_zc[length(ipos)+(1:length(ineg)),ineg]=diag(-1., length(ineg))
      }
      rownames(mi_zc)=nm_izc
      li_zc=rep(zc, length(nm_izc)) # that's ok for both pos and neg constraints
      ui_zc=cBind(mi_zc%*%(md%*%invAfl%*%p2bfl+mf),
         Matrix(0., nrow=nrow(mi_zc), ncol=nb_sc))
      if (nb_fgr > 0) {
         ui_zc=cBind(ui_zc, mi_zc%*%((md%*%invAfl%*%g2bfl)+mg*nb_f$mu))
      } else if (nb_poolf > 0) {
         ui_zc=cBind(ui_zc, Matrix(0., nrow=nrow(mi_zc), ncol=nb_poolf))
      }
      ci_zc=li_zc-mi_zc%*%mic
      # remove constant inequalities
      if (ncol(ui_zc)) {
         zi=apply(ui_zc,1,function(v){return(max(abs(v))<=1.e-14)})
      } else {
         # remove all flux inequalities as there is no free params
         zi=rep(TRUE, nrow(ui_zc))
      }

      if (all(ci_zc[zi]<=1.e-10)) {
         ui_zc=ui_zc[!zi,,drop=F]
         ci_zc=ci_zc[!zi]
         nm_izc=nm_izc[!zi]
         mi_zc=mi_zc[!zi,,drop=F]
      } else {
         cat("The following constant inequalities are not satisfied:\\n", file=fcerr)
         cat(nm_izc[zi][ci_zc[zi]>1.e-14], sep="\\n", file=fcerr)
      }

      # remove redundant/contradictory inequalities
      nb_zc=nrow(ui_zc)
      nb_i=nrow(ui)
      ired=c()
      if (nb_zc > 0) {
         for (i in 1:nb_zc) {
            nmqry=nm_izc[i]
            for (j in 1:nb_i) {
               if ((max(abs(ui[j,]-ui_zc[i,]))<1.e-10 ||
                  max(abs(ui[j,]+ui_zc[i,]))<1.e-10) &&
                  abs(abs(ci[j])-abs(ci_zc[i]))<=1.e-2) {
#browser()
                  # redundancy
                  cat("inequality '", nmqry, "' redundant or contradictory with '", nm_i[j], "' is removed.\\n", sep="", file=fclog)
                  ired=c(ired, i)
                  break
               }
            }
         }
      }
      if (!is.null(ired)) {
         # remove all ired inequalities
         ui_zc=ui_zc[-ired,,drop=F]
         ci_zc=ci_zc[-ired]
         nm_izc=nm_izc[-ired]
         mi_zc=mi_zc[-ired,,drop=F]
      }
      if (nrow(ui_zc)) {
         # add zc inequalities
         ui=rBind(ui, ui_zc)
         ci=c(ci, ci_zc)
         nm_i=c(nm_i, nm_izc)
         mi=rBind(mi, mi_zc)
      }
   #print(ui)
   #print(ci)
   }

   # set initial scale values to sum(measvec*simvec/dev**2)/sum(simvec**2/dev**2)
   # for corresponding measurements
   vr=param2fl_x(param, cjac=FALSE, jx_f, nb_f, nm_list, nb_x, invAfl, p2bfl, g2bfl, bp, fc, xi, spa, emu, pool, measurements, ipooled)
   if (!is.null(vr$err) && vr$err) {
      cat("param2fl_x", runsuf, ": ", vr$mes, "\\n", file=fcerr, sep="")
      close(fkvh)
      next
   }
   jx_f=vr$jx_f
   if (nb_sc > 0) {
      if (optimize) {
         simvec=jx_f$usimcumom
         measinvvar=1./measurements$dev$labeled**2
         ms=measvec*simvec*measinvvar
         ss=simvec*simvec*measinvvar
         # get only valid measurements
         iva=!is.na(ms)
         for (i in nb_ff+1:nb_sc) {
            im=(ir2isc==(i+1)) & iva
            if (sum(im) < 2) {
               cat(sprintf("scaling: no sufficient valid data for scaling factor '%s'", nm_par[i]), "\\n", sep="", file=fcerr)
               stop()
            }
            param[i]=sum(ms[im])/sum(ss[im])
         }
      } else {
         # if no optimization, set all scaling params to 1.
         param[nb_ff+1:nb_sc]=1.
      }
   }

   # prepare flux index conversion
   ifwrv=1:nb_fwrv
   names(ifwrv)=nm_fwrv
   ifl_in_fw=ifwrv[paste("fwd", substring(c(nm_fln, nm_flx), 4), sep="")]
   iff_in_fw=ifwrv[paste("fwd", substring(c(nm_ffn, nm_ffx), 4), sep="")]
   if (nb_fgr > 0) {
      ifg_in_fw=ifwrv[paste("fwd", substring(nm_fgr, 4), sep="")]
   } else {
      ifg_in_fw=numeric(0)
   }
   # index couples for jacobian df_dfl, df_dffd
   cfw_fl=crv_fl=cBind(ifl_in_fw, 1:nb_fl)
   cfw_ff=crv_ff=cBind(iff_in_fw, 1:nb_ff)
   cfw_fg=crv_fg=cBind(ifg_in_fw, tail(1:(nb_ff+nb_fgr), nb_fgr))
   crv_fl[,1]=(nb_fwrv/2)+crv_fl[,1]
   crv_ff[,1]=(nb_fwrv/2)+crv_ff[,1]
   crv_fg[,1]=(nb_fwrv/2)+crv_fg[,1]
   
   # store it in nb_f
   nb_f=append(nb_f, list(cfw_fl=cfw_fl, crv_fl=crv_fl, cfw_ff=cfw_ff,
      crv_ff=crv_ff, cfw_fg=cfw_fg, crv_fg=crv_fg))

   # see if there are any active inequalities at starting point
   ineq=as.numeric(ui%*%param-ci)
   names(ineq)=rownames(ui)
   if (any(abs(ineq)<=1.e-10)) {
      cat("The following ", sum(abs(ineq)<=1.e-10), " ineqalitie(s) are active at starting point", runsuf, ":\\n",
         paste(names(ineq[abs(ineq)<=1.e-10]), collapse="\\n"), "\\n", sep="", file=fclog)
   }
""")

    f.write("""
   if (TIMEIT) {
      cat("kvh init: ", date(), "\\n", sep="", file=fclog)
   }

"""%{
    "fullorg": escape(fullorg, "\\"),
})
    # main part: call optimization
    f.write("""
   cat("influx\\n", file=fkvh)
   cat("\\tversion\\t", vernum, "\\n", file=fkvh, sep="")
   # save options of command line
   cat("\\truntime options\\n", file=fkvh)
   cat("\\t\\t%s\\n", file=fkvh)
   """%join("\n\t\t", ropts)[1:-1])
    f.write("""
   obj2kvh(R.Version(), "R.Version", fkvh, indent=1)
   cat("\\tR command line\\n", file=fkvh)
   obj2kvh(opts, "opts", fkvh, indent=2)

   # resume system sizes
   obj2kvh(nb_sys, "system sizes", fkvh)

   # save initial flux and cumomer distribution
   cat("starting point\\n", file=fkvh)
   names(param)=nm_par
   obj2kvh(param, "starting free parameters", fkvh, indent=1)

   x=vr$x
   names(x)=nm_x
   obj2kvh(cumo2mass(x), "starting MID vector", fkvh, indent=1)
   # replace :i by #x1's
   nm_mask=t(sapply(strsplit(names(x), "[:+]"), function(v) {
      n=length(v)
      metab=v[1]
      mask=int2bit(v[2], len=clen[metab])
      c(metab, mask, if (n==3) v[3] else NULL)
   }))
   o=order(nm_mask[,1], nm_mask[,2])
   # replace 0 by x in mask
   nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T)
   tmp=paste(nm_mask[,1], nm_mask[,2], sep="#")
   if (emu) {
      nm_mask=paste(tmp, nm_mask[,3], sep="+")
   } else {
      nm_mask=tmp
   }
   names(nm_mask)=names(x)
   names(x)=nm_mask
   obj2kvh(x[o], "starting cumomer vector", fkvh, indent=1)

   fwrv=vr$fwrv
   n=length(fwrv)
   names(fwrv)=nm_fwrv
   obj2kvh(fwrv, "starting fwd-rev flux vector", fkvh, indent=1)

   f=vr$fallnx
   n=length(f)
   names(f)=nm_fallnx
   obj2kvh(f, "starting net-xch01 flux vector", fkvh, indent=1)

   rres=cumo_resid(param, cjac=TRUE, jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
   jx_f=rres$jx_f
   names(rres$res)=nm_resid
   o=order(names(rres$res))
   obj2kvh(rres$res[o], "starting residuals (simulated-measured)/sd_exp", fkvh, indent=1)

   rcost=cumo_cost(param, jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
   obj2kvh(rcost, "starting cost value", fkvh, indent=1)

   obj2kvh(Afl, "flux system (Afl)", fkvh, indent=1)
   btmp=as.numeric(p2bfl%*%param[1:nb_f$nb_ff]+bp)
   names(btmp)=dimnames(Afl)[[1]]
   obj2kvh(btmp, "flux system (bfl)", fkvh, indent=1)
   names(measvec)=nm_meas
   obj2kvh(measvec, "labeled measurement vector", fkvh, indent=1)

   #cat("mass vector:\\n", file=fclog)
   #print_mass(x)

   names(param)=nm_par
""")
    f.write("""
   control_ftbl=list(%(ctrl_ftbl)s)
"""%{
    "ctrl_ftbl": join(", ", (k[8:]+"="+str(v) for (k,v) in netan["opt"].iteritems() if k.startswith("optctrl_"))),
})
    f.write("""
   if (optimize) {
      # check if at starting position all fluxes can be resolved
#browser()
      qrj=qr(jx_f$dr_dff, LAPACK=T)
      d=diag(qrj$qr)
      qrj$rank=sum(abs(d)>abs(d[1])*1.e-14)
      if (qrj$rank) {
         nm_uns=nm_ff[qrj$pivot[-(1:qrj$rank)]]
      } else {
         nm_uns=nm_ff
      }
      if (qrj$rank < nb_ff && !(least_norm || method!="nlsic")) {
         library(MASS)
         # Too bad. The jacobian of free fluxes is not of full rank.
         dimnames(jx_f$dr_dff)[[2]]=c(nm_ffn, nm_ffx)
         fname="dbg_dr_dff_singular" %s+% runsuf %s+% ".csv"
         write.matrix(formatC(jx_f$dr_dff, 15), file=fname, sep="\\t")
         cat(paste("Provided measurements (isotopomers and fluxes) are not sufficient to resolve all free fluxes.\\nUnsolvable fluxes may be:
", paste(nm_uns, sep=", ", collapse=", "),
            "\\nJacobian dr_dff is dumped in " %s+% fname, sep=""),
            "\\n", file=fcerr)
         close(fkvh)
         next
      }
      if (TIMEIT) {
         cat("optim   : ", date(), "\\n", sep="", file=fclog)
      }
      # pass control to the chosen optimization method
      capture.output(res <- opt_wrapper(measurements, jx_f), file=fclog)
      if ((!is.null(res$err) && res$err) || is.null(res$par)) {
         cat("first optimization pass", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
         res$par=rep(NA, length(param))
         res$cost=NA
      } else if (!is.null(res$mes) && nchar(res$mes)) {
         cat("first optimization pass", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
      }
      jx_f=res$retres$jx_f
      if (any(is.na(res$par))) {
         res$retres$jx_f=NULL # to avoid writing of huge data
         obj2kvh(res, "failed first pass optimization process information", fkvh)
         cat("Optimization failed", runsuf, "\\n", file=fcerr, sep="")
         close(fkvh)
         next
      }
      param=res$par
      res_save=res
#browser()
      if (zerocross && !is.null(mi_zc)) {
         # inverse active "zc" inequalities
         nm_inv=names(which((ui%*%res$par-ci)[,1]<=1.e-10))
         i=grep("^zc ", nm_inv, v=T)
         if (length(i) > 0) {
            i=str2ind(i, nm_i)
            cat("The following inequalities are active after first pass
of zero crossing strategy and will be inverted", runsuf, ":\\n", paste(nm_i[i], collapse="\\n"), "\\n", sep="", file=fclog)
            ipos=grep(">=", nm_i[i], v=T)
            ineg=grep("<=", nm_i[i], v=T)
            ui[i,]=-ui[i,,drop=F]
            if (length(ipos)) {
               ipzc=str2ind(ipos, nm_izc)
               ipos=str2ind(ipos, nm_i)
               ci[ipos]=as.numeric(zc+mi_zc[ipzc,,drop=F]%*%mic)
               nm_i[ipos]=sub(">=", "<=-", nm_i[ipos])
            }
            if (length(ineg)) {
               inzc=str2ind(ineg, nm_izc)
               ineg=str2ind(ineg, nm_i)
               ci[ineg]=as.numeric(zc+mi_zc[inzc,,drop=F]%*%mic)
               nm_i[ineg]=sub("<=-", ">=", nm_i[ineg])
            }
            rownames(ui)=nm_i
            names(ci)=nm_i
            # enforce new inequalities
            reopt=TRUE
            capture.output(param <- put_inside(res$par, ui, ci), file=fclog)
            if (any(is.na(param))) {
               if (!is.null(attr(param, "err")) && attr(param, "err")!=0) {
                  # fatal error occured, don't reoptimize
                  cat(paste("put_inside", runsuf, ": ", attr(param, "mes"), "\\n", collapse=""), file=fcerr)
                  param=res_save$par
                  reopt=FALSE
               }
            } else if (!is.null(attr(param, "err")) && attr(param, "err")==0){
               # non fatal problem
               cat(paste("put_inside", runsuf, ": ", attr(param, "mes"), "\\n", collapse=""), file=fcerr)
            }
            # reoptimize
            if (reopt) {
               cat("Second zero crossing pass", runsuf, "\\n", sep="", file=fclog)
               capture.output(res <- opt_wrapper(measurements, jx_f), file=fclog)
               jx_f=res$retres$jx_f
               if (res$err || is.null(res$par)) {
                  cat("second zero crossing pass: ", res$mes, "\\n", sep="", file=fcerr)
                  res$par=rep(NA, length(param))
                  res$cost=NA
               } else if (!is.null(res$mes) && nchar(res$mes)) {
                  cat("second zero crossing pass", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
               }
               if(!res$err && !is.null(res$par) && !any(is.na(res$par))) {
                  param=res$par
                  res_save=res
               }
               if (any(is.na(res$par))) {
                  res$retres$jx_f=NULL # to avoid writing of huge data
                  obj2kvh(res, "failed second pass optimization process information", fkvh)
                  cat("Second zero crossing pass failed. Keep free parameters from previous pass", runsuf, "\\n", file=fcerr, sep="")
               }
            }
            # last pass, free all zc constraints
            i=grep("^zc ", nm_i)
            if (length(i) > 0) {
               ui=ui[-i,,drop=F]
               ci=ci[-i]
               nm_i=nm_i[-i]
               cat("Last zero crossing pass (free of zc constraints)", runsuf, "\\n", sep="", file=fclog)
               capture.output(res <- opt_wrapper(measurements, jx_f), file=fclog)
               jx_f=res$retres$jx_f
               if (res$err || is.null(res$par)) {
                  cat("last zero crossing (free of zc)", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
                  res$par=rep(NA, length(param))
                  res$cost=NA
               } else if (!is.null(res$mes) && nchar(res$mes)) {
                  cat("last zero crossing (free of zc)", runsuf, ": ", res$mes, "\\n", sep="", file=fcerr)
               }
               if(!res$err && !is.null(res$par) && !any(is.na(res$par))) {
                  param=res$par
                  res_save=res
               }
               if (any(is.na(res$par))) {
                  res$retres$jx_f=NULL # to avoid writing of huge data
                  obj2kvh(res, "failed last pass optimization process information", fkvh)
                  cat("Last zero crossing pass failed. Keep free parameters from previous passes", runsuf, "\\n", file=fcerr, sep="")
               }
            }
         } else {
            cat("After the first optimization, no zero crossing equality was activated. So no reptimization", runsuf, "\\n", sep="", file=fclog)
         }
      } # end if zero crossing
      res=res_save
      param=res$par
      names(param)=nm_par
      if (excl_outliers != F) {
         # detect outliers
         iva=!is.na(res$res)
         iout=which(rz.pval.bi(res$res) <= excl_outliers | !iva)
         #cat("iout=", iout, "\\n", file=fclog)
         if (length(iout)) {
            measurements$outlier=iout
            cat("Excluded outliers at p-value ", excl_outliers, ":\\n",
               paste(nm_resid[iout], res$res[iout], sep="\\t", collapse="\\n"), "\\n", sep="", file=fclog)
            jx_f_save=jx_f
            capture.output(resout <- opt_wrapper(measurements, jx_f), file=fclog)
            if (resout$err || is.null(resout$par)) {
               cat("wo outliers: ", resout$mes, "\\n", sep="", file=fcerr)
            } else if (!is.null(resout$mes) && nchar(resout$mes)) {
               cat("wo outliers: ", resout$mes, "\\n", sep="", file=fcerr)
            }
            if (any(is.na(resout$par))) {
               cat("Optimization with outliers excluded has failed", runsuf, "\\n", file=fcerr, sep="")
               # continue without outlier exclusion
               measurements$outlier=NULL
               jx_f=jx_f_save
            } else {
               res=resout
               param=res$par
               names(param)=nm_par
               obj2kvh(nm_resid[iout], "excluded outliers", fkvh)
               jx_f=resout$retres$jx_f
            }
         } else {
            cat("Outlier exclusion at p-value "%s+%excl_outliers%s+%" has been requested but no outlier was detected at this level.", "\\n", sep="", file=fcerr)
         }
      }
      res$jacobian=res$retres$jx_f=NULL # to avoid writing huge data
      obj2kvh(res, "optimization process information", fkvh)
   }
   if (TIMEIT) {
      cat("postopt : ", date(), "\\n", sep="", file=fclog)
   }
   # active constraints
   ine=as.numeric(abs(ui%*%param-ci))<1.e-10
   if (any(ine)) {
      obj2kvh(nm_i[ine], "active inequality constraints", fkvh)
   }
#expp   poolall[nm_poolf]=exp(param[nm_poolf])
   poolall[nm_poolf]=param[nm_poolf]

#browser()
   rres=cumo_resid(param, cjac=T, jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
   jx_f=rres$jx_f
   rcost=cumo_cost(param, jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)

   pres[,irun]=param
   costres[irun]=rcost
   obj2kvh(rcost, "final cost", fkvh)
   if (is.null(measurements$outlier) || length(measurements$outlier)==0) {
      names(rres$res)=nm_resid
   } else {
      names(rres$res)=nm_resid[-measurements$outlier]
   }
   o=order(names(rres$res))
   obj2kvh(rres$res[o], "(simulated-measured)/sd_exp", fkvh)

   # simulated measurements -> kvh
   obj2kvh(cBind(value=jx_f$usimcumom,sd=measurements$dev$labeled), "simulated unscaled labeling measurements", fkvh)
   if (nb_sc > 0) {
      val=jx_f$usimcumom*c(1.,param)[ir2isc]
      names(val)=nm_meas
      obj2kvh(cBind(value=val,sd=measurements$dev$labeled), "simulated scaled labeling measurements", fkvh)
   }
   if (nb_fmn) {
      obj2kvh(cBind(value=jx_f$fallnx[nm_fmn], sd=measurements$dev$flux), "simulated flux measurements", fkvh)
   }
   if (nb_poolm) {
      obj2kvh(cBind(value=measurements$mat$pool%*%poolall, sd=measurements$dev$pool), "simulated metabolite pool measurements", fkvh)
   }

   # gradient -> kvh
   gr=2*as.numeric(crossprod(jx_f$res, jx_f$jacobian))
   names(gr)=nm_par
   obj2kvh(gr, "gradient vector", fkvh)
   colnames(jx_f$udr_dp)=nm_par
   obj2kvh(jx_f$udr_dp, "jacobian dr_dp (without 1/sd_exp)", fkvh)
   
   # generalized inverse of non reduced jacobian
   svj=svd(jx_f$udr_dp)
   invj=svj$v%*%(t(svj$u)/svj$d)
   dimnames(invj)=rev(dimnames(jx_f$udr_dp))
   obj2kvh(invj, "generalized inverse of jacobian dr_dp (without 1/sd_exp)", fkvh)

   if (fullsys) {
      nm_flist=nm_list
      nm_flist$rcumo=nm_cumo
      nm_flist$rcumo_in_cumo=match(nm_rcumo, nm_cumo)
      nb_f$cumos=nb_cumos""")
    f.write("""
      nm_xi_f=c(%s)
      xi_f=c(%s)"""%(join(", ", netan["cumo_input"].keys(), '"', '"'),
      join(", ", netan["cumo_input"].values())))
    f.write("""
      names(xi_f)=nm_xi_f
      nm_flist$xi=nm_xi_f
      v=param2fl_x(param, cjac=F, jx_f, nb_f, nm_flist, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi_f, spAbr_f, emu=F, pool, measurements, ipooled)
   } else {
      v=param2fl_x(param, cjac=F, jx_f, nb_f, nm_list, nb_x, invAfl, p2bfl, g2bfl, bp, fc, xi, spa, emu, pool, measurements, ipooled)
   }
   jx_f=v$jx_f
   # set final variable depending on param
   x=v$x
   if (fullsys) {
      names(x)=nm_cumo
   } else {
      names(x)=nm_x
   }

   # write some info in result kvh
   obj2kvh(cumo2mass(x), "MID vector", fkvh)
   # replace decimal :i by binary :x1's
   nm_mask=t(sapply(strsplit(names(x), "[:+]"), function(v) {
      metab=v[1]
      v[2]=int2bit(v[2], len=clen[metab])
      v
   }))
   o=order(nm_mask[,1], nm_mask[,2])
   # replace 0 by x in mask
   nm_mask[,2]=gsub("0", "x", nm_mask[,2], fixed=T)
   if (emu) {
      nm_mask=paste(paste(nm_mask[,1], nm_mask[,2], sep="#"), nm_mask[,3], sep="+M")
   } else {
      nm_mask=paste(nm_mask[,1], nm_mask[,2], sep="#")
   }
   names(x)=nm_mask
   obj2kvh(x[o], "labeling vector", fkvh)

   fwrv=v$fwrv
   names(fwrv)=nm_fwrv
   #obj2kvh(fwrv, "fwd-rev flux vector", fkvh)

   fallnx=v$fallnx
   names(fallnx)=nm_fallnx

   flnx=v$flnx
   names(flnx)=c(nm_fl)

   fgr=fallnx[nm_fgr]

   # keep last jx_f in jx_f_last
   jx_f_last=jx_f
   while (sensitive=="mc") {
      if (TIMEIT) {
         cat("monte-ca: ", date(), "\\n", sep="", file=fclog)
      }
      if(set_seed) {
         set.seed(seed)
      }
      # Monte-Carlo simulation in parallel way (if asked and possible)
      simcumom=c(1.,param)[ir2isc]*jx_f$usimcumom
      simfmn=f[nm_fmn]
      simpool=as.numeric(measurements$mat$pool%*%poolall)
      #.Platform$OS.type="bidon"
      if (np > 1L) {
         # parallel execution
         # prepare cluster
         if (.Platform$OS.type=="unix") {
            type="FORK"
            nodes=np
         } else {
            type="PSOCK"
            nodes=rep("localhost", np)
         }
         #seeds=sample(1L:10000L, nmc)
         cl=makeCluster(nodes, type)
         if (.Platform$OS.type!="unix") {
            if (TIMEIT) {
               cat("cl init : ", date(), "\\n", sep="", file=fclog)
            }
            clusterEvalQ(cl, c(require(bitops), require(nnls), require(Matrix)))
""")
    f.write("""
            if (TIMEIT) {
               cat("cl sourc: ", date(), "\\n", sep="", file=fclog)
            }
            clusterEvalQ(cl, c(source("%(dirx)s/tools_ssg.R"), source("%(dirx)s/nlsic.R")))
"""%{"dirx": escape(dirx, "\\")})
    f.write("""
            if (TIMEIT) {
               cat("cl expor: ", date(), "\\n", sep="", file=fclog)
            }
            clusterExport(cl, c("nb_ff", "fcerr", "lsi_fun", "nm_ff", "nm_fmn", "dfm_dff", "cumo_jacob", "ind_bx", "fx2jr", "trisparse_solv", "fwrv2Abr", "pool", "ir2isc", "ipooled", "emu", "Heaviside", "nm_fwrv", "df_dffp", "DEBUG", "fallnx2fwrv", "fc", "dfcg2fallnx", "g2bfl", "bp", "p2bfl", "c2bfl", "invAfl", "param2fl", "nb_rcumos", "nm_list", "nb_f", "xi", "spa", "param2fl_x", "is.diff", "cumo_resid", "ui", "ci", "nlsic", "control_ftbl", "param", "norm2", "method", "sln", "nb_meas", "simcumom", "nb_fmn", "simfmn", "nb_poolm", "simpool", "measurements", "opt_wrapper", "dfl_dffg"))
            if (TIMEIT) {
               cat("cl optim: ", date(), "\\n", sep="", file=fclog)
            }
         }
         #mc_res=mclapply(1L:nmc, mc_sim)
         clusterSetRNGStream(cl)
         mc_res=parLapply(cl, 1L:nmc, mc_sim)
         stopCluster(cl)
      } else {
         mc_res=lapply(1L:nmc, mc_sim)
      }
      free_mc=sapply(mc_res, function(l) {if (class(l)=="character" || is.na(l$cost) || l$err) { ret=rep(NA, nb_param+3) } else { ret=c(l$cost, l$it, l$normp, l$par) }; ret })
      if (length(free_mc)==0) {
         cat("Parallel exectution of Monte-Carlo simulations has failed.", "\\n", sep="", file=fcerr)
         free_mc=matrix(NA, nb_param+2, 0)
      }
      cost_mc=free_mc[1,]
      nmc_real=nmc-sum(is.na(free_mc[4,]))
      cat("monte-carlo\\n", file=fkvh)
      cat("\\tsample\\n", file=fkvh)
      indent=2
      obj2kvh(nmc, "requested number", fkvh, indent)
      obj2kvh(nmc_real, "calculated number", fkvh, indent)
      obj2kvh(nmc-nmc_real, "failed to calculate", fkvh, indent)
      # convergence section in kvh
      indent=1
      mout=rbind(round(free_mc[1:2,,drop=F], 2),
         format(free_mc[3,,drop=F], di=2, sci=T))
      dimnames(mout)=list(c("cost", "it.numb", "normp"), iseq(ncol(free_mc)))
      obj2kvh(mout, "convergence per sample", fkvh, indent)
#expp      if (nb_poolf) {
#expp         # from log to natural concentraion
#expp         free_mc[nb_ff+nb_sc+1:nb_poolf,]=exp(free_mc[nb_ff+nb_sc+1:nb_poolf,])
#expp      }
      # remove failed m-c iterations
      free_mc=free_mc[-(1:3),,drop=F]
      ifa=which(is.na(free_mc[1,]))
      if (length(ifa)) {
         if (ncol(free_mc) > length(ifa)) {
            cat("Some Monte-Carlo iterations failed.", "\\n", sep="", file=fcerr)
         }
         free_mc=free_mc[,-ifa,drop=F]
         cost_mc=cost_mc[-ifa]
      }
      if (nmc_real <= 1) {
         cat("No sufficient calculated monter-carlo samples for statistics.", "\\n", sep="", file=fcerr)
         break
      }
#browser()
      rownames(free_mc)=nm_par
      indent=1
      obj2kvh(avaco, "detected cores", fkvh, indent)
      avaco=max(1, avaco, na.rm=T)
      obj2kvh(min(avaco, options()$mc.cores, na.rm=T), "used cores", fkvh, indent)
      
      # cost section in kvh
      cat("\\tcost\\n", file=fkvh)
      indent=2
      obj2kvh(mean(cost_mc), "mean", fkvh, indent)
      obj2kvh(median(cost_mc), "median", fkvh, indent)
      obj2kvh(sd(cost_mc), "sd", fkvh, indent)
      obj2kvh(sd(cost_mc)*100/mean(cost_mc), "rsd (%)", fkvh, indent)
      obj2kvh(quantile(cost_mc, c(0.025, 0.975)), "ci", fkvh, indent)
      
      # free parameters section in kvh
      cat("\\tStatistics\\n", file=fkvh)
      mout=c()
      indent=2
      # param stats
      # mean
      parmean=apply(free_mc, 1, mean)
      # median
      parmed=apply(free_mc, 1, median)
#browser()
      # covariance matrix
      covmc=cov(t(free_mc))
      obj2kvh(covmc, "covariance", fkvh, indent)
      # sd
      sdmc=sqrt(diag(covmc))
      # confidence intervals
      ci_mc=t(apply(free_mc, 1, quantile, probs=c(0.025, 0.975)))
      ci_mc=cBind(ci_mc, t(diff(t(ci_mc))))
      colnames(ci_mc)=c("CI 2.5%", "CI 97.5%", "CI length")
      mout=cBind(mout, mean=parmean, median=parmed, sd=sdmc,
         "rsd (%)"=sdmc*100/abs(parmean), ci_mc)
      obj2kvh(mout, "free parameters", fkvh, indent)

      # net-xch01 stats
      fallnx_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fallnx)
      fallnx=param2fl(param, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fallnx
      if (length(fallnx_mc)) {
         dimnames(fallnx_mc)[[1]]=nm_fallnx
         # form a matrix output
         fallout=matrix(0, nrow=nrow(fallnx_mc), ncol=0)
         # mean
#browser()
         parmean=apply(fallnx_mc, 1, mean)
         # median
         parmed=apply(fallnx_mc, 1, median)
         # covariance matrix
         covmc=cov(t(fallnx_mc))
         dimnames(covmc)=list(nm_fallnx, nm_fallnx)
         # sd
         sdmc=sqrt(diag(covmc))
         # confidence intervals
         ci_mc=t(apply(fallnx_mc, 1, quantile, probs=c(0.025, 0.975)))
         ci_mc=cBind(ci_mc, t(diff(t(ci_mc))))
         ci_mc=cBind(ci_mc, ci_mc[,3]*100/abs(parmean))
         colnames(ci_mc)=c("CI 2.5%", "CI 97.5%", "CI 95% length", "relative CI (%)")
         fallout=cBind(fallout, mean=parmean, median=parmed, sd=sdmc,
            "rsd (%)"=sdmc*100/abs(fallnx), ci_mc)
         o=order(nm_fallnx)
         obj2kvh(fallout[o,,drop=F], "all net-xch01 fluxes", fkvh, indent)
         obj2kvh(covmc[o,o], "covariance of all net-xch01 fluxes", fkvh, indent)

         # fwd-rev stats
         fwrv_mc=apply(free_mc, 2, function(p)param2fl(p, nb_f, nm_list, invAfl, p2bfl, g2bfl, bp, fc)$fwrv)
         dimnames(fwrv_mc)[[1]]=nm_fwrv
         fallout=matrix(0, nrow=nrow(fwrv_mc), ncol=0)
         # mean
         parmean=apply(fwrv_mc, 1, mean)
         # median
         parmed=apply(fwrv_mc, 1, median)
         # covariance matrix
         covmc=cov(t(fwrv_mc))
         dimnames(covmc)=list(nm_fwrv, nm_fwrv)
         # sd
         sdmc=sqrt(diag(covmc))
         # confidence intervals
         ci_mc=t(apply(fwrv_mc, 1, quantile, probs=c(0.025, 0.975)))
         ci_mc=cBind(ci_mc, t(diff(t(ci_mc))))
         ci_mc=cBind(ci_mc, ci_mc[,3]*100/abs(fwrv))
         dimnames(ci_mc)[[2]]=c("CI 2.5%", "CI 97.5%", "CI 95% length", "relative CI (%)")
         fallout=cBind(fallout, mean=parmean, median=parmed, sd=sdmc,
            "rsd (%)"=sdmc*100/abs(parmean), ci_mc)
         o=order(nm_fwrv)
         obj2kvh(fallout[o,,drop=F], "forward-reverse fluxes", fkvh, indent)
         obj2kvh(covmc[o,o], "covariance of forward-reverse fluxes", fkvh, indent)
      }
      break
   }
   if (length(sensitive) && nchar(sensitive) && sensitive != "mc") {
      cat(paste("Unknown sensitivity '", sensitive, "' method chosen.", sep=""), "\\n", sep="", file=fcerr)
   }

   if (TIMEIT) {
      cat("linstats: ", date(), "\\n", sep="", file=fclog)
   }
   # Linear method based on jacobian x_f
   # reset fluxes and jacobians according to param
   jx_f=jx_f_last
   if (is.null(jx_f$jacobian)) {
      rres=cumo_resid(param, cjac=T, jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spa, emu, pool, ipooled)
      jx_f=rres$jx_f
   } # else use last calculated jacobian

   # covariance matrix of free fluxes
   svj=svd(jx_f$jacobian)
   i=svj$d/svj$d[1]<1.e-10
   if (all(!i) && svj$d[1]<1.e-10) {
      # we could not find very small d, take just the last
      i[length(i)]=T
   }
   ibad=apply(svj$v[, i, drop=F], 2, which.contrib)
   ibad=unique(unlist(ibad))
   if (length(ibad) > 0) {
      cat(paste(if (nchar(runsuf)) runsuf%s+%": " else "", "Inverse of covariance matrix is numerically singular.\\nStatistically undefined parameter(s) seems to be:\\n",
         paste(nm_par[ibad], collapse="\\n"), "\\nFor more complete list, see sd columns in '/linear stats'\\nin the result file.", sep=""), "\\n", sep="", file=fcerr)
   }
   # "square root" of covariance matrix (to preserve numerical positive definitness)
   rtcov=(svj$u)%*%(t(svj$v)/svj$d)
   # standart deviations of free fluxes
   cat("linear stats\\n", file=fkvh)

   # sd free+dependent+growth net-xch01 fluxes
   nm_flfd=c(nm_ff, nm_fgr, nm_fl)
   if (nb_ff > 0 || nb_fgr > 0) {
      i=1:nb_param
      i=c(head(i, nb_ff), tail(i, nb_fgr))
      covfl=crossprod(rtcov[, i, drop=F]%mmt%(rBind(diag(nb_ff+nb_fgr), dfl_dffg)%mrv%c(rep.int(1., nb_ff), fgr)))
      dimnames(covfl)=list(nm_flfd, nm_flfd)
      sdfl=sqrt(diag(covfl))
   } else {
      sdfl=rep(0., nb_fl)
   }
   fl=c(head(param, nb_ff), fgr, flnx)
   mtmp=cBind("value"=fl, "sd"=sdfl, "rsd"=sdfl/abs(fl))
   rownames(mtmp)=nm_flfd
   o=order(nm_flfd)
   obj2kvh(mtmp[o,,drop=F], "net-xch01 fluxes (sorted by name)", fkvh, indent=1)
   obj2kvh(covfl[o, o], "covariance net-xch01 fluxes", fkvh, indent=1)

   # sd of all fwd-rev
   if (nb_ff > 0 || nb_fgr > 0) {
      i=1:nb_param
      i=c(head(i, nb_ff), tail(i, nb_fgr))
      covf=crossprod(tcrossprod(rtcov[,i, drop=F], jx_f$df_dffp%mrv%c(rep.int(1., nb_ff), head(poolall[nm_poolf], nb_fgr))))
      dimnames(covf)=list(nm_fwrv, nm_fwrv)
      sdf=sqrt(diag(covf))
   } else {
      sdf=rep(0., length(fwrv))
   }
   mtmp=cBind(fwrv, sdf, sdf/abs(fwrv))
   dimnames(mtmp)[[2]]=c("value", "sd", "rsd")
   o=order(nm_fwrv)
   obj2kvh(mtmp[o,], "fwd-rev fluxes (sorted by name)", fkvh, indent=1)
   if (nb_ff > 0 || nb_fgr > 0) {
      obj2kvh(covf, "covariance fwd-rev fluxes", fkvh, indent=1)
   }
   # pool -> kvh
   sdpf=poolall
   sdpf[]=0.

   if (nb_poolf > 0) {
      # covariance matrix of free pools
      # "square root" of covariance matrix (to preserve numerical positive definitness)
#expp      poolall[nm_poolf]=exp(param[nm_poolf])
      poolall[nm_poolf]=param[nm_poolf]
#expp      # cov wrt exp(pf)=pool
#expp      covpf=crossprod(rtcov[,nb_ff+nb_sc+1:nb_poolf, drop=F]%mrv%poolall[nm_poolf])
      # cov poolf
      covpf=crossprod(rtcov[,nb_ff+nb_sc+1:nb_poolf, drop=F])
      dimnames(covpf)=list(nm_poolf, nm_poolf)
      sdpf[nm_poolf]=sqrt(diag(covpf))
   }
   if (length(poolall) > 0) {
      mtmp=cBind("value"=poolall, "sd"=sdpf, "rsd"=sdpf/poolall)
      rownames(mtmp)=nm_poolall
      o=order(nm_poolall)
      obj2kvh(mtmp[o,,drop=F], "metabolite pools (sorted by name)", fkvh, indent=1)
      if (nb_poolf > 0) {
         o=order(nm_poolf)
         obj2kvh(covpf[o, o], "covariance free pools", fkvh, indent=1)
      }
   }

   # khi2 test for goodness of fit
   # goodness of fit (khi2 test)
   nvres=sum(!is.na(jx_f$res))
   if (nvres >= nb_param) {
      khi2test=list("khi2 value"=rcost, "data points"=nvres,
         "fitted parameters"=nb_param, "degrees of freedom"=nvres-nb_param)
      khi2test$`khi2 reduced value`=khi2test$`khi2 value`/khi2test$`degrees of freedom`
      khi2test$`p-value, i.e. P(X^2<=value)`=pchisq(khi2test$`khi2 value`, df=khi2test$`degrees of freedom`)
      khi2test$conclusion=if (khi2test$`p-value, i.e. P(X^2<=value)` > 0.95) "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD" else "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD"
      obj2kvh(khi2test, "goodness of fit (khi2 test)", fkvh, indent=1)
   } else {
      cat(sprintf("khi2: Measurment number %d is lower than parameter number %d. Khi2 test cannot be done.\\n", nvres, nb_param), sep="", file=fcerr)
   }
   if (prof) {
      Rprof(NULL)
   }
   close(fkvh)
""")
    f.write("""
   # write edge.netflux property
   fedge=file("%(d)s/edge.netflux.%(org)s" %%s+%% runsuf, "w")
   cat("netflux (class=Double)\\n", sep="", file=fedge)
   nm_edge=names(edge2fl)
   cat(paste(nm_edge, fallnx[edge2fl], sep=" = "), sep="\\n" , file=fedge)
   close(fedge)

   # write edge.xchflux property
   fedge=file("%(d)s/edge.xchflux.%(org)s" %%s+%% runsuf, "w")
   flxch=paste(".x", substring(edge2fl, 4), sep="")
   ifl=charmatch(flxch, substring(names(fallnx), 2))
   cat("xchflux (class=Double)\\n", sep="", file=fedge)
   cat(paste(nm_edge, fallnx[ifl], sep=" = "), sep="\\n" , file=fedge)
   close(fedge)

   # write node.log2pool property
   if (length(poolall)> 0) {
      fnode=file("%(d)s/node.log2pool.%(org)s" %%s+%% runsuf, "w")
      cat("log2pool (class=Double)\\n", sep="", file=fnode)
      nm_node=substring(names(poolall), 4)
      cat(paste(nm_node, log2(poolall), sep=" = "), sep="\\n" , file=fnode)
      close(fnode)
   }
}

#expp # back from log in pres
#expp i=grep("^pf:", nm_par)
#expp pres[i,]=exp(pres[i,])
pres=rBind(cost=costres, pres)
fco=file("%(d)s/%(org)s.pres.csv", open="w")
cat("row_col\t", file=fco)
write.table(file=fco, pres, row.n=T, quot=F, sep="\\t")
close(fco)
if (TIMEIT) {
   cat("rend    : ", date(), "\\n", sep="", file=fclog)
}
close(fclog)
close(fcerr)
"""%{
    "org": escape(org, "\\"),
    "d": escape(dirorg, "\\"),
})

    f.close()
    # try to make output files just readable to avoid later casual edition
    try:
	os.chmod(n_R, stat.S_IREAD)
    except:
        pass
