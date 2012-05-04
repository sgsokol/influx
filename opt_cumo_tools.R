#DEBUG=0; # 1 to enable debugging information, 0=disable
#TIMEIT=0; # 1 to enable time printing at some stages
if (length(find("TIMEIT")) && TIMEIT) {
   cat("load    : ", date(), "\n", sep="");
}
jx_f=list();
suppressPackageStartupMessages(library(bitops));
#suppressPackageStartupMessages(library(MASS)); # for generalized inverse
#suppressPackageStartupMessages(library(fUtilities)); # for Heaviside function
suppressPackageStartupMessages(library(nnls)); # for non negative least square
#suppressPackageStartupMessages(library(lattice)); # to keep Matrix silent
suppressPackageStartupMessages(library(Matrix, warn=F, verbose=F)); # for sparse matrices
suppressPackageStartupMessages(library(expm, warn=F, verbose=F)); # for sparse matrices
#library(inline); # for inline fortran compilation

trisparse_solv=function(A, b, w, method="dense") {
   # solve A*x=b where A=tridiag(Al,Ac,Au)+s*e^t and b is dense
   if (method=="dense") {
      n=ncol(A);
      if (DEBUG) {
         write.matrix(A, file=paste("dbg_Acumo_d_",w,".txt", sep=""),sep="\t");
      }
      # factorize the matrix
      fA=qr(A);
      d=diag(fA$qr)
      fA$rank=sum(abs(d)>abs(d[1])*1.e-14)
      if (fA$rank < n) {
          mes=paste("trisparse_solv: Cumomer matrix ",
             n, "x", n, " at weight ", w,
             " is singular.\nUnsolvable cumomers are:\n",
             paste(fA$pivot[-(1:fA$rank)], sep="", collapse="\n"),
             "\nThe matrix is dumped in dbg_Acumo_singular.txt\n",
             sep="", collapse="");
          write.matrix(A, file="dbg_Acumo_singular.txt");
          return(list(x=NULL, fA=fA, err=1, mes=mes));
      }
      x=solve(fA,b);
      return(list(x=x, fA=fA, err=0, mes=NULL));
   } else if (method=="sparse") {
      # sparse
      # fulfill a matrix
      #require(Matrix);
      #A=Matrix(A, sparse=T);
#browser();
      x=try(solve(A,b)); # A has its factorized form
      if (inherits(x, "try-error")) {
         # find 0 rows if any
         izc=apply(A, 1, function(v)sum(abs(v))<=1.e-10)
         izf=names(which(abs(jx_f$fwrv)<1.e-7))
         if (sum(izc) || length(izf)) {
            mes=paste("Cumomer matrix is singular. Try '--clownr N' option with small N, say 1.e-3\nor constrain the flux(es) (see below) to be non zero\n",
               "Zero rows in cumomer matrix A at weight ", w, ":\n",
               paste(rownames(A)[izc], collapse="\n"), "\n",
               "Zero fluxes are:\n",
               paste(izf, collapse="\n"), "\n",
               sep="")
         } else {
            mes="Cumomer matrix is singular. Try '--clownr N' option with small N, say 1.e-3."
         }
#browser()
         stop(mes)
      }
      if (DEBUG) {
         cat("A=", str(A), "\n", sep="");
      }
      return(list(x=as.numeric(x), fA=A, err=0, mes=NULL));
   } else if (method=="smw") {
      # Sherman-Morrison-Woodbury for low rank matrix modification
      require(matrid, lib.loc="/home/sokol/R/lib");
      atrim=new("matridm", A);
      if (DEBUG) {
         cat(paste("dim A at weight ", w, ":\n", sep=""));
         print(dim(A));
         write.matrix(cbind(A,b=b),file=paste("dbg_tridmA_",w,".txt", sep=""),sep="\t");
      #   print(A);
      }
      x=qr.solve(atrim,b);
      return(x);
   } else {
      stop(paste("Unknown method '", method, "'", sep=""));
   }
}

dfc2fallnx=function(nb_f, flnx, param, fc) {
   # produce complete flux (net,xch)*(dep,free,constr) vector
   # from dep,free,constr
   f=numeric(0);
   if (nb_f$nb_fln) {
      f=c(f, flnx[1:nb_f$nb_fln]);
   }
   if (nb_f$nb_ffn) {
      f=c(f, param[1:nb_f$nb_ffn]);
   }
   if (nb_f$nb_fcn) {
      f=c(f, fc[1:nb_f$nb_fcn]);
   }
   if (nb_f$nb_flx) {
      f=c(f, flnx[(nb_f$nb_fln+1):nb_f$nb_fl]);
   }
   if (nb_f$nb_ffx) {
      f=c(f, param[(nb_f$nb_ffn+1):nb_f$nb_ff]);
   }
   if (nb_f$nb_fcx) {
      f=c(f, fc[(nb_f$nb_fcn+1):nb_f$nb_fc]);
   }
   return(f);
}

cumo_resid=function(param, cjac=TRUE, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb) {
   if (length(measinvvar)) {
      sqm=sqrt(measinvvar);
   } else {
      sqm=c();
   }
   if (length(invfmnvar)) {
      sqf=sqrt(invfmnvar);
   } else {
      sqf=c();
   }
   # find x for all weights
   if (!identical(param, jx_f$param) ||
         (cjac && is.null(jx_f$x_f))) {
      lres=param2fl_x(param, cjac, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAb);
      if (!is.null(lres$err) && lres$err) {
         return(list(err=1, mes=lres$mes));
      }
      # store x, fluxes, ... in global space
      jx_f$x<<-lres$x;
      jx_f$fallnx<<-lres$fallnx;

      # find simulated scaled measure vector scale*(measmat*x)
      jx_f$usimcumom<<-(measmat%*%c(jx_f$x[imeas],1.))[,1];
      simvec=jx_f$usimcumom*c(1.,param)[ir2isc];

      # diff between simulated and measured
      res=c((simvec-measvec), (jx_f$fallnx[ifmn]-fmn));
      jx_f$ures<<-res;
      jx_f$res<<-res*c(sqm, sqf);
      # invalidate old jacobian as x_f was recalculated
      jx_f$jacobian<<-NULL
   } # else we have everything in jx_f
   # jacobian
   if (cjac) {
      if (!identical(param, jx_f$param) || is.null(jx_f$jacobian)) {
         # recalculate it
         cumo_jacob(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl,
            bp, fc, xi, imeas, measmat, measvec, ir2isc);
         jacobian=jx_f$udr_dp*c(sqm,sqf);
         jx_f$jacobian <<- jacobian;
      } else {
         jacobian=jx_f$jacobian
      }
   } else {
      jacobian=NULL;
   }

   return(list(res=jx_f$res, fallnx=jx_f$fallnx,
      jacobian=jacobian));
}
icumo_resid=function(param, cjac=TRUE, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvecti, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb, pool, ti, tifull) {
   # claculates residual vector of labeling propagation corresponding to param
   #cat("icumo_resid: param=", param, ", cjac=", cjac, "\n")
   if (length(measinvvar)) {
      sqm=sqrt(measinvvar);
   } else {
      sqm=c();
   }
   if (length(invfmnvar)) {
      sqf=sqrt(invfmnvar);
   } else {
      sqf=c();
   }
   nb_ti=length(ti)
   nb_sc=nb_f$nb_sc
   nb_poolf=nb_f$nb_poolf
   
   jacobian=NULL;
   # find usimvec
   recalcjac=(cjac && !identical(param, jx_f$param)) ||
         (cjac && (is.null(jx_f$ujaclab) || is.null(jx_f$uujac)))
   recalcx=!identical(param, jx_f$param) ||
         (cjac && is.null(jx_f$ujaclab))
   if (recalcx) {
      lres=param2fl_usm(param, cjac, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAb, pool, ti, measmat, imeas, tifull);
      if (!is.null(lres$err) && lres$err) {
         return(list(err=1, mes=lres$mes));
      }
      # scale simulated measurements scale*(usm)
      if (nb_sc) {
         simvec=jx_f$usm*c(1.,param)[ir2isc];
      } else {
         simvec=jx_f$usm
      }
      jx_f$simvec <<- simvec
      # diff between simulated and measured
      if (is.null(measvecti)) {
         jx_f$res <<- NULL
      } else {
         jx_f$ureslab <<- simvec-measvecti
         jx_f$uresflu <<- jx_f$fallnx[nm$fmn]-fmn
         jx_f$reslab <<- jx_f$ureslab*rep(sqm, nb_ti-1)
         jx_f$resflu <<- jx_f$uresflu*sqf
         jx_f$res <<- c(jx_f$reslab, jx_f$resflu);
         names(jx_f$res) <<- c(outer(rownames(jx_f$reslab), paste("ti=", colnames(jx_f$reslab), sep=""), paste), names(jx_f$resflu))
         jx_f$ures <<- c(jx_f$ureslab, jx_f$uresflu);
         names(jx_f$ures) <<- names(jx_f$res)
      }
   }
   if (recalcjac) {
      # add measured fluxes part of jacobian
      if (nb_f$nb_fmn) {
         mdfm_dff=dfm_dff()
         jx_f$uujac <<- rbind(jx_f$ujaclab, cbind(mdfm_dff, matrix(0., nrow=nrow(mdfm_dff), ncol=nb_sc+nb_poolf)))
      }
      # scale jacobian by sd
      jacobian=c(rep(sqm, nb_ti-1), sqf)*jx_f$uujac
      #colnames(jacobian)=nm$par
      if (nb_ff > 0) {
         jx_f$dr_dff <<- jx_f$uujac[,1:nb_ff,drop=F]
      } else {
         jx_f$dr_dff <<- jx_f$uujac[,0,drop=F]
      }
      jx_f$udr_dp <<- jx_f$uujac
      jx_f$jacobian <<- jacobian
   } # else we have everything in jx_f
   #if (cjac) {
   #   cat("icumo_resid: jx_f$param=", jx_f$param, "\n")
   #}
   return(list(res=jx_f$res, jacobian=jx_f$jacobian, ures=jx_f$ures, usm=jx_f$usm, simvec=jx_f$simvec, fallnx=jx_f$fallnx));
}
cumo_cost=function(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb) {
   resl=cumo_resid(param, cjac=FALSE, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
   if (!is.null(resl$err) && resl$err) {
      return(NULL);
   }
   res=resl$res;
   fn=sum(res*res);
   if (DEBUG) {
      write.matrix(fn, file="dbg_cost.txt", sep="\t");
   }
   return(fn);
}
icumo_cost=function(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb, pool, ti, tifull) {
   resl=icumo_resid(param, cjac=FALSE, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb, pool, ti, tifull);
   if (!is.null(resl$err) && resl$err) {
      return(NULL);
   }
   res=resl$res;
   fn=sum(res*res);
   if (DEBUG) {
      write.matrix(fn, file="dbg_cost.txt", sep="\t");
   }
   return(fn);
}
cumo_grad=function(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb) {
   # calculate gradient of cost function for cumomer minimization problem
   # method: forward finite differences f(x+h)-f(x)/h
   # x+h is taken as (1+fact)*x
   fact=1.e-7;
   grad=param; # make place for gradient
   # f(x)
   f=cumo_cost(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
   for (i in 1:length(param)) {
      x=param[i];
      h=x*fact;
      param[i]=x+h;
      if (param[i]==x) {
         # we are too close to zero here
         param[i]=fact;
      }
      fh=cumo_cost(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
      # restore modified param
      param[i]=x;
      grad[i]=(fh-f)/h;
   }
   return(grad);
}

param2fl=function(param, nb_f, nm, invAfl, p2bfl, bp, fc) {
   # claculate all fluxes from free fluxes
   flnx=c(invAfl%*%(p2bfl%*%head(param, nb_f$nb_ff)+c(bp)));
   names(flnx)=nm$flnx;
   fallnx=c(dfc2fallnx(nb_f, flnx, param, fc));
   names(fallnx)=nm$fallnx;
   fwrv=c(fallnx2fwrv(fallnx));
   names(fwrv)=nm$fwrv;
   if (DEBUG) {
      write.matrix(p2bfl%*%head(param, nb_f$nb_ff)+bp, file="dbg_bfl.txt", sep="\t");
      n=length(fwrv);
      names(fwrv)=nm_fwrv;
      write.matrix(fwrv, file="dbg_fwrv.txt", sep="\t");
      write.matrix(cbind(1:n,nm_fallnx,fallnx), file="dbg_fallnx.txt", sep="\t");
   }
   return(list(fallnx=fallnx, fwrv=fwrv, flnx=flnx));
}

p2f=function(param, nb_f, nm, invAfl, p2bfl, bp, fc, xi) {
   # translate param to fwd-rev fluxes (for num deriv purpose only)
   return(param2fl(param, nb_f, nm, invAfl, p2bfl, bp, fc, xi)$fwrv);
}

param2fl_x=function(param, cjac=TRUE, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAb) {
   # translate free params (fluxes+scales) to fluxes and cumomers
   nbw=length(spAb);
   if (!is.null(jx_f$param) && identical(param, jx_f$param) &&
      (length(jx_f$x)==sum(spAb[[nbw]]$tb_x))) {
      if (cjac) {
          if (!is.null(jx_f$x_f)) {
             return(list(x=jx_f$x, x_f=jx_f$x_f, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx));
          } # else recalculate it here
      } else {
         # just x, fallnx, ... that are already calculated
         return(list(x=jx_f$x, x_f=jx_f$x_f, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx));
      }
   }
   # cumulated sum
   nbc_cumos=c(0, cumsum(nb_cumos))
   # calculate all fluxes from free fluxes
   lf=param2fl(param, nb_f, nm, invAfl, p2bfl, bp, fc);
   jx_f$fallnx<<-lf$fallnx;
   jx_f$fwrv<<-lf$fwrv
   jx_f$flnx<<-lf$flnx
   # construct the system A*x=b from fluxes
   # and find x for every weight
   # if fj_rhs is not NULL, calculate jacobian x_f
   nb_fwrv=length(lf$fwrv);
   nb_xi=length(xi);
   x=c(1, xi);
   if (cjac) {
      x_f=matrix(0., nrow=sum(nb_cumos), ncol=nb_fwrv);
   } else {
      x_f=NULL;
   }
   if (DEBUG) {
      tmp=lf$fwrv;
      names(tmp)=nm$fwrv;
      conct=file("dbg_fwrv.txt", "wb");
      obj2kvh(tmp, "fwrv", conct);
      tmp=lf$fallnx;
      names(tmp)=nm$fallnx
      obj2kvh(tmp, "net-xch", conct);
   }
#browser();
   ba_x=0;
   for (iw in 1:nb_w) {
      nx=length(x);
      nxl=nx-1-nb_xi; # number of lighter cumomers
      #ncumow=nb_cumos[iw];
      #A=matrix(0.,ncumow,ncumow);
      #b=double(ncumow);
      #res<-.Fortran(fortfun, fl=as.double(lf$fwrv), nf=nb_fwrv, x=as.double(x[]), iw=as.integer(iw), n=as.integer(ncumow), A=as.matrix(A), b=as.double(b), calcA=as.integer(TRUE), calcb=as.integer(TRUE), NAOK=TRUE, DUP=FALSE);
      # old usage of cumo.f (for tests only)
      #lAb=fwrv2Ab(lf$fwrv, spAb[[iw]], x, nm$rcumo[(nbc_cumos[iw]+1):nbc_cumos[iw+1]]);
      lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], x, nm$rcumo[(nbc_cumos[iw]+1):nbc_cumos[iw+1]]);
      A=lAb$A;
      b=lAb$b;
      #if (any(A != lAbr$A) || any(b != lAbr$b)) {
      #   browser()
      #}
      if (DEBUG) {
         write.matrix(cbind(as.matrix(A), b=b), file=paste("dbg_cumoAb_",iw,".txt", sep=""), sep="\t");
         if (iw==1) {
            jx_f$wA<<-list(A);
            jx_f$wb<<-list(b);
         } else {
            jx_f$wA<<-append(jx_f$wA, list(A));
            jx_f$wb<<-append(jx_f$wb, list(b));
         }
      }
      #solve the system A*x=b;
      #lsolv=trisparse_solv(A, b, iw, method=ifelse(ncumow>200, "sparse", "dense"));
      lsolv=trisparse_solv(lAb$A, lAb$b, iw, method="sparse");
      #lsolv=lsolvb;
      if (!is.null(lsolv$err) && lsolv$err) {
         return(list(err=1, mes=lsolv$mes));
      }
      #stopifnot(isTRUE(all.equal(lsolv$x, lsolvb$x)));
      xw=lsolv$x;
      nxw=length(xw);
      if (cjac) {
         # calculate jacobian x_f
         # first, calculate right hand side for jacobian calculation
         #j_rhsw=matrix(0., nxw, nb_fwrv);
         #b_x=matrix(0., nxw, nxl);
         #res<-.Fortran(fj_rhs, fl=as.double(lf$fwrv), nf=nb_fwrv, x=as.double(x), xw=as.double(xw), x_f=as.double(x_f), nx=as.integer(length(x)), nxw=as.integer(nxw), iw=as.integer(iw), j_rhs=as.matrix(j_rhsw), b_x=as.matrix(b_x), NAOK=TRUE, DUP=FALSE);
         # j_rhsw, b_x from sparse matrices
         # bind cumomer vector
#browser();
         x=c(x,xw);
         j_b_x=fx2jr(lf$fwrv, spAb[[iw]], x);
         j_rhsw=j_b_x$j_rhsw;
         b_x=j_b_x$b_x;
         if (DEBUG) {
            write.matrix(cbind(j_rhsw, b_x),
               file=paste("dbg_j_rhs_",iw,".txt", sep=""), sep="\t");
            if (iw==1) {
               jx_f$wbx<<-list(b_x);
               jx_f$wjr<<-list(j_rhsw);
            } else {
               jx_f$wbx<<-append(jx_f$wbx, list(b_x));
               jx_f$wjr<<-append(jx_f$wjr, list(j_rhsw));
            }
         }
#browser()
         if (iw > 1) {
            x_f[ba_x+(1:nb_cumos[iw]),]=
               as.matrix(solve(lsolv$fA, j_rhsw+b_x%*%x_f[1:ba_x,]));
         } else {
            x_f[ba_x+(1:nb_cumos[iw]),]=
               as.matrix(solve(lsolv$fA, j_rhsw));
         }
      } else {
         # bind cumomer vector
         x=c(x,xw);
         if (!is.null(jx_f$param) && !identical(param, jx_f$param)) {
            x_f=NULL;
         }
      }
      ba_x=ba_x+nb_cumos[iw];
   }
#print(x);
   # store usefull information in global list jx_f
   jx_f$param<<-param;
   jx_f$x_f<<-x_f;
   return(append(list(x=x[(2+nb_xi):length(x)], x_f=x_f), lf));
}

param2fl_usm=function(param, cjac=TRUE, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, spAb, pool, ti, measmat, irmeas, tifull=ti) {
   # translate free params (fluxes+scales) to fluxes and
   # unscaled simulated measurements (usm) for labeling propagation
   # tifull may be more fine grained than ti. All ti must be in tifull
   # only ti moments are reported in usm and jacobian
   
   #print(jx_f$param-param)
   
   # recalculate or not the labeling?
   calcx=is.null(jx_f$param) ||
      !identical(param, jx_f$param) ||
      (length(jx_f$x)!=dim(spAb[[nb_w]]$tb_x)[1]+nb_cumos[nb_w])
   if (!calcx) {
      if (cjac) {
         if (!is.null(jx_f$ujaclab)) {
            return(list(usm=jx_f$usm, x=jx_f$x, ujaclab=jx_f$ujaclab, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx));
         } # else recalculate it here
      } else {
         # just x, fallnx, ... that are already calculated
         return(list(usm=jx_f$usm, x=jx_f$x, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv));
      }
   }
   # here calcx==T or (calcx==F && cjac==T && is.null(jx_f$ujaclab))
   # so recalculate what is not in the cache in jx_f
   nb_ti=length(ti)
   nb_tifu=length(tifull)
   if (nb_ti < 2) {
      return(list(err=1, mes="Number of time points is less than 2"))
   }
   if (!all(ti %in% tifull)) {
      return(list(err=1, mes="Not all time moments in ti are present in tifull vector"))
   }
   if (calcx) {
      cat("param2fl_usm: recalculate labprop\n")
   }
   
   dt=diff(tifull)
   # cumulated sum
   nbc_cumos=c(0, cumsum(nb_cumos))
   # calculate all fluxes from free fluxes
   lf=param2fl(param, nb_f, nm, invAfl, p2bfl, bp, fc);
   jx_f$fallnx <<- lf$fallnx
   jx_f$fwrv <<- lf$fwrv
   jx_f$flnx <<- lf$flnx
   nb_fwrv=length(lf$fwrv)
   nb_xi=length(xi);
   nb_poolf=nb_f$nb_poolf
   nb_sc=nb_f$nb_sc
   vsc=c(1.,param)[ir2isc]
   if (!is.null(measmat)) {
      nb_mcol=ncol(measmat)
      mema1=measmat[,-nb_mcol,drop=F]
      memaone=measmat[,nb_mcol]
      irmeas_xi=irmeas+nb_xi+1
   }
   # fullfill pool with free pools
   if (nb_poolf > 0) {
      pool[nm$poolf]=param[nm$poolf]
   }
   # prepare inverse of pool vectors
   invpool=1./pool
   # invm has the same length as full cumomer vector
   invm=invpool[nb_f$ip2ix]
   invmw=lapply(1:nb_w, function(iw)invm[nbc_cumos[iw]+(1:nb_cumos[iw])])
#browser()
   # prepare vectors at t1=0 with zero labeling
   if (calcx) {
      # incu, xi is supposed to be in [0; 1]
      x1=c(1., xi, rep(0., nbc_cumos[nb_w+1]))
      # unscaled simulated measurements
      usm=matrix(0., nrow=nrow(measmat), ncol=nb_ti-1)
      # construct the matrices invm*A in the systems pool*dx_dt=A*x+s from fluxes
      lwA=lapply(1:nb_w, function(iw) {fwrv2Abr(lf$fwrv, spAb[[iw]], x1, nm$rcumo[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], getb=F)$A*invmw[[iw]]})
      lwinva=lapply(lwA, Matrix::solve)
      #lwinva=lapply(lwA, function(a) {
      #   qa=qr(as.matrix(a), LAPACK=T)
      #   return(qr.solve(qa, tol=1e-14))
      #})
      # s
      s1=lapply(1:nb_w, function(iw) {as.double(fwrv2sp(lf$fwrv, spAb[[iw]], x1)$s)*invmw[[iw]]})
      # xp first derivative of x
      xp1=c(rep(0., nb_xi+1), unlist(s1)) # A*x1 is omited as x1==0
      # sp first derivative of s
      sp1=lapply(1:nb_w, function(iw) as.double(fwrv2sp(lf$fwrv, spAb[[iw]], x1, xp1, gets=F)$sp)*invmw[[iw]])
      expadt=list()
      xsim=matrix(0., nrow=length(x1), ncol=nb_tifu)
      xpsim=matrix(0., nrow=length(xp1), ncol=nb_tifu)
      xsim[,1]=x1
      xpsim[,1]=xp1
      spx=list(sp1)
   } else {
      xsim=jx_f$xsim
      xpsim=jx_f$xpsim
      expadt=jx_f$expadt
      lwA=jx_f$lwA
      lwinva=jx_f$lwinva
      spx=jx_f$spx
      usm=jx_f$usm
      
      x1=xsim[,1]
      xp1=xpsim[,1]
      sp1=spx[[1]]
   }
   if (cjac) {
      cat("param2fl_usm: recalculate jacobian\n")
      mdf_dff=df_dff(param, lf$flnx)
      jx_f$df_dff <<- mdf_dff
      nb_ff=ncol(mdf_dff)
      nb_poolf=nb_f$nb_poolf
      nb_meas=nrow(measmat)
      xff1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      xffp1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      if (nb_ff > 0) {
         # jacobian ff rhs
         sj1=lapply(1:nb_w, function(iw) {
            sj=fx2jr(lf$fwrv, spAb[[iw]], x1, xp1)
            sj$jrhs=(-invmw[[iw]])*(sj$j_rhsw%*%mdf_dff)
            sj$jrhsp=(-invmw[[iw]])*(sj$j_rhswp%*%mdf_dff)
            xffp1[nbc_cumos[iw]+1:nb_cumos[iw],] <<- as.matrix(sj$jrhs)
            return(sj)
         })
      }
      xpf1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
      xpfp1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
      if (nb_poolf > 0) {
         # jacobian pf rhs
         spf1=lapply(1:nb_w, function(iw) {
            spf=matrix(0., nrow=nb_cumos[iw], ncol=nb_poolf)
            i2x=nb_f$iparpf2ix[[iw]]
            icw=nbc_cumos[iw]+1:nb_cumos[iw]
            ixw=1+nb_xi+icw
            xp=xp1[ixw]
            spf[i2x]=(-invmw[[iw]]*xp)[i2x[,1]]
            xpfp1[icw,] <<- spf
            spfp=matrix(0., nrow=nb_cumos[iw], ncol=nb_poolf)
            xpp=lwA[[iw]]%*%xp+sp1[[iw]]
            spfp[i2x]=(-invmw[[iw]]*xpp)[i2x[,1]]
            return(list(s=spf, sp=spfp))
         })
      }
      dur_dsc=matrix(0., nrow=nb_meas, ncol=nb_sc)
      jacobian=array(0., dim=c(nb_meas, nb_ff+nb_sc+nb_poolf, 0))
   } else {
      jacobian=jx_f$ujaclab
   }

   s2=list()
   sp2=list()
   sj2=list()
   spf2=list()
   spfp2=list()
   if (cjac) {
      xff2=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      xffp2=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      xpf2=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
      xpfp2=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
   } else {
      xff2=NULL
      xpf2=NULL
   }
   # time run from time t1 to time t2
   for (iti in 2:nb_tifu) {
      dtr=as.character(round(dt[iti-1], 10))
      if (calcx) {
         x2=c(1., xi)
         xp2=rep(0., nb_xi+1)
      } else {
         x2=xsim[,iti]
         xp2=xpsim[,iti]
         sp2=spx[[iti]]
      }
      # prepare expm(a*dt)
      if (calcx) {
         lapply(1:nb_w, function(iw) {
            if (is.null(expadt[[dtr]])) {
               expadt[[dtr]] <<- list()
            }
            if (length(expadt[[dtr]]) < iw) {
               expadt[[dtr]][[iw]] <<- expm(as.matrix(lwA[[iw]]*dt[iti-1]))
            }
         })
      }
      for (iw in 1:nb_w) {
         # prepare xw1 and xwp1
         icw=nbc_cumos[iw]+(1:nb_cumos[iw])
         ixw=1+nb_xi+icw
         xw1=x1[ixw]
         xwp1=xp1[ixw]
         # prepare s and sp
         if (calcx) {
            sw1=s1[[iw]]
            swp1=sp1[[iw]]
            if (iw > 1) {
               ssp=fwrv2sp(lf$fwrv, spAb[[iw]], x2, xp2)
               sw2=as.double(ssp$s)*invmw[[iw]]
               swp2=as.double(ssp$sp)*invmw[[iw]]
            } else {
               # for the first weight the source is constant in pulse experiment
               sw2=sw1
               swp2=swp1
            }
            s2[[iw]]=sw2
            sp2[[iw]]=swp2

            # make a time step for xw
            if (T && iti > 2) {
               # in the hope that the s curve is smooth enough
               xw2=as.double(expm_cub_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sw1, sw2, swp1, swp2, xw1))
            } else {
               # for a steep step, low order to preserve monotony
               xw2=as.double(expm_lin_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sw1, sw2, xw1))
            }
            # bring to interval [0; 1]
            #xw2[xw2<0.]=0.
            #xw2[xw2>1.]=1.
            xwp2=as.double(lwA[[iw]]%*%xw2)+sw2
            x2=c(x2, xw2)
            xp2=c(xp2, xwp2)
         } else {
            xw2=x2[ixw]
            xwp2=xp2[ixw]
         }
         if (cjac) {
            # prepare source 2 for jacobian ff
            if (nb_ff > 0) {
               sj=fx2jr(lf$fwrv, spAb[[iw]], x2, xp2)
               sj$jrhs=(-invmw[[iw]])*(sj$j_rhsw%*%mdf_dff)
               sj$jrhsp=(-invmw[[iw]])*(sj$j_rhswp%*%mdf_dff)
               if (iw > 1) {
                  # add lighter cumomers to jacobian source term
                  sj$jrhs=sj$jrhs-invmw[[iw]]*(sj$b_x%*%xff2[1:nbc_cumos[iw],])
                  sj$jrhsp=sj$jrhsp-invmw[[iw]]*(sj$b_xp%*%xff2[1:nbc_cumos[iw],]+sj$b_x%*%xffp2[1:nbc_cumos[iw],])
               }
               sj2[[iw]]=sj
               # make a time step for jacobian xff
               if (T && iti > 2) {
                  xff2[icw,]=as.matrix(expm_cub_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sj1[[iw]]$jrhs, sj2[[iw]]$jrhs, sj1[[iw]]$jrhsp, sj2[[iw]]$jrhsp, xff1[icw,,drop=F]))
               } else {
                  xff2[icw,]=as.matrix(expm_lin_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sj1[[iw]]$jrhs, sj2[[iw]]$jrhs, xff1[icw,,drop=F]))
               }
               xffp2[icw,]=as.matrix(lwA[[iw]]%*%xff2[icw,,drop=F])+as.matrix(sj2[[iw]]$jrhs)
            }
            # prepare source 2 for jacobian poolf
            if (nb_poolf > 0) {
               if (nb_ff==0 && iw > 1) {
                  sj=fx2jr(lf$fwrv, spAb[[iw]], x2, xp2)
               }
               spf=matrix(0., nrow=nb_cumos[iw], ncol=nb_poolf)
               i2x=nb_f$iparpf2ix[[iw]]
               spf[i2x]=(-invmw[[iw]]*xwp2)[i2x[,1]]
               spfp=matrix(0., nrow=nb_cumos[iw], ncol=nb_poolf)
               xpp=lwA[[iw]]%*%xwp2+sp2[[iw]]
               spfp[i2x]=(-invmw[[iw]]*xpp)[i2x[,1]]
               if (iw > 1) {
                  # add lighter cumomers to jacobian source term
                  spf=spf-invmw[[iw]]*(sj$b_x%*%xpf2[1:nbc_cumos[iw],,drop=F])
                  spfp=spfp-invmw[[iw]]*(sj$b_xp%*%xpf2[1:nbc_cumos[iw],,drop=F]+sj$b_x%*%xpfp2[1:nbc_cumos[iw],,drop=F])
               }
               spf2[[iw]]=list(s=spf, sp=spfp)
               # make a time step for jacobian xpf
               if (iti > 2) {
                  xpf2[icw,]=as.matrix(expm_cub_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], spf1[[iw]]$s, spf2[[iw]]$s, spf1[[iw]]$sp, spf2[[iw]]$sp, xpf1[icw,,drop=F]))
               } else {
                  xpf2[icw,]=as.matrix(expm_lin_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], spf1[[iw]]$s, spf2[[iw]]$s, xpf1[icw,,drop=F]))
               }
               xpfp2[icw,]=as.matrix(lwA[[iw]]%*%xpf2[icw,,drop=F])+as.matrix(spf)
            }
         }
      }
      it=which(tifull[iti] == ti)
      if (calcx) {
         xsim[,iti]=x2
         xpsim[,iti]=xp2
         spx=append(spx, list(sp2))
         
         if (length(it) > 0) {
            m2=mema1%*%x2[irmeas_xi]+memaone
            usm[,it-1]=m2
         }
      } else {
         if (length(it) > 0) {
            m2=usm[,it-1]
         }
      }
      if (cjac && length(it) > 0) {
         # scale part of jacobian
         dur_dsc[]=0.
         if (nb_f$nb_sc > 0) {
            is2m=nb_f$is2m
            dur_dsc[is2m]=m2[is2m[,1]]
         }
         # append ff part of jacobian
         if (nb_ff > 0) {
            m2=mema1%*%xff2[irmeas,,drop=F]
            if (nb_sc > 0) {
               m2=vsc*m2
            }
         } else {
            m2=matrix(0., nrow=nb_meas, ncol=0)
         }
         if (nb_f$nb_poolf > 0) {
            mpf=mema1%*%xpf2[irmeas,,drop=F]
            if (nb_sc > 0) {
               mpf=vsc*mpf
            }
         } else {
            mpf=matrix(0., nrow=nb_meas, ncol=0)
         }
         jacobian=array(c(jacobian, m2, dur_dsc, mpf), dim=dim(jacobian)+c(0,0,1))
      }
      # prepare next step if any
      if (iti < nb_tifu) {
         s1=s2
         sp1=sp2
         x1=x2
         xp1=xp2
         if (cjac) {
            sj1=sj2
            xff1=xff2
            xffp1=xffp2
            spf1=spf2
            xpf1=xpf2
            xpfp1=xpfp2
         }
      }
   }
   dimnames(usm)=list(rownames(measmat), ti[-1])
   dimnames(xsim)=list(c("one", nm$xi, nm$rcumo), tifull)
   dimnames(xpsim)=list(c("one", nm$xi, nm$rcumo), tifull)
#print(x);
   # store usefull information in global list jx_f
   jx_f$param<<-param;
   jx_f$usm<<-usm;
   jx_f$xsim<<-xsim;
   jx_f$xpsim<<-xpsim;
   jx_f$spx<<-spx;
   jx_f$expadt<<-expadt;
   jx_f$lwA<<-lwA;
   jx_f$lwinva<<-lwinva;
   if (cjac) {
      # unscaled and permuted jacobian
      # index run measures at time 1 then 2, ... for first free flux
      # then the same for second free flux and so on
      jacobian=matrix(aperm(jacobian, c(1,3,2)), ncol=length(param))
      rownames(jacobian)=outer(nm$meas, ti[-1], paste, sep=", ti=")
      colnames(jacobian)=names(param)
      jx_f$ujaclab <<- jacobian
      jx_f$uujac <<- NULL # invalidate old ujac
   } else if (calcx) {
      # invalidate old jacobian as x were recalculated
      jx_f$ujaclab <<- NULL
   }
   x=x2[-(1:(nb_xi+1))]
   names(x)=nm$rcumo
   jx_f$x <<- x
   jx_f$xff <<- xff2
   jx_f$xpf <<- xpf2
   return(append(list(usm=usm, x=x, xff=xff2, xpf=xpf2, ujaclab=jx_f$ujaclab), lf));
}
num_jacob=function(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo") {
   # numerical calculation of jacobian dx_df
   # we variate fvrw one by one and recalculate the whole x
   # for each variation. The result is returned as matix.
   
   # calculate all fluxes from free fluxes
   lf=param2fl(param, nb_f, nm, invAfl, p2bfl, bp, fc);
   f0=lf$fwrv;
   dfl=0.00001;
   x_f=matrix(0., ncol=0, nrow=sum(nb_cumos));
   nb_fwrv=length(f0);
   for (i in 0:length(f0)) {
      f1=f0;
      if (i > 0) {
         f1[i]=f0[i]+dfl;
      }
      # take the system A*x=b from the global variables jx_f$wA jx_f$wb
      # and find x for every weight
      # if fj_rhs is not NULL, calculate jacobian x_f
      x=c(1, xi);
      for (iw in 1:nb_w) {
         ncumow=nb_cumos[iw];
         A=matrix(0.,ncumow,ncumow);
         b=double(ncumow);
         #fwrv2Abcumo(fl, nf, x, nx, iw, n, A, b)
         res<-.Fortran(fortfun,
            fl=as.double(f1),
            nf=nb_fwrv,
            x=as.double(x),
            iw=as.integer(iw),
            n=as.integer(ncumow),
            A=as.matrix(A),
            b=as.double(b),
            calcA=as.integer(TRUE),
            calcb=as.integer(TRUE),
            NAOK=TRUE,
            DUP=FALSE);
         # solve the system A*x=b;
         lsolv=trisparse_solv(A, b, iw, method=ifelse(ncumow>200, "sparse", "dense"));
         xw=lsolv$x;
         nxw=length(xw);
         # bind vectors and matrices
         x=c(x,xw);
      }
      if (i > 0) {
         x_f=cbind(x_f, (x-x0)/dfl);
      } else {
         x0=x;
         #x_f=cbind(x_f, x0);
      }
   }
   return(x_f);
}
num_fjacob=function(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb) {
   # numerical calculation of full jacobian dres_dparam
   # we variate param one by one and recalculate the whole residual
   # for each variation. The result is returned as matix.
   np=length(param);
   p0=param;
   dp=0.00001;
   for (i in 0:np) {
      p1=p0;
      if (i>0) {
         p1[i]=p0[i]+dp;
      }
      rres=cumo_resid(p1, cjac=F, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, ir2isc, ifmn, fmn, measinvvar, invfmnvar, spAb);
      if (i > 0) {
         dr_dp=cbind(dr_dp, (rres$res-r0)/dp);
      } else {
         r0=rres$res;
         dr_dp=matrix(0, nrow=length(r0), ncol=0);
      }
   }
   return(dr_dp);
}
Tiso2cumo=function(len) {
   if (len<0) {
      return(FALSE);
   }
   if (len==0) {
      return(matrix(1,1,1));
   }
   # recursive call for len>1
   T=Tiso2cumo(len-1);
   return(rbind(cbind(T,T),cbind(diag(0,NROW(T)),T)));
}
Tcumo2iso=function(len) {
   if (len<0) {
      return(FALSE);
   }
   if (len==0) {
      return(matrix(1,1,1));
   }
   # recursive call for len>1
   T=Tcumo2iso(len-1);
   return(rbind(cbind(T,-T),cbind(diag(0,NROW(T)),T)));
}
Tiso2mass=function(len) {
   mass=matrix(0, len+1, 2**len);
   for (i in 0:(2**len-1)) {
      s=sumbit(i);
      mass[s+1,i+1]=1;
   }
   return(mass);
}
Vcumo2iso0=function(len) {
   # coefficients of first row of matrix Tcumo2iso
   # giving the conversion to isotopomer of weight 0
   if (len<0) {
      return(FALSE);
   }
   if (len==0) {
      return(c(1));
   }
   # recursive call for len>1
   V=Vcumo2iso0(len-1);
   return(c(V,-V));
}
sumbit=function(i) {
   i=as.integer(i);
   res=0;
   movi=1;
   while (movi<=i) {
      res=res+(bitAnd(i,movi)>0);
      movi=movi*2;
   }
   return(res);
}
cumo2mass=function(x, sep=":") {
   # convert cumomer vector to mass vectors

   # separate cumos by name and order by weight
   res=c();
   n=length(x);
   nm_x=names(x);
   if (length(nm_x)!=n) {
      return();
   }
   metabs=c(); # unique metab names
   spl=unlist(strsplit(nm_x, sep, fixed=T));
   if (length(spl) != 2*n) {
      # badly formed names
      return(NULL);
   }
   i=1:n;
   icumo=as.integer(spl[2*i]);
   metabs=spl[2*i-1];
   umetabs=union(metabs, NULL);
#cat("metabs:\n");
#print(metabs);
#cat("tbl:\n");
#print(tbl);
   # extract, order and convert each metab vector
   for (metab in umetabs) {
#      cat(paste(metab,":\n",sep=""));
      im=metabs==metab;
#print(d);
      o=order(icumo[im]);
      # ordered cumomer vector with #0==1 component
      vcumo=c(1,x[im][o]);
      clen=log2(length(vcumo));
      # check that we have all components
      sb=sumbit(max(icumo[im]));
      if (!isTRUE(all.equal(sb, clen))) {
         next;
      }
      # mass vector
      mass=c(Tiso2mass(clen)%*%(Tcumo2iso(clen)%*%vcumo));
      names(mass)=paste(metab, "+", 0:clen, sep="");
      res=c(res, mass);
   }
   return(res);
}
cumo2lab=function(x) {
   # converts cumomer vector to fraction of labeled isotopomer 1-i#0
   # separate cumos by name and order by weight
   n=length(x);
   nm_x=names(x);
   if (length(nm_x)!=n) {
      return();
   }
   metabs=c(); # unique metab names
   spl=unlist(strsplit(nm_x,":"));
   i=1:n;
   icumo=as.integer(spl[2*i]);
   metabs=spl[2*i-1];
   umetabs=union(metabs, NULL);
#cat("metabs:\n");
#print(metabs);
#cat("tbl:\n");
#print(tbl);
   # extract, order and convert each metab vector
   res=c();
   for (metab in umetabs) {
#      cat(paste(metab,":\n",sep=""));
      im=metabs==metab;
#print(d);
      o=order(icumo[im]);
      # ordered cumomer vector with #0==1 component
      vcumo=c(1,x[im][o]);
      clen=log2(length(vcumo));
      # labeled fraction
      lab=1-Vcumo2iso0(clen)%*%vcumo;
      names(lab)=metab;
      res=c(res, lab);
   }
   return(res);
}
cumo_gradj=function(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb) {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: mult jacobian by residual 2*jac*resid*invvar

   # grad=c(2*t(dr_df%*%df_dff)*(measinvvar*rescumo)+
   #    t(drf_dff)*(invfmnvar*resfl), 2*(t(dr_dw)%*%rescumo))

   # gradient
   grad=2*t(jx_f$ures*c(measinvvar, invfmnvar))%*%jx_f$udr_dp;
   return(c(grad));
}
# cost function for donlp2 solver
cumo_fn=function(p) {
   return(cumo_cost(p, nb_f, nm, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb))
}
cumo_dfn=function(p) {
   return(cumo_gradj(p, nb_f, nm, nb_rw, nb_rcumos, invAfl, p2bfl, bp, fc, xi, irmeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb));
}
attr(cumo_fn, "gr")=cumo_dfn;
#cumo_fn@gr=cumo_dfn;
cumo_jacob=function(param, nb_f, nm, nb_w, nb_cumos,
   invAfl, p2bfl, bp, fc, xi,
   imeas, measmat, measvec, ir2isc) {
   # calculate jacobian dmeas_dparam and some annexe matricies
   # without applying invvar matrix
   # The result is stored in a global list jx_f.
   if (!all(param==jx_f$param)) {
      # Normally it must not be. The cost function must be already
      # called by this moment.
      # But for some strange reason it didn't happen.
      # So let recalculate the cost and some matricies
      #res=cumo_cost(param, nb_f, nm, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, xi, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, spAb);
      stopifnot(identical(param, jx_f$param));
   }
   ah=1.e-10; # a heavyside parameter to make it derivable in [-ah; ah]
   names(param)=nm_par;
   # grad=c(2*t(dr_df%*%df_dff)*(measinvvar*rescumo)+
   #    t(drf_dff)*(invfmnvar*resfl), 2*(t(dr_dw)%*%rescumo))
   
   # store it for all flux covariance
   jx_f$df_dff<<-df_dff(param, jx_f$flnx);
   
   # measured fluxes derivation
   #dfm_dff=matrix(0., length(nm_fmn), length(nm_ff));
   #dimnames(dfm_dff)=list(nm_fmn, nm_ff);
   #for (nm_y in nm_fmn) {
   #   nm_arr=strsplit(nm_y, "\\.")[[1]];
   #   if (nm_arr[1] == "f") {
   #      # measured flux is free => trivial derivation
   #      dfm_dff[nm_y, nm_y]=1.;
   #   } else if (nm_arr[1] == "d") {
   #      # measured flux is dependent flux
   #      dfm_dff[nm_y,]=dfl_dff[nm_y,];
   #   }
   #}
   # store flux part of jacobian for sensitivity matrix
   jx_f$dfm_dff<<-dfm_dff();
   if (nb_sc > 0) {
      # scale factor part is just zero
      jx_f$dfm_dp<<-cbind(jx_f$dfm_dff, matrix(0, nrow=nrow(jx_f$dfm_dff), ncol=nb_sc));
   } else {
      jx_f$dfm_dp<<-jx_f$dfm_dff;
   }
   
   dmx_df=measmat[,-dim(measmat)[2],drop=F]%*%jx_f$x_f[imeas,];
   drc_dff=c(1.,param)[ir2isc]*(dmx_df%*%jx_f$df_dff);
   
   # store isotope part of jacobian for sensitivity matrix
   jx_f$drc_dff<<-drc_dff;
   
   # add flux measures
   jx_f$dr_dff<<-rbind(drc_dff, jx_f$dfm_dff);
   
   # scale factor part
   sm=jx_f$usimcumom; # just cumomer measure part
   z=rep(0., NROW(sm));
   # each column is fullfilled with a part of residual vector
   # corresponding to a given scale parameter
   if (nb_sc > 0) {
      jx_f$drc_sc<<-apply(t(nb_ff+1+(1:nb_sc)), 2, function(isc){i=ir2isc==isc; v=z; v[i]=sm[i]; v;});
      jx_f$udr_dp<<-cbind(jx_f$dr_dff, rbind(jx_f$drc_sc, matrix(0, nrow=nb_fmn, ncol=nb_sc)));
   } else {
      jx_f$udr_dp<<-jx_f$dr_dff;
   }
   #dimnames(jx_f$udr_dp)=list(names(jx_f$res), nm_par);
   return(NULL);
}
fwrv2Ab=function(fwrv, spAb, incu, nm_rcumo) {
   # calculate sparse A and its rhs b from fields of the list spAb
   # according to conventions explained in comments to python function
   # netan2Abcumo_sp() generating spAb
   # return a list A, b
   # 2010-10-17 sokol
   a=spAb$tA;
   a@x[]=0;
   n=nrow(a);
   res<-.Fortran("f2a",
      fi=as.matrix(spAb$f2ta),
      islot=a@i,
      n=as.integer(ncol(spAb$f2ta)),
      fwrv=as.double(fwrv),
      x=a@x,
      nx=as.integer(length(a@x)),
      NAOK=TRUE,
      DUP=FALSE
   );
   # sum off-diagonal terms
   # and add fluxes from b
   fb=rep(0., n);
   b=rep(0., n);
   res<-.Fortran("f2b",
      bfpr=as.matrix(spAb$bfpr),
      n=as.integer(ncol(spAb$bfpr)),
      fwrv=as.double(fwrv),
      incu=as.double(incu),
      fb=as.double(fb),
      b=as.double(b),
      NAOK=TRUE,
      DUP=FALSE
   );
   diag(a)=-apply(a, 2, sum)-fb;
   A=-t(a)
   dimnames(A)=list(nm_rcumo, nm_rcumo);
   b=-b
   names(b)=nm_rcumo;
   return(list(A=A, b=b));
}
fwrv2Abr=function(fwrv, spAbr, incu, nm_rcumo, getb=T) {
   # calculate sparse A and b (in A*x=b where x is cumomer vector)
   # from fields of the list spAbr
   # according to conventions explained in comments to python function
   # netan2Abcumo_spr() generating spAbr
   # return a list A, b
   # 2012-02-21 sokol
   
   # construct off-diagonal terms of a
   a_pre=spAbr$a_pre
   a=spAbr$tA;
   if (nrow(a_pre) > 0) {
      a_pre@x=fwrv[spAbr$ind_fa]
      a@x=Matrix::colSums(a_pre)
   }
   
   # get just fluxes in b
   b_pre=spAbr$b_pre;
   b_pre@x=fwrv[spAbr$ind_fb]
   b=spAbr$b
   b@x=Matrix::colSums(b_pre)
   
   # sum off-diagonal terms
   # and add fluxes from b
   diag(a)=as.numeric(-Matrix::colSums(a))-as.numeric(b);
   A=t(a)
   dimnames(A)=list(nm_rcumo, nm_rcumo);
   
   # construct a complete b vector
   if (getb) {
      b_pre@x=b_pre@x*incu[spAbr$ind_x1]*incu[spAbr$ind_x2]
      b@x=-colSums(as.matrix(b_pre))

      b=as.numeric(b)
      names(b)=nm_rcumo;
      return(list(A=A, b=b));
   } else {
      return(list(A=A, b=NULL));
   }
}
fx2j=function(fwrv, spAb, incu) {
   # calculate sparse j_rhs and b_x from fields of the list spAb
   # according to conventions explained in comments
   # to python function netan2Abcumo_sp() generating spAb
   # Return a list j_rhs, b_x
   # 2010-10-22 sokol
   fwrv=as.double(fwrv);
   incu=as.double(incu);
   a_fx=spAb$ta_fx;
   a_fx@x[]=0.;
   res<-.Fortran("x2a_fx",
      xi=as.matrix(spAb$x2ta_fx),
      n=as.integer(ncol(spAb$x2ta_fx)),
      islot=a_fx@i,
      x0=c(0., incu),
      x=a_fx@x,
      nx=as.integer(length(a_fx@x)),
      NAOK=TRUE,
      DUP=FALSE
   );
   b_f=spAb$tb_f;
   b_f@x[]=0.;
   res<-.Fortran("x2b_f",
      xi=as.matrix(spAb$x2tb_f),
      n=as.integer(ncol(spAb$x2tb_f)),
      islot=b_f@i,
      incu=incu,
      x=b_f@x,
      nx=as.integer(length(b_f@x)),
      NAOK=TRUE,
      DUP=FALSE
   );
   b_x=spAb$tb_x;
   b_x@x[]=0.;
   res<-.Fortran("fx2b_x",
      fxi=as.matrix(spAb$fx2tb_x),
      n=as.integer(ncol(spAb$fx2tb_x)),
      islot=b_x@i,
      fwrv=fwrv,
      incu=incu,
      x=b_x@x,
      nx=as.integer(length(b_x@x)),
      NAOK=TRUE,
      DUP=FALSE
   );
   return(list(j_rhsw=-t(b_f-a_fx), b_x=-t(b_x)));
}

fx2jr=function(fwrv, spAbr, incu, incup=NULL) {
   # calculate sparse j_rhs and b_x from fields of the list spAbr
   # according to conventions explained in comments
   # to python function netan2Abcumo_spr() generating spAbr
   # Return a list j_rhs, b_x
   # 2012-02-22 sokol
   
   x0=c(0., incu)
   a_fx=spAbr$ta_fx
   # form a_fx_pre
   x2ta_fx=spAbr$x2ta_fx
   a_fx_pre=Matrix(0, ncol=length(a_fx@i), nrow=max(x2ta_fx[,1])+1)
   a_fx_pre@i=as.integer(x2ta_fx[,1])
   inz=which(diff(c(a_fx_pre@i, 0))<=0)
   a_fx_pre@p=as.integer(c(0, cumsum(a_fx_pre@i[inz]+1)))
   a_fx_pre@x=x0[x2ta_fx[,2]]-x0[x2ta_fx[,3]]
   # calculate @x slot of a_fx
   a_fx@x=colSums(as.matrix(a_fx_pre))
   
   # prepare b_f
   b_f=spAbr$tb_f
   x2tb_f=spAbr$x2tb_f
   b_f@x=incu[x2tb_f[,1]]*incu[x2tb_f[,2]]
   
   # prepare b_x
   b_x=spAbr$tb_x
   fx2tb_x=spAbr$fx2tb_x
   if (nrow(b_x) > 0 && all(dim(fx2tb_x) > 0)) {
      # form b_x_pre
      b_x_pre=Matrix(0, ncol=length(b_x@i), nrow=max(fx2tb_x[,1])+1)
      b_x_pre@i=as.integer(fx2tb_x[,1])
      inz=which(diff(c(b_x_pre@i, 0))<=0)
      b_x_pre@p=as.integer(c(0, cumsum(b_x_pre@i[inz]+1)))
      b_x_pre@x=fwrv[fx2tb_x[,2]]*incu[fx2tb_x[,3]]
      # calculate @x slot of b_x
      b_x@x=Matrix::colSums(b_x_pre)
   }
   if (!is.null(incup)) {
      # calculate first derivative in time
      xp0=c(0., incup)
      # form a_fx_pre
      a_fx_pre@x=xp0[x2ta_fx[,2]]-xp0[x2ta_fx[,3]]
      # calculate @x slot of a_fxp
      a_fxp=a_fx
      a_fxp@x=colSums(as.matrix(a_fx_pre))

      # prepare b_fp
      b_fp=spAbr$tb_f
      b_fp@x=incup[x2tb_f[,1]]*incu[x2tb_f[,2]]+incu[x2tb_f[,1]]*incup[x2tb_f[,2]]

      # prepare b_xp
      b_xp=spAbr$tb_x
      if (nrow(b_x) > 0 && all(dim(fx2tb_x) > 0)) {
         # form b_x_pre
         b_x_pre@x=fwrv[fx2tb_x[,2]]*incup[fx2tb_x[,3]]
         # calculate @x slot of b_x
         b_xp@x=Matrix::colSums(b_x_pre)
      }
      j_rhswp=-t(b_fp+a_fxp)
      b_xp=-t(b_xp)
   } else {
      j_rhswp=NULL
      b_xp=NULL
   }
   return(list(j_rhsw=-t(b_f+a_fx), b_x=-t(b_x), j_rhswp=j_rhswp, b_xp=b_xp))
}

put_inside=function(param, ui, ci) {
   # put param inside of feasible domain delimited by u%*%param >= ci
   mes=""
   ineq=ui%*%param-ci;
   if (all(ineq>1.e-10)) {
      # nothing to do, already inside and well inside
      return(param)
   }
   dp=ldp(ui, -ineq);
   if (!is.null(dp)) {
      # get new active inequalities
      ineqd=ui%*%(param+dp)-ci
      # check that we are at least at the border and not outside
      if (any(ineqd < -1.e-7)) {
         param=NA
         attr(param, "mes")="Inequality system is ill-conditionned. Failed to solve."
         attr(param, "err")=1
         return(param)
      }
      iact=ineqd<=1.e-10
#print(ineqd[iact])
      # solve an ldp pb to find non decrising direction
      # for active inequalities
      ma=ui[iact,,drop=F]
      na=sum(iact)
      # find closest vector to c(1,1,...) making the direction increasing
      tma=tcrossprod(ma)
      bet=ldp(tma, 1.e-3 - apply(tma, 1, sum))
      if (is.null(bet)) {
         param=param+dp
         attr(param, "mes")="Infeasible constraints for inside direction."
         attr(param, "err")=0
         return(param)
      }
      vn=ma%tmm%(1.+bet)
      vn=vn/norm(vn)
      decr=(ui%*%vn)<0.
      alpha=((-ineqd)/(ui%*%vn))[decr]
      alpha=alpha[alpha>0]
      alpha=0.5*min(0.001, alpha)
      dpn=dp+alpha*vn
      # check that new dp is still inside
      if (any(ui%*%(param+dpn)-ci < 0.)) {
         attr(param, "err")=0 # just a warning
         attr(param, "mes")="Failed to put free parameters strictly inside of the feasible domain. They are left on the border."
         dpn=dp
      }
      names(param)=nm_par;
      if (nb_ff > 0) {
         i=abs(dpn[1:nb_ff])>=1.e-10
         if (sum(i) > 0) {
            tmp=cbind(param[1:nb_ff], param[1:nb_ff]+dpn[1:nb_ff], dpn[1:nb_ff]);
            dimnames(tmp)=list(nm_par[1:nb_ff], c("outside", "inside", "delta"));
            obj2kvh(tmp[i,,drop=F], "Free parameters put inside of feasible domain");
         }
      }
      # move starting point slightly inside of feasible domain
      param=param+c(dpn);
   } else {
      param=NA
      mes="Infeasible inequalities."
      if (!is.null(rownames(ui))) {
         mes=join("\n", c(mes, rownames(ui)))
      }
      attr(param, "mes")=mes;
      attr(param, "err")=1
   }
   return(param)
}
fwrv2sp=function(fwrv, spAbr, x, xp=NULL, gets=TRUE) {
   # calculate s and ds/dt in A*x+s (where x is cumomer vector
   # and xp its first derivative in time)
   # from fields of the spAbr
   # according to conventions explained in comments to python function
   # netan2Abcumo_spr() generating spAbr
   # return a list s and sp
   # 2012-03-07 sokol
   
   # construct the sources s and its derivatives (if xp not null)
   # for this weight
   if (gets) {
      s_pre=spAbr$b_pre;
      s_pre@x=fwrv[spAbr$ind_fb]*x[spAbr$ind_x1]*x[spAbr$ind_x2]
      s=spAbr$b
      s@x=Matrix::colSums(s_pre)
   } else {
      s=NULL
   }
   # sp
   if (!is.null(xp)) {
      s_pre=spAbr$b_pre;
      s_pre@x=fwrv[spAbr$ind_fb]*(xp[spAbr$ind_x1]*x[spAbr$ind_x2]+x[spAbr$ind_x1]*xp[spAbr$ind_x2])
      sp=spAbr$b
      sp@x=Matrix::colSums(s_pre)
   } else {
      sp=NULL
   }
   return(list(s=s, sp=sp));
}
expm_cub_step=function(inva, dt, expadt, s1, s2, sp1, sp2, x1) {
   # make a single time step for integration of ode xp=a*x+s
   # according to cubic scheme with exponential matrix
   # return a vector x
   # 2012/03/07 sokol
   
   # auxiliary terms for integration
   dt2=2*dt
   invdt2=1./dt/dt
   tmp1=6*(s2-s1)
   tmp2=(sp1+sp2)*dt2
   # second derivative on the left (at t1)
   sddl=tmp1-tmp2; #+sd[,-nb_ti,drop=F])*(2*dt)
   # second derivative on the right (at t2)
   sddr=(-sddl+sp2*dt2)*invdt2
   sddl=(sddl-sp1*dt2)*invdt2
   # third derivative*dt**3
   sddd=(-2*tmp1+tmp2*3)*invdt2/dt
   tmp=inva%*%sddd
   # term t
   tmp1=s1+inva%*%(sp1+inva%*%(sddl+tmp))
   # term t2
   tmp2=s2+inva%*%(sp2+inva%*%(sddr+tmp))
   # increment
   dx=inva%*%(expadt%*%tmp1-tmp2)
   res=expadt%*%x1+dx
   return(res)
}
expm_lin_step=function(inva, dt, expadt, s1, s2, x1) {
   # make a single time step for integration of ode xp=a*x+s
   # according to linear scheme with exponential matrix
   # return a vector x
   # 2012/04/25 sokol
   
   # auxiliary terms for integration
   tmp=inva%*%((s2-s1)/dt)
   # term t
   tmp1=s1+tmp
   # term t2
   tmp2=s2+tmp
   # increment
   dx=inva%*%(expadt%*%tmp1-tmp2)
   res=expadt%*%x1+dx
   return(res)
}
expm_quad_step=function(inva, dt, expadt, s1, s2, sp1, sp2, x1) {
   # make a single time step for integration of ode xp=a*x+s
   # according to quadratic scheme with exponential matrix
   # return a vector x
   # 2012/04/30 sokol
   
   # auxiliary terms for integration
   tmp=inva%*%((sp2-sp1)/dt)
   # term t
   tmp1=s1+inva%*%(sp1+tmp)
   # term t2
   tmp2=s2+inva%*%(sp2+tmp)
   # increment
   dx=inva%*%(expadt%*%tmp1-tmp2)
   res=expadt%*%x1+dx
   return(res)
}
df_dff=function(param, flnx) {
   ah=1.e-10; # a heavyside parameter to make it derivable in [-ah; ah]
   nb_fwrv=length(nm_fwrv)
   i_fln=grep("d.n.", names(flnx), fixed=T)
   i_flx=grep("d.x.", names(flnx), fixed=T)
   i_ffn=grep("f.n.", nm_par, fixed=T)
   i_ffx=grep("f.x.", nm_par, fixed=T)
   nb_ff=length(i_ffn)+length(i_ffx)
   # dfl_dff=invAfl%*%p2bfl; is already calculated
   # df_dfl=fw-rv derived by dependent n-x where x in 0;1
   df_dfl=matrix(0., length(nm_fwrv), length(flnx));
   #dimnames(df_dfl)=list(nm_fwrv, nm_fl);
   # df_dffd=fw-rv derived directly by free n-x where x in 0;1
   df_dffd=matrix(0., length(nm_fwrv), nb_ff);
   #dimnames(df_dffd)=list(nm_fwrv, c(nm_ffn, nm_ffx));
#browser();
   # derivation by dependent fluxes
   # net part
   net=flnx[i_fln];
   hnet=Heaviside(net);
   i=abs(net)<ah;
   hnet[i]=net[i]/ah;
   # xch part
   xch=flnx[i_flx];
   xch=1./(1.-xch)**2;
   
   # forward fluxes
   df_dfl[cfw_fl]=c(hnet, xch);
   # reverse fluxes
   df_dfl[crv_fl]=c(hnet-1., xch);
   
   # derivation by free fluxes
   # net part
   net=param[i_ffn];
   hnet=Heaviside(net);
   i=abs(net)<ah;
   hnet[i]=net[i]/ah;
   # xch part
   xch=param[i_ffx];
   xch=1./(1.-xch)**2;
   
   # forward fluxes
   df_dffd[cfw_ff]=c(hnet, xch);
   # reverse fluxes
   df_dffd[crv_ff]=c(hnet-1., xch);
   
   res=df_dfl%*%dfl_dff+df_dffd;
   dimnames(res)=list(nm_fwrv, names(param)[c(i_ffn, i_ffx)]);
   return(res)
}
dfm_dff=function() {
   # measured fluxes derivation
   res=matrix(0., length(nm_fmn), length(nm_ff));
   dimnames(res)=list(nm_fmn, nm_ff);
   # derivate free measured fluxes (trivial)
   i=grep("f.n.", nm_fmn, fixed=T, value=T)
   res[cbind(i, i)]=1.
   # derivate dependent measured fluxes
   i=grep("d.n.", nm_fmn, fixed=T, value=T)
   res[i,]=dfl_dff[i,];
   return(res)
}
plot_ti=function(ti, x, m=NULL, ...) {
   # plot time curse curves x[icurve, itime] and points from m
   # x and m are supposed to have the same dimension and organization
   nm=rownames(x)
   nb_curve=nrow(x)
   if (is.null(nm)) {
      nm=1:nb_curve
   } else {
      o=order(nm)
      nm=nm[o]
      x=x[o,,drop=F]
      if (!is.null(m) && nrow(x)==nrow(m)) {
         m=m[o,,drop=F]
      }
   }
   # x and m may have different time moments
   if (is.null(m)) {
      tim=ti
   } else {
      tim=as.numeric(colnames(m))
   }
   plot(range(ti, tim), range(c(x,m)), t="n", ylab="Labeling", xlab="Time", ...)
   matplot(ti, t(x), t="l", lty=1:nb_curve, col=1:nb_curve, lwd=2, add=T)
   legend("topright", legend=nm, lty=1:nb_curve, col=1:nb_curve, lwd=2, cex=0.75)
   if (!is.null(m)) {
      # plot measured points
      for (i in 1:nrow(m)) {
         points(tim, m[i,], col=i, cex=0.5, t="b", lwd=0.5)
      }
   }
}
get_usm=function(f) {
   # return list of ti, usm from a _res.kvh file f

   # get skip and end number in the kvh
   res=system(sprintf("grep -n 'simulated unscaled labeling measurements' '%s'", f), intern=T)
   nstart=as.numeric(strsplit(res, ":")[[1]][1])
   res=system(sprintf("grep -n 'simulated scaled labeling measurements' '%s'", f), intern=T)
   nend=as.numeric(strsplit(res, ":")[[1]][1])
   d=as.matrix(read.table(f, header=T, skip=nstart, row.names=2, sep="\t", nrows=nend-nstart-2, check=F, comment=""))[,-1]
   ti=as.numeric(colnames(d))
   o=order(rownames(d))
   return(list(ti=ti, usm=d[o,,drop=F]))
}
get_labcin=function(f, nm_meas) {
   # get labeling cinetic data form file f
   # with rows matching at the best nm_meas
   d=as.matrix(read.table(f, header=T, row.names=1, sep="\t", check=F))
   # strip the last field (row id) in nm_meas and make it correspond
   # to rownames of d
   nm_r=rownames(d)
   nm=sapply(nm_meas, function(s) {
      v=strsplit(s, ":")[[1]]
      v[length(v)]=""
      pmatch(paste(v, sep=":", collapse=":"), nm_r)
   })
   return(d[nm,])
}
get_hist=function(f, v) {
   # return list of ti, usm from a _res.kvh file f

   # get skip and end number in the kvh
   res=system(sprintf("grep -m 1 -n 'history' '%s'", f), intern=T)
   nstart=as.numeric(strsplit(res, ":")[[1]][1])
   res=system(sprintf("tail -n +%d '%s' 2>&1| grep -m 1 -n '%s$'", nstart, f, v), intern=T, ignore.stderr=T)
   nstart=nstart+as.numeric(strsplit(res, ":")[[1]][1])
   res=system(sprintf("tail -n +%d '%s' 2>&1| grep -m 1 -n -E ^$'\\t?\\t\\w'", nstart, f), intern=T, ignore.stderr=T)
   nend=nstart+as.numeric(strsplit(res, ":")[[1]][1])
   d=as.matrix(read.table(f, header=F, skip=nstart, row.names=4, sep="\t", nrows=nend-nstart-2, check=F))[,-(1:3)]
   return(d)
}
get_matrix=function(f, v) {
   # return a matrix stored in the field v of a kvh file f
   
   # get start line number by grep successively all fields in v
   # and the indent
   indent=0
   nstart=0
   res=""
   lapply(v, function(it) {
      res<<-system(sprintf("tail -n +%d '%s' 2>&1| grep -m 1 -n '\\t*%s$'", nstart+1, f, it), intern=T, ignore.stderr=T)
      spl=strsplit(res, ":")[[1]]
      nstart<<-nstart+as.numeric(spl[1])
      indent<<-which(nchar(strsplit(spl[2], "\t")[[1]])>0)[1]-1
      return(NULL)
   })
   # get end number in the kvh
   wop=options("warn")$warn
   options(warn=-1)
   res=system(sprintf("tail -n +%d '%s' 2>&1| grep -m 1 -n -E ^'%s\\w'", nstart+1, f, join("", rep("\t", indent))), intern=T)
   options(warn=wop)
   if (length(res) == 0) {
      nrows=-1
   } else {
      nrows=as.numeric(strsplit(res, ":")[[1]][1])-1
   }
   d=as.matrix(read.table(f, header=F, skip=nstart, row.names=indent+2, sep="\t", nrows=nrows, check=F, comment=""))[,-(1:(indent+1))]
   if (rownames(d)[1]=="row_col") {
      colnames(d)=d[1,]
      d=d[-1,,drop=F]
      dn=dimnames(d)
      options(warn=-1)
      num<-as.numeric(d)
      options(warn=wop)
      if (!any(is.na(num))) {
         d=matrix(num, nrow=nrow(d))
         dimnames(d)=dn
      }
   }
   return(d)
}
