#DEBUG=0; # 1 to enable debugging information, 0=disable
#TIMEIT=0; # 1 to enable time printing at some stages
if (length(find("TIMEIT")) && TIMEIT) {
   cat("load    : ", date(), "\n", sep="");
}
jx_f=list();
suppressPackageStartupMessages(library(bitops));
suppressPackageStartupMessages(library(MASS)); # for generalized inverse
suppressPackageStartupMessages(library(fUtilities)); # for Heaviside function
suppressPackageStartupMessages(library(nnls)); # for non negative least square
#suppressPackageStartupMessages(library(lattice)); # to keep Matrix silent
suppressPackageStartupMessages(library(Matrix, warn=F, verbose=F)); # for sparse matrices
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
         jx_f$jacobian<<-jacobian;
      } else {
         jacobian=jx_f$jacobian
      }
   } else {
      jacobian=NULL;
   }

   return(list(res=jx_f$res, fallnx=jx_f$fallnx,
      jacobian=jacobian));
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
   flnx=c(invAfl%*%(p2bfl%*%param[1:nb_f$nb_ff]+c(bp)));
   names(flnx)=nm$flnx;
   fallnx=c(dfc2fallnx(nb_f, flnx, param, fc));
   names(fallnx)=nm$fallnx;
   fwrv=c(fallnx2fwrv(fallnx));
   names(fwrv)=nm$fwrv;
   if (DEBUG) {
      write.matrix(p2bfl%*%param[1:nb_f$nb_ff]+bp, file="dbg_bfl.txt", sep="\t");
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
   n_2=nb_fwrv/2;
   names(param)=nm_par;
   # grad=c(2*t(dr_df%*%df_dff)*(measinvvar*rescumo)+
   #    t(drf_dff)*(invfmnvar*resfl), 2*(t(dr_dw)%*%rescumo))
   # dfl_dff=invAfl%*%p2bfl; is already calculated
   # df_dfl=fw-rv derived by dependent n-x where x in 0;1
   df_dfl=matrix(0., length(jx_f$fwrv), length(jx_f$flnx));
   #dimnames(df_dfl)=list(nm_fwrv, nm_fl);
   # df_dffd=fw-rv derived directly by free n-x where x in 0;1
   df_dffd=matrix(0., length(jx_f$fwrv), nb_ff);
   #dimnames(df_dffd)=list(nm_fwrv, c(nm_ffn, nm_ffx));
#browser();
   # derivation by dependent fluxes
   # net part
   net=jx_f$flnx[nm_fln];
   hnet=Heaviside(net);
   i=abs(net)<ah;
   hnet[i]=net[i]/ah;
   # xch part
   xch=jx_f$flnx[nm_flx];
   xch=1./(1.-xch)**2;
   
   # forward fluxes
   df_dfl[cfw_fl]=c(hnet, xch);
   # reverse fluxes
   df_dfl[crv_fl]=c(hnet-1., xch);
   
   # derivation by free fluxes
   # net part
   net=param[nm_ffn];
   hnet=Heaviside(net);
   i=abs(net)<ah;
   hnet[i]=net[i]/ah;
   # xch part
   xch=param[nm_ffx];
   xch=1./(1.-xch)**2;
   
   # forward fluxes
   df_dffd[cfw_ff]=c(hnet, xch);
   # reverse fluxes
   df_dffd[crv_ff]=c(hnet-1., xch);
   
   df_dff=df_dfl%*%dfl_dff+df_dffd;
   dimnames(df_dff)=list(nm_fwrv, nm_ff);
   # store it for all flux covariance
   jx_f$df_dff<<-df_dff;
   
   # measured fluxes derivation
   dfm_dff=matrix(0., length(nm_fmn), length(nm_ff));
   dimnames(dfm_dff)=list(nm_fmn, nm_ff);
   for (nm_y in nm_fmn) {
      nm_arr=strsplit(nm_y, "\\.")[[1]];
      if (nm_arr[1] == "f") {
         # measured flux is free => trivial derivation
         dfm_dff[nm_y, nm_y]=1.;
      } else if (nm_arr[1] == "d") {
         # measured flux is dependent flux
         dfm_dff[nm_y,]=dfl_dff[nm_y,];
      }
   }
   # store flux part of jacobian for sensitivity matrix
   jx_f$dfm_dff<<-dfm_dff;
   if (nb_sc > 0) {
      # scale factor part is just zero
      jx_f$dfm_dp<<-cbind(jx_f$dfm_dff, matrix(0, nrow=nrow(jx_f$dfm_dff), ncol=nb_sc));
   } else {
      jx_f$dfm_dp<<-jx_f$dfm_dff;
   }
   
   dmx_df=measmat[,-dim(measmat)[2],drop=F]%*%jx_f$x_f[imeas,];
   drc_dff=c(1.,param)[ir2isc]*(dmx_df%*%df_dff);
   
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
   a_pre@x=fwrv[spAbr$ind_fa]
   a=spAbr$tA;
   a@x=apply(a_pre, 2, sum)
   
   # get fluxes in b
   b_pre=spAbr$b_pre;
   b_pre@x=fwrv[spAbr$ind_fb]
   b=spAbr$b
   b@x=apply(b_pre, 2, sum)
   
   # sum off-diagonal terms
   # and add fluxes from b
   diag(a)=as.numeric(-apply(a, 2, sum))-as.numeric(b);
   A=t(a)
   dimnames(A)=list(nm_rcumo, nm_rcumo);
   
   # construct a complete b vector
   if (getb) {
      b_pre@x=b_pre@x*incu[spAbr$ind_x1]*incu[spAbr$ind_x2]
      b@x=-apply(b_pre, 2, sum)

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

fx2jr=function(fwrv, spAbr, incu) {
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
   a_fx@x=apply(a_fx_pre, 2, sum)
   
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
      b_x@x=apply(b_x_pre, 2, sum)
   }
   return(list(j_rhsw=-t(b_f+a_fx), b_x=-t(b_x)))
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
      tmp=cbind(param[1:nb_ff], param[1:nb_ff]+dpn[1:nb_ff], dpn[1:nb_ff]);
      i=abs(dpn[1:nb_ff])>=1.e-10
      dimnames(tmp)=list(nm_par[1:nb_ff], c("outside", "inside", "delta"));
      obj2kvh(tmp[i,,drop=F], "Free parameters put inside of feasible domain");
      # move starting point slightly inside of feasible domain
      param=param+c(dpn);
   } else {
      param=NA
      attr(param, "mes")="Infeasible inequalities.";
      attr(param, "err")=1
   }
   return(param)
}
