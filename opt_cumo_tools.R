DEBUG=0;
jx_f=list();
library(bitops);
library(MASS); # for generalized inverse
library(fUtilities); # for Heaviside function

trisparse_solv=function(A, b, w, method="dense") {
   # solve A*x=b where A=tridiag(Al,Ac,Au)+s*e^t and b is dense
   if (method=="dense") {
#cat("trisparse: A,b\n");
#print(A);
#print(b);
      if (DEBUG) {
         write.matrix(A,file=paste("dbg_Acumo_d_",w,".txt", sep=""),sep="\t");
      }
      qrA=qr(A);
      x=try(solve(qrA,b));
      if (inherits(x, "try-error")) {
         # matrix seems to be singular
         # switch to Moore-Penrose inverse
         if ((exists("control") && control$trace) || DEBUG) {
            cat("switch to generalized inverse at weight ", w, "\n", sep="");
         }
         x=ginv(A)%*%b;
      }
      return(list(x=x, qrA=qrA));
   } else if (method=="sparse") {
      # sparse
      # fulfill a matrix
      require(Matrix);
      As=Matrix(A);
      x=solve(As,b);
#q=qr(A);
#print("qr");
#print(q@V);
#print(q@R);
#print(q@p);
      return(x);
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
cumo_resid=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo", fj_rhs="frj_rhs") {
#cat("resid: \n")
#print(nb_f);
#print(nb_w);
#print(param);
#print(p2bfl);
   # find x for all weights
   lres=param2fl_x(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun, fj_rhs)
#print(imeas);
   # find simulated scaled measure vector scale*(measmat*x)
   simvec=c(1.,param)[ir2isc]*(measmat%*%c(lres$x[imeas],1.));
   # diff between simulated and measured
   return(list(res=(simvec-measvec), fallnx=lres$fallnx));
}
cumo_cost=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo", fj_rhs="frj_rhs") {
#cat("cost: ");
#cat("list nb_f\n")
#print(nb_f);
#print(nb_w);
#cat("param\n");
#print(param);
    resl=cumo_resid(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun, fj_rhs);
   res=resl$res;
   fallnx=resl$fallnx;
   # flux residuals
   resfl=fallnx[ifmn]-fmn;
   fn=sum(res*res*measinvvar)+sum(resfl*resfl*invfmnvar);
   if (DEBUG) {
      write.matrix(fn, file="dbg_cost.txt", sep="\t");
   }
   # complete usefull information in global list jx_f
   jx_f$res<<-resl$res;
   jx_f$resfl<<-resfl;
   return(fn);
}
cumo_grad=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo", fj_rhs=NULL) {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: forward finite differences f(x+h)-f(x)/h
   # x+h is taken as (1+fact)*x
   fact=1.e-7;
   grad=param; # make place for gradient
   # f(x)
   f=cumo_cost(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun, fj_rhs);
   for (i in 1:length(param)) {
      x=param[i];
      h=x*fact;
      param[i]=x+h;
      if (param[i]==x) {
         # we are too close to zero here
         param[i]=fact;
      }
      fh=cumo_cost(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun, fj_rhs);
      # restore modified param
      param[i]=x;
      grad[i]=(fh-f)/h;
   }
   return(grad);
}
param2fl=function(param, nb_f, invAfl, p2bfl, bp, fc) {
      # claculate all fluxes from free fluxes
#cat("resid: \n")
#cat("nb_f", str(nb_f), "\n");
#cat("nb_w", nb_w, "\n");
#cat("param", param, "\n");
#cat("bp", bp, "\n");
#cat("p2bfl before", str(p2bfl), "\n");
#   p2bfl=matrix(p2bfl, ncol=nb_f$nb_ff);
#cat("p2bfl", str(p2bfl), "\n");
#   invAfl=matrix(invAfl, nrow=nb_f$nb_fl);
#cat("invAfl", invAfl, "\n");
   flnx=invAfl%*%(p2bfl%*%param[1:nb_f$nb_ff]+bp);
#cat("flnx");
#print(flnx);
   fallnx=dfc2fallnx(nb_f, flnx, param, fc);
   fwrv=fallnx2fwrv(fallnx);
   if (DEBUG) {
      write.matrix(p2bfl%*%param[1:nb_f$nb_ff]+bp, file="dbg_bfl.txt", sep="\t");
      n=length(fwrv);
      names(fwrv)=nm_fwrv;
      write.matrix(fwrv, file="dbg_fwrv.txt", sep="\t");
      write.matrix(cbind(1:n,nm_fallnx,fallnx), file="dbg_fallnx.txt", sep="\t");
#cat("fwrv");
#print(fwrv);
   }
   return(list(fallnx=fallnx, fwrv=fwrv, flnx=flnx));
}

param2fl_x=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo", fj_rhs="frj_rhs") {
   # calculate all fluxes from free fluxes
   lf=param2fl(param, nb_f, invAfl, p2bfl, bp, fc);
   # construct the system A*x=b from fluxes
   # and find x for every weight
   # if fj_rhs is not NULL, calculate jacobian x_f
   nb_fwrv=length(lf$fwrv);
   x=numeric(0);
   x_f=matrix(0., nrow=0, ncol=nb_fwrv);
   if (DEBUG) {
      tmp=lf$fwrv;
      names(tmp)=nm_fwrv;
      conct=file("dbg_fwrv.txt", "wb");
      obj2kvh(tmp, "fwrv", conct);
      tmp=lf$fallnx;
      names(tmp)=nm_fallnx
      obj2kvh(tmp, "net-xch", conct);
#      write.matrix(tmp, file="dbg_fwrv.txt", sep="\t");
   }
   for (iw in 1:nb_w) {
      nx=length(x);
      ncumow=nb_cumos[iw];
      A=matrix(0.,ncumow,ncumow);
      b=double(ncumow);
      #fwrv2Abcumo(fl, nf, x, nx, iw, n, A, b)
      res<-.Fortran(fortfun,
         fl=as.double(lf$fwrv),
         nf=nb_fwrv,
         x=as.double(x),
         nx=as.integer(nx),
         iw=as.integer(iw),
         n=as.integer(ncumow),
         A=as.matrix(A),
         b=as.double(b),
         calcA=as.integer(TRUE),
         calcb=as.integer(TRUE),
         NAOK=TRUE,
         DUP=FALSE);
      # solve the system A*x=b;
      if (DEBUG) {
         write.matrix(cbind(A, b=b), file=paste("dbg_cumoAb_",iw,".txt", sep=""), sep="\t");
      }
      lsolv=trisparse_solv(A, b, iw, method="dense");
      xw=lsolv$x;
      nxw=length(xw);
#fj_rhs(fl, nf, x, x_f, nx, iw, n, j_rhs)
      if (length(fj_rhs) && nchar(fj_rhs)) {
         # calculate jacobian x_f
         # first, calculate right hand side for jacobian solve
         j_rhsw=matrix(0., nxw, nb_fwrv);
         b_x=matrix(0., nxw, nx);
         res<-.Fortran(fj_rhs,
         fl=as.double(lf$fwrv),
         nf=nb_fwrv,
         x=as.double(x),
         xw=as.double(xw),
         x_f=as.double(x_f),
         nx=as.integer(nx),
         nxw=as.integer(nxw),
         iw=as.integer(iw),
         j_rhs=as.matrix(j_rhsw),
         b_x=as.matrix(b_x),
         NAOK=TRUE,
         DUP=FALSE);
         if (DEBUG) {
            write.matrix(cbind(j_rhsw,b_x), file=paste("dbg_j_rhs_",iw,".txt", sep=""), sep="\t");
            browser();
         }
         if (iw > 1) {
            x_f=rbind(x_f, solve(lsolv$qrA, j_rhsw+b_x%*%x_f));
         } else {
            x_f=rbind(x_f, solve(lsolv$qrA, j_rhsw));
         }
      }
      # bind vectors and matrices
      x=c(x,xw);
   }
#print(x);
   # store usefull information in global list jx_f
   jx_f$param<<-param;
   jx_f$x_f<<-x_f;
   jx_f$fallnx<<-lf$fallnx;
   jx_f$fwrv<<-lf$fwrv;
   jx_f$flnx<<-lf$flnx;
   return(append(list(x=x, x_f=x_f), lf));
}
num_jacob=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo") {
   # numerical calculation of jacobian dx_df
   # we variate fvrw one by one and recalculate the whole x
   # each variation. The result is returned as matix.
   # The first column of matrix is just cumomer vector corresponding
   # to non perturbed fluxes.
   
   # calculate all fluxes from free fluxes
   lf=param2fl(param, nb_f, invAfl, p2bfl, bp, fc);
   f0=lf$fwrv;
   dfl=0.0001;
   x_f=matrix(0., ncol=0, nrow=sum(nb_cumos));
   nb_fwrv=length(f0);
   for (i in 0:length(f0)) {
      f1=f0;
      if (i > 0) {
         f1[i]=f0[i]+dfl;
      }
      # construct the system A*x=b from fluxes
      # and find x for every weight
      # if fj_rhs is not NULL, calculate jacobian x_f
      x=numeric(0);
      for (iw in 1:nb_w) {
         nx=length(x);
         ncumow=nb_cumos[iw];
         A=matrix(0.,ncumow,ncumow);
         b=double(ncumow);
         #fwrv2Abcumo(fl, nf, x, nx, iw, n, A, b)
         res<-.Fortran(fortfun,
            fl=as.double(f1),
            nf=nb_fwrv,
            x=as.double(x),
            nx=as.integer(nx),
            iw=as.integer(iw),
            n=as.integer(ncumow),
            A=as.matrix(A),
            b=as.double(b),
            calcA=as.integer(TRUE),
            calcb=as.integer(TRUE),
            NAOK=TRUE,
            DUP=FALSE);
         # solve the system A*x=b;
         lsolv=trisparse_solv(A, b, iw, method="dense");
         xw=lsolv$x;
         nxw=length(xw);
         # bind vectors and matrices
         x=c(x,xw);
      }
      if (i > 0) {
         x_f=cbind(x_f, (x-x0)/dfl);
      } else {
         x0=x;
         x_f=cbind(x_f, x0);
      }
   }
   return(x_f);
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
cumo2mass=function(x) {
   # convert cumomer vector to mass vectors

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
cumo_gradj=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo", fj_rhs="frj_rhs") {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: jacobian by implicite derivation (cf. Wiechert)
   # a global variable jx_f is used to store the jacobian x_f
   # during the call to cost function (to take advantage of matrix
   # factorisation made at this time).
   
   stopifnot(param==jx_f$param); # must not be. The cost function must be alread called at this time
   # grad=c(2*t(dr_df%*%df_dff)*(measinvvar*rescumo)+
   #    t(drf_dff)*(invfmnvar*resfl), 2*(t(dr_dw)%*%rescumo))
   # dfl_dff=invAfl%*%p2bfl; is already calculated
   # df_dfl=fw-rv derived by dependent n-x where x in 0;1
   df_dfl=matrix(0., length(jx_f$fwrv), length(jx_f$flnx));
   dimnames(df_dfl)=list(nm_fwrv, nm_fl);
   # df_dffd=fw-rv derived directly by free n-x where x in 0;1
   df_dffd=matrix(0., length(jx_f$fwrv), nb_ff);
   dimnames(df_dffd)=list(nm_fwrv, c(nm_ffn, nm_ffx));
   for (nm_y in nm_fwrv) {
      nm_arr=strsplit(nm_y, "\\.")[[1]];
      reac=nm_arr[2];
      if (nm_arr[1]=="fwd") {
         # derive forward and reverse fluxes
         # by dependent net part
         nm_x=join(".", c("d", "n", reac));
         nm_yr=join(".", c("rev", reac))
         if (nm_x %in% nm_fl) {
            net=jx_f$flnx[nm_x,1];
            hnet=Heaviside(net);
            df_dfl[nm_y, nm_x]=hnet;
            df_dfl[nm_yr, nm_x]=1-hnet;
         }
         # by free net part
         nm_x=join(".", c("f", "n", reac));
         if (nm_x %in% nm_ffn) {
            net=param[nm_x];
            hnet=Heaviside(net);
            df_dffd[nm_y, nm_x]=hnet;
            df_dffd[nm_yr, nm_x]=1-hnet;
         }
         # by dependent xch part
         nm_x=join(".", c("d", "x", reac));
         if (nm_x %in% nm_fl) {
            tmp=1./(1.-jx_f$flnx[nm_x,1])**2;
            df_dfl[nm_y, nm_x]=tmp;
            df_dfl[nm_yr, nm_x]=tmp;
         }
         # by free xch part
         nm_x=join(".", c("f", "x", reac));
         if (nm_x %in% nm_ffx) {
            tmp=1./(1.-param[nm_x])**2;
            df_dffd[nm_y, nm_x]=tmp;
            df_dffd[nm_yr, nm_x]=tmp;
         }
      } else {
         # already derived skip it
         next;
      }
   }
   df_dff=df_dfl%*%dfl_dff+df_dffd;
   
   # measured fluxes derivation
   dfm_dff=matrix(0., length(nm_fmn), length(nm_ffn)+length(nm_ffx));
   dimnames(dfm_dff)=list(nm_fmn, c(nm_ffn, nm_ffx));
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
   dmx_df=measmat[,-dim(measmat)[2]]%*%jx_f$x_f[imeas,];
   dr_dff=c(1.,param)[ir2isc]*(dmx_df%*%df_dff);
   # part of gradient due to free fluxes
   r=measinvvar*jx_f$res;
   grad=2*(t(dr_dff)%*%r+
      t(dfm_dff)%*%(invfmnvar*jx_f$resfl));
   # part of gradient due to scale factors
   if (nb_sc > 0) {
      v=rep(0., nb_sc);
      for (isc in 1:nb_sc) {
         v[isc]=sum(r[ir2isc==(isc+1)]);
      }
      grad=c(grad, 2*v);
   }
   return(grad);
}
