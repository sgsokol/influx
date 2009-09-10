DEBUG=0;
library(bitops);
library(MASS); # for generalized inverse

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
         if (exists("control") && control$trace) {
            cat("switch to generalized inverse\n");
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

dfc2flcnx=function(nb_f, flnx, param, fc) {
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
cumo_resid=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo") {
#cat("resid: \n")
#print(nb_f);
#print(nb_w);
#print(param);
#print(p2bfl);
   # find x for all weights
   lres=param2fl_x(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun)
#print(imeas);
   # find simulated scaled measure vector scale*(measmat*x)
   simvec=c(1.,param)[ir2isc]*(measmat%*%c(lres$x[imeas],1.));
   # diff between simulated and measured
   return(list(res=(simvec-measvec), flcnx=lres$flcnx));
}
cumo_cost=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo") {
#cat("cost: ");
#cat("list nb_f\n")
#print(nb_f);
#print(nb_w);
#cat("param\n");
#print(param);
    resl=cumo_resid(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun);
   res=resl$res;
   flcnx=resl$flcnx;
   # flux residuals
   resfl=flcnx[ifmn]-fmn;
   fn=sum(res*res*measinvvar)+sum(resfl*resfl*invfmnvar);
   if (DEBUG) {
      write.matrix(fn, file="dbg_cost.txt", sep="\t");
   }
   return(fn);
}
cumo_grad=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo") {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: forward finite differences f(x+h)-f(x)/h
   # x+h is taken as (1+fact)*x
   fact=1.e-7;
   grad=param; # make place for gradient
   # f(x)
   f=cumo_cost(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun);
   for (i in 1:length(param)) {
      x=param[i];
      h=x*fact;
      param[i]=x+h;
      if (param[i]==x) {
         # we are too close to zero here
         param[i]=fact;
      }
      fh=cumo_cost(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun);
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
   flcnx=dfc2flcnx(nb_f, flnx, param, fc);
   fwrv=flcnx2fwrv(flcnx);
   if (DEBUG) {
      write.matrix(p2bfl%*%param[1:nb_f$nb_ff]+bp, file="dbg_bfl.txt", sep="\t");
      n=length(fwrv);
      nms=paste(nm_fwrv,c(rep("fwd", n/2),rep("rev", n/2)),sep="_");
      write.matrix(cbind(1:n,nms,fwrv), file="dbg_fwrv.txt", sep="\t");
      write.matrix(cbind(1:n,nm_flcnx,flcnx), file="dbg_flcnx.txt", sep="\t");
#cat("fwrv");
#print(fwrv);
   }
   return(list(flcnx=flcnx, fwrv=fwrv));
}

param2fl_x=function(param, nb_f, nb_w, nb_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo", fj_rhs=NULL) {
   # calculate all fluxes from free fluxes
   lf=param2fl(param, nb_f, invAfl, p2bfl, bp, fc);
   # construct the system A*x=b from fluxes
   # and find x for every weight
   # if fj_rhs is not NULL, calculate jacobian x_f
   nb_fwrv=length(lf$fwrv);
   x=numeric(0);
   x_f=matrix(0., nrow=0., ncol=nb_fwrv);
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
#fj_rhs(fl, nf, x, x_f, nx, iw, n, j_rhs)
      if (length(fj_rhs)) {
         # calculate jacobian x_f
         # first, calculate right hand side for jacobian solve
         j_rhsw=matrix(0., nx, nb_fwrv);
         res<-.Fortran(fj_rhs,
         fl=as.double(lf$fwrv),
         nf=nb_fwrv,
         x=as.double(x),
         x_f=as.double(x_f),
         nx=as.integer(nx),
         iw=as.integer(iw),
         n=as.integer(ncumow),
         j_rhs=as.matrix(j_rhsw),
         NAOK=TRUE,
         DUP=FALSE);
         if (DEBUG) {
            write.matrix(j_rhsw, file=paste("dbg_j_rhs_",iw,".txt", sep=""), sep="\t");
         }
         if (iw > 0) {
            x_f=rbind(x_f, solve(lsolv$qrA, j_rhsw));
         }
      }
      # bind vectors and matrices
      x=c(x,xw);
   }
#print(x);
   return(append(list(x=x, x_f=x_f), lf));
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
