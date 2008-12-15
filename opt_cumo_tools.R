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
      x=try(solve(A,b));
      if (inherits(x, "try-error")) {
         # matrix seems to be singular
         # switch to Moore-Penrose inverse
         if (exists("control") && control$trace) {
            cat("switch to generalized inverse\n");
         }
         x=ginv(A)%*%b;
      }
      return(x);
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

dfc2flcnx=function(no_f, flnx, param, fc) {
   # produce complete flux (net,xch)*(dep,free,constr) vector
   # from dep,free,constr
   f=numeric(0);
   if (no_f$no_fln) {
      f=c(f, flnx[1:no_f$no_fln]);
   }
   if (no_f$no_ffn) {
      f=c(f, param[1:no_f$no_ffn]);
   }
   if (no_f$no_fcn) {
      f=c(f, fc[1:no_f$no_fcn]);
   }
   if (no_f$no_flx) {
      f=c(f, flnx[(no_f$no_fln+1):no_f$no_fl]);
   }
   if (no_f$no_ffx) {
      f=c(f, param[(no_f$no_ffn+1):no_f$no_ff]);
   }
   if (no_f$no_fcx) {
      f=c(f, fc[(no_f$no_fcn+1):no_f$no_fc]);
   }
   return(f);
}
cumo_resid=function(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo") {
#cat("resid: \n")
#print(no_f);
#print(no_w);
#print(param);
#print(p2bfl);
   # find x for all weights
   lres=param2fl_x(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun)
#print(imeas);
   # find simulated scaled measure vector scale*(measmat*x)
   simvec=c(1.,param)[ir2isc]*(measmat%*%c(lres$x[imeas],1.));
   # diff between simulated and measured
   return(list(res=(simvec-measvec), flcnx=lres$flcnx));
}
cumo_cost=function(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo") {
#cat("cost: ");
#cat("list no_f\n")
#print(no_f);
#print(no_w);
#cat("param\n");
#print(param);
    resl=cumo_resid(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun);
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
cumo_grad=function(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun="fwrv2rAbcumo") {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: forward finite differences f(x+h)-f(x)/h
   # x+h is taken as (1+fact)*x
   fact=1.e-7;
   grad=param; # make place for gradient
   # f(x)
   f=cumo_cost(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun);
   for (i in 1:length(param)) {
      x=param[i];
      h=x*fact;
      param[i]=x+h;
      if (param[i]==x) {
         # we are too close to zero here
         param[i]=fact;
      }
      fh=cumo_cost(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn, fortfun);
      # restore modified param
      param[i]=x;
      grad[i]=(fh-f)/h;
   }
   return(grad);
}
param2fl_x=function(param, no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc, fortfun="fwrv2rAbcumo") {
   # claculate all fluxes from free fluxes
#cat("resid: \n")
#print(no_f);
#print(no_w);
#print(param);
#print(p2bfl);
   flnx=invAfl%*%(p2bfl%*%param[1:no_f$no_ff]+bp);
#cat("flnx");
#print(flnx);
   flcnx=dfc2flcnx(no_f, flnx, param, fc);
   fwrv=flcnx2fwrv(flcnx);
   if (DEBUG) {
      write.matrix(p2bfl%*%param[1:no_f$no_ff]+bp, file="dbg_bfl.txt", sep="\t");
      n=length(fwrv);
      nms=paste(nm_fwrv,c(rep("fwd", n/2),rep("rev", n/2)),sep="_");
      write.matrix(cbind(1:n,nms,fwrv), file="dbg_fwrv.txt", sep="\t");
      write.matrix(cbind(1:n,nm_flcnx,flcnx), file="dbg_flcnx.txt", sep="\t");
#cat("fwrv");
#print(fwrv);
   }
   # construct the system A*x=b from fluxes
   # and find x for every weight
   x=numeric(0);
   for (iw in 1:no_w) {
      nx=length(x);
      ncumow=no_cumos[iw];
      A=matrix(0.,ncumow,ncumow);
      b=double(ncumow);
      #fwrv2Abcumo(fl, nf, x, nx, iw, n, A, b)
      res<-.Fortran(fortfun,
         fl=as.double(fwrv),
         nf=length(fwrv),
         x=as.double(x),
         nx=as.integer(nx),
         iw=as.integer(iw),
         n=as.integer(ncumow),
         A=as.matrix(A),
         b=as.double(b),
         NAOK=TRUE,
         DUP=FALSE);
      # solve the system A*x=b;
if (DEBUG) {
   write.matrix(cbind(A, b=b), file=paste("dbg_cumoAb_",iw,".txt", sep=""), sep="\t");
}
      xw=trisparse_solv(A, b, iw, method="dense");
      x=c(x,xw);
if (DEBUG) {
   write.matrix(xw, file=paste("dbg_cumox_",iw,".txt", sep=""), sep="\t");
}
   }
#print(x);
   return(list(x=x, flcnx=flcnx, fwrv=fwrv));
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
print_mass=function(x) {
   # separate cumos by name and order by weight
   n=length(x);
   nm_x=names(x);
   tbl=matrix(0,0,3); # metab,icumo,value_cumo
   metabs=list();
   if (length(nm_x)!=n) {
      return();
   }
   for (i in 1:n) {
      nm=nm_x[i];
      tmp=strsplit(nm,":");
      metab=tmp[[1]][1];
      icumo=tmp[[1]][2];
      tbl=rbind(tbl,c(metab,as.integer(icumo),x[i]));
      metabs[metab]="";
   }
#cat("metabs:\n");
#print(metabs);
#cat("tbl:\n");
#print(tbl);
   # extract, order and print each metab vector
   for (metab in names(metabs)) {
      cat(paste(metab,":\n",sep=""));
      im=tbl[,1]==metab;
      d=matrix(tbl[im,],nrow=sum(im), ncol=3);
#print(d);
      o=order(as.integer(d[,2]));
      # ordered cumomer vector with #0==1 component
      vcumo=c(1,as.double(d[o,3]));
      clen=log2(length(vcumo));
      # mass vector
      mass=Tiso2mass(clen)%*%(Tcumo2iso(clen)%*%vcumo);
      dimnames(mass)=list(paste("m+", 0:clen, sep=""),NULL);
      print(mass);
   }
}
