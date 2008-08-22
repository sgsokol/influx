DEBUG=1;
library(bitops);
trisparse_solv=function(Al, Ac, Au, spA, b, method="dense") {
   # solve A*x=b where A=tridiag(Al,Ac,Au)+spA and b is dense
   # temporary solution by qr()
   n=length(Ac);
   if (method=="dense") {
      # fulfill a matrix
      A=matrix(0., n , n);
      for (i in 1:n) {
         if (i > 1) {
            A[i,i-1]=Al[i];
         }
         A[i,i]=Ac[i];
         if (i < n) {
            A[i,i+1]=Au[i];
         }
      }
      A=A+as.matrix(spA);
#cat("trisparse: A,b\n");
#print(A);
#print(b);
      if (DEBUG) {
         library(MASS);
         write.matrix(A,file="dbg_Acumo.txt",sep="\t");
         write.matrix(Ac,file="dbg_Ac_cumo.txt",sep="\t");
      }
      x=solve(A,b);
      return(x);
   } else if (method=="sparse") {
      # sparse
      # fulfill a matrix
      require(Matrix);
      A=Matrix(0., n , n);
      for (i in 1:n) {
         if (i > 1) {
            A[i,i-1]=Al[i];
         }
         A[i,i]=Ac[i];
         if (i < n) {
            A[i,i+1]=Au[i];
         }
      }
      A=A+spA;
      x=solve(A,b);
#q=qr(A);
#print("qr");
#print(q@V);
#print(q@R);
#print(q@p);
      return(x);
   } else if (method=="smw") {
      # Sherman-Morrison-Woodbury for low rank matrix modification
      require(matrid, lib.loc="/home/sokol/R/lib");
      atri=new("matrid", Al, Ac, Au);
      # extract non zero columns from spA
      e=integer(0);
      m=0;
      if (class(spA)=="Matrix") {
         for (j in 1:n) {
            if (spA@p[j+1] > spA@p[j]) {
               e=c(e,j);
               m=m+1;
            }
         }
         s=as.matrix(spA[,e]);
      } else {
         # spA is treated as a dense matrix
         for (j in 1:n) {
            if (any(spA[,j])) {
               e=c(e,j);
               m=m+1;
            }
         }
         s=as.matrix(spA[,e]);
      }
      A=new("matridm", atri, s, e);
      x=qr.solve(A,b);
      return(x);
   } else {
      stop(paste("Unknown method '", method, "'", sep=""));
   }
#print("A");
#print(A);
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
cumo_resid=function(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc) {
#cat("resid: \n")
#print(no_f);
#print(no_w);
#print(param);
#print(p2bfl);
   # find x for all weights
   lres=param2fl_x(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc)
#print(imeas);
   # find simulated scaled measure vector scale*(measmat*x)
   simvec=c(1.,param)[ir2isc]*(measmat%*%c(lres$x[imeas],1.));
   # diff between simulated and measured
   return(list(res=(simvec-measvec), flcnx=lres$flcnx));
}
cumo_cost=function(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn) {
#cat("cost: ");
#cat("list no_f\n")
#print(no_f);
#print(no_w);
#cat("param\n");
#print(param);
    resl=cumo_resid(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc);
   res=resl$res;
   flcnx=resl$flcnx;
   # flux residuals
   resfl=flcnx[ifmn]-fmn;
   return(sum(res*res*measinvvar)+sum(resfl*resfl*invfmnvar));
}
cumo_grad=function(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn) {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: forward finite differences f(x+h)-f(x)/h
   # x+h is taken as (1+10**-7)*x
   fact=1.e-7;
   grad=param; # make place for gradient
   # f(x)
   f=cumo_cost(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
   for (i in 1:length(param)) {
      x=param[i];
      h=x*fact;
      param[i]=x+h;
      if (param[i]==x) {
         # we are too close to zero here
         param[i]=fact;
      }
      fh=cumo_cost(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
      # restore modified param
      param[i]=x;
      grad[i]=(fh-f)/h;
   }
   return(grad);
}
param2fl_x=function(param, no_f, no_w, invAfl, p2bfl, bp, fc, imeas, measmat, measvec, ir2isc) {
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
      library(MASS);
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
      A=fwrv2Acumo(fwrv, iw);
      b=fwrv_x2bcumo(fwrv, x, iw);
      # solve the system A*x=b;
      x=c(x,trisparse_solv(A$Al, A$Ac, A$Au, A$spA, b));
      if (DEBUG && iw==1) {
         library(MASS);
         write.matrix(cbind(A$Al, A$Ac, A$Au, A$spA), file="dbg_cumoA.txt", sep="\t");
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
   tbl=matrix(0,0,3);
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
