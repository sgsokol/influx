#DEBUG=0; # 1 to enable debugging information, 0=disable
#TIMEIT=0; # 1 to enable time printing at some stages
if (length(find("TIMEIT")) && TIMEIT) {
   cat("load    : ", date(), "\n", sep="", file=fclog)
}
jx_f=list()
trisparse_solv=function(A, b, w, method="dense") {
   # solve A*x=b where A=tridiag(Al,Ac,Au)+s*e^t and b is dense
   if (method=="dense") {
      n=ncol(A)
      if (DEBUG) {
         write.matrix(A, file=paste("dbg_Acumo_d_",w,".txt", sep=""),sep="\t")
      }
      # factorize the matrix
      fA=qr(A)
      d=diag(fA$qr)
      fA$rank=sum(abs(d)>abs(d[1])*1.e-14)
      if (fA$rank < n) {
          mes=paste("trisparse_solv: Cumomer matrix ",
             n, "x", n, " at weight ", w,
             " is singular.\nUnsolvable cumomers are:\n",
             paste(fA$pivot[-(1:fA$rank)], sep="", collapse="\n"),
             "\nThe matrix is dumped in dbg_Acumo_singular.txt\n",
             sep="", collapse="")
          write.matrix(A, file="dbg_Acumo_singular.txt")
          return(list(x=NULL, fA=fA, err=1, mes=mes))
      }
      x=solve(fA,b)
      return(list(x=x, fA=fA, err=0, mes=NULL))
   } else if (method=="sparse") {
      # sparse
      x=try(solve(A,b)); # A has its factorized form
      if (inherits(x, "try-error")) {
         # find 0 rows if any
         izc=apply(A, 1, function(v)sum(abs(v))<=1.e-10)
         izf=names(which(abs(jx_f$fwrv)<1.e-7))
         if (sum(izc) || length(izf)) {
            mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
               "Zero rows in cumomer matrix A at weight ", w, ":\n",
               paste(rownames(A)[izc], collapse="\n"), "\n",
               "Zero fluxes are:\n",
               paste(izf, collapse="\n"), "\n",
               sep="")
         } else {
            mes="Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3."
         }
#browser()
         stop(mes)
      }
      if (DEBUG) {
         cat("A=", str(A), "\n", sep="", file=fclog)
      }
      return(list(x=as.matrix(x), fA=A, err=0, mes=NULL))
   } else if (method=="smw") {
      # Sherman-Morrison-Woodbury for low rank matrix modification
      require(matrid, lib.loc="/home/sokol/R/lib")
      atrim=new("matridm", A)
      if (DEBUG) {
         cat(paste("dim A at weight ", w, ":\n", sep="", file=fclog))
         print(dim(A))
         write.matrix(cbind(A,b=b),file=paste("dbg_tridmA_",w,".txt", sep=""),sep="\t")
      #   print(A)
      }
      x=qr.solve(atrim,b)
      return(x)
   } else {
      stop(paste("Unknown method '", method, "'", sep=""))
   }
}

dfcg2fallnx=function(nb_f, flnx, param, fc, fg) {
   # produce complete flux (net,xch)*(dep,free,constr,growth) vector
   # from dep,free,constr,growth
   f=numeric(0)
   if (nb_f$nb_fln) {
      f=c(f, flnx[1:nb_f$nb_fln])
   }
   if (nb_f$nb_ffn) {
      f=c(f, param[1:nb_f$nb_ffn])
   }
   if (nb_f$nb_fcn) {
      f=c(f, fc[1:nb_f$nb_fcn])
   }
   if (nb_f$nb_fgr) {
      f=c(f, fg)
   }
   if (nb_f$nb_flx) {
      f=c(f, flnx[(nb_f$nb_fln+1):nb_f$nb_fl])
   }
   if (nb_f$nb_ffx) {
      f=c(f, param[(nb_f$nb_ffn+1):nb_f$nb_ff])
   }
   if (nb_f$nb_fcx) {
      f=c(f, fc[(nb_f$nb_fcn+1):nb_f$nb_fc])
   }
   if (nb_f$nb_fgr) {
      f=c(f, rep.int(0., nb_f$nb_fgr))
   }
   return(f)
}

cumo_resid=function(param, cjac=TRUE, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, emu, pool, ipooled) {
   # find x for all weights
   if (is.diff(param, jx_f$param) || is.diff(param, jx_f$paramx, 1.e-14) ||
         (cjac && is.null(jx_f$x_f))) {
      lres=param2fl_x(param, cjac, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, spAb, emu, pool, measurements, ipooled)
      if (!is.null(lres$err) && lres$err) {
         return(list(err=1, mes=lres$mes))
      }
      jx_f=lres$jx_f

      # find simulated scaled measure vector scale*(measmat*x)
      simvec=jx_f$usimcumom*c(1.,param)[ir2isc]

      # diff between simulated and measured
#expp      pool[nm$poolf]=exp(param[nm$poolf])
      pool[nm$poolf]=param[nm$poolf]
      res=c((simvec-measurements$vec$labeled), (jx_f$fallnx[measurements$mat$flux]-measurements$vec$flux), as.numeric(measurements$mat$pool%*%pool)-measurements$vec$pool)
      names(res)=nm$resid
      if (is.null(measurements$outlier) || length(measurements$outlier)==0) {
         jx_f$ures <- res
         jx_f$res <- res/c(measurements$dev$labeled, measurements$dev$flux, measurements$dev$pool)
      } else {
         # exclude outliers
         jx_f$ures <- res[-measurements$outlier]
         jx_f$res <- (res/c(measurements$dev$labeled, measurements$dev$flux, measurements$dev$pool))[-measurements$outlier]
      }
   } # else we have everything in jx_f
   # jacobian
   if (cjac) {
      if (is.diff(param, jx_f$param) || is.null(jx_f$jacobian)) {
         # recalculate it
         jx_f=cumo_jacob(param, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl,
            bp, fc, xi, measurements, ir2isc)
         jacobian=as.matrix(jx_f$udr_dp/c(measurements$dev$labeled, measurements$dev$flux, measurements$dev$pool))
         dimnames(jacobian)=list(nm$resid, nm$par)
         
         if (!is.null(measurements$outlier) && !length(measurements$outlier)==0) {
            # exclude outliers
            jacobian=jacobian[-measurements$outlier,,drop=F]
         }
         jx_f$jacobian <- jacobian
         if (nb_f$nb_ff > 0) {
            jx_f$dr_dff <- jacobian[,1:nb_f$nb_ff,drop=F]
         } else {
            jx_f$dr_dff <- matrix(0., nrow=nrow(jacobian), ncol=0)
         }
      } else {
         jacobian=jx_f$jacobian
      }
   } else {
      jacobian=NULL
   }

   return(list(res=jx_f$res, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, x=jx_f$x,
      jacobian=jacobian, jx_f=jx_f))
}
cumo_cost=function(param, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, emu, pool, ipooled) {
   # NB! This function does not return jx_f
   resl=cumo_resid(param, cjac=FALSE, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, emu, pool, ipooled)
   if (!is.null(resl$err) && resl$err) {
      return(NULL)
   }
   res=resl$res
   iva=!is.na(res)
   vres=res[iva]
   fn=crossprod(vres)[1]
   if (DEBUG) {
      write.matrix(fn, file="dbg_cost.txt", sep="\t")
   }
   return(fn)
}

param2fl=function(param, nb_f, nm, invAfl, p2bfl, g2bfl, bp, fc) {
   # claculate all fluxes from free fluxes
   fg=numeric(nb_f$nb_fg)
   names(fg)=nm$fgr
   if (nb_f$nb_fg > 0) {
#expp      fg[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*exp(param[nm$poolf])
      fg[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   }
   flnx=as.numeric(invAfl%*%(p2bfl%*%head(param, nb_f$nb_ff)+c(bp)+g2bfl%*%fg))
   names(flnx)=nm$flnx
   fallnx=c(dfcg2fallnx(nb_f, flnx, param, fc, fg))
   names(fallnx)=nm$fallnx
   fwrv=c(fallnx2fwrv(fallnx))
   names(fwrv)=nm$fwrv
   if (DEBUG) {
      write.matrix(p2bfl%*%head(param, nb_f$nb_ff)+bp, file="dbg_bfl.txt", sep="\t")
      n=length(fwrv)
      names(fwrv)=nm_fwrv
      write.matrix(fwrv, file="dbg_fwrv.txt", sep="\t")
      write.matrix(cbind(1:n, nm_fallnx, fallnx), file="dbg_fallnx.txt", sep="\t")
   }
   return(list(fallnx=fallnx, fwrv=fwrv, flnx=flnx))
}

param2fl_x=function(param, cjac=TRUE, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl,
   g2bfl, bp, fc, xi, spAb, emu, pool, measurements, ipooled) {
   # translate free params (fluxes+scales) to fluxes and cumomers
   # or emus
   nb_w=length(spAb)
   nb_xi=length(xi)
   calcx=TRUE
   if (!is.diff(param, jx_f$paramx, tol=1.e-14) &&
      (length(jx_f$x)==(spAb[[nb_w]]$nb_cl+spAb[[nb_w]]$nb_c+nb_xi+1))) {
      calcx=FALSE
      if (cjac) {
          if (!is.null(jx_f$x_f) && !is.diff(param, jx_f$param)) {
             return(list(jx_f=jx_f, x=tail(jx_f$x, -nb_xi-1), x_f=jx_f$x_f, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx))
          } # else recalculate it here
      } else {
         # just x, fallnx, ... that are already calculated
         return(list(jx_f=jx_f, x=tail(jx_f$x, -nb_xi-1), x_f=jx_f$x_f, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx))
      }
   }
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   nb_fg=(if (nb_f$include_growth_flux) nb_poolf else 0)
   fg=numeric(nb_fg)
   names(fg)=nm$fgr
#expp   fg[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*exp(param[nm$poolf])
   fg[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   # cumulated sum
   nbc_cumos=c(0, cumsum(nb_cumos))
   # calculate all fluxes from free fluxes
   if (calcx) {
      #cat("param2fl_x: recalc x\n")
      jx_f$paramx <- param
      jx_f$lA <- list()
      lf=param2fl(param, nb_f, nm, invAfl, p2bfl, g2bfl, bp, fc)
      jx_f$fallnx <- lf$fallnx
      jx_f$fwrv <- lf$fwrv
      jx_f$flnx <- lf$flnx
      jx_f$lf <- lf
      # construct the system A*x=b from fluxes
      # and find x for every weight
      # if fj_rhs is not NULL, calculate jacobian x_f
      x=c(1, xi)
   }
   if (cjac) {
      #cat("param2fl_x: recalc jac\n")
      # derivation of fwrv fluxes by free parameters: free fluxes+concentrations
      mdf_dffp=df_dffp(param, jx_f$flnx, nb_f)
      jx_f$df_dffp <- mdf_dffp
      if (emu) {
         x_f=matrix(0., nrow=sum(nb_emus), ncol=nb_ff+nb_fg)
      } else {
         x_f=matrix(0., nrow=sum(nb_cumos), ncol=nb_ff+nb_fg)
      }
   } else {
      x_f=NULL
   }
   # prepare measurement pooling operations
   nb_meas=length(ipooled$ishort)
   nb_sc=nb_f$nb_sc
   vsc=c(1.,param)[ir2isc]
   # fullfill pool with free pools
   if (nb_poolf > 0) {
#expp      pool[nm$poolf]=exp(param[nm$poolf])
      pool[nm$poolf]=param[nm$poolf]
   }
   measmat=measurements$mat$labeled
   memaone=measurements$one$labeled
   nb_mcol=ncol(measmat)
   # ponderate measmat by relative pool concentrations for pooled measurements
   collect_pools=c()
   lapply(names(ipooled), function(nmp) {
      if (nmp=="ishort") {
         return(NULL)
      }
      metabs=strsplit(nmp, ":", fix=T)[[1]][2]
      collect_pools <<- c(collect_pools, metabs)
   })
   collect_pools=unique(collect_pools)
   pwei=list()
   dpwei=list() # derivatives for jacobian
   nm_metabs=matrix(sapply(names(pool), function(nm) {
      strsplit(nm, ":")[[1]]
   }), nrow=2)[2,]
   if (cjac && nb_poolf > 0) {
      nm_metabf=matrix(sapply(nm$poolf, function(nm) {
         strsplit(nm, ":")[[1]]
      }), nrow=2)[2,]
      # matrix of ponderation derivation
      fpw2m=Matrix(0., nrow=nrow(measmat), ncol=nb_poolf)
   }
   for (metabs in collect_pools) {
      metabv=strsplit(metabs, "+", fix=T)[[1]]
      ime=match(metabv, nm_metabs)
      vp=pool[ime]
      vs=sum(vp)
      pwei[[metabs]]=vp/vs # pool is assumed non negative, non zero vector
      # auxiliary matrix for jacobian part depending on ponderation by pools
      if (cjac && nb_poolf > 0) {
         in_pf=match(metabv, nm_metabf)
         imef=which(!is.na(in_pf))
         if (length(imef)==0) {
            next
         }
         in_pf=in_pf[!is.na(in_pf)]
         icoupl=cbind(imef, in_pf)
         vd=-vp%o%rep.int(1., nb_poolf)
         vd[,-in_pf]=0.
         vd[icoupl]=vd[icoupl]+vs
         dpwei[[metabs]]=vd/(vs*vs)
      }
   }
   # ponderation itself
   measmatp=measmat
   memaonep=memaone
   for (nmp in names(ipooled)) {
      if (nmp=="ishort") {
         next
      }
      i=ipooled[[nmp]]
      metabs=strsplit(nmp, ":", fix=T)[[1]][2]
      measmatp[i,]=measmat[i,,drop=F]*pwei[[metabs]]
      memaonep[i]=memaone[i]*pwei[[metabs]]
      # auxiliary matrix for jacobian itself
      if (cjac && nb_poolf > 0 && !is.null(dpwei[[metabs]])) {
         fpw2m[i,]=dpwei[[metabs]]
      }
   }
   mema1=measmatp

   # simulate labeling weight by weight
   ba_x=0
   for (iw in 1:nb_w) {
      nb_c=spAb[[iw]]$nb_c
      if (nb_c == 0) {
         next
      }
      if (calcx) {
         #if (iw ==1) {
         #   cat("recalc x\n")
         #}
         if (emu) {
            lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], x, nm$emu[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], emu=emu)
         } else {
            lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], x, nm$rcumo[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], emu=emu)
         }
         A=lAb$A
         b=lAb$b; # may have several columns if emu is TRUE
         #solve the system A*x=b
         lsolv=trisparse_solv(lAb$A, lAb$b, iw, method="sparse")
         jx_f$lA[[iw]] <- lsolv$fA
         if (!is.null(lsolv$err) && lsolv$err) {
            return(list(err=1, mes=lsolv$mes))
         }
         if (emu) {
            xw=c(lsolv$x, 1.-rowSums(lsolv$x))
         } else {
            xw=lsolv$x
         }
         x=c(x,xw)
      } else {
         x=head(jx_f$x, 1+nb_xi+nbc_cumos[iw+1])
      }
      if (cjac) {
         #if (iw ==1) {
         #   cat("recalc jac\n")
         #}
         # calculate jacobian x_f
         # first, calculate right hand side for jacobian calculation
         # j_rhsw, b_x from sparse matrices
         # bind cumomer vector
         j_b_x=fx2jr(jx_f$fwrv, spAb[[iw]], nb_f, x)
         j_rhsw=j_b_x$j_rhsw%*%mdf_dffp
         b_x=j_b_x$b_x
         if (iw > 1 && ba_x > 0) {
            j_rhsw=j_rhsw+b_x%*%x_f[1:ba_x,,drop=F]
            if (emu) {
               dim(j_rhsw)=c(nb_c, iw*ncol(x_f))
               xf=solve(jx_f$lA[[iw]], j_rhsw)
               dimnames(xf)=list(NULL, NULL)
               xf=array(as.numeric(xf), c(nb_c, iw, nb_ff+nb_fg))
               x_f[ba_x+(1:(iw*nb_c)),]=xf
               # m+N component
               x_f[ba_x+iw*nb_c+(1:nb_c),]= -apply(xf, c(1,3), sum)
                  #-rowSums(aperm(xf, c(1,3,2)), dims=2)
            } else {
               x_f[ba_x+(1:nb_c),]=
                  as.matrix(solve(jx_f$lA[[iw]], j_rhsw))
            }
         } else if (iw == 1) {
            x_f[(1:nb_c),]=
               as.matrix(solve(jx_f$lA[[iw]], j_rhsw))
            if (emu) {
               x_f[nb_c+(1:nb_c),]=-x_f[(1:nb_c),]
            }
         }
      } else {
         if (!is.null(jx_f$param) && is.diff(param, jx_f$param)) {
            x_f=NULL
         }
      }
      if (emu) {
         ba_x=ba_x+nb_emus[iw]
      } else {
         ba_x=ba_x+nb_c
      }
   }
   if (emu) {
      names(x)=c("one", nm$xi, nm$emu)
   } else {
      names(x)=c("one", nm$xi, nm$rcumo)
   }
   jx_f$x <- x
   x=tail(x, -nb_xi-1)
   
   # calculate unreduced and unscaled measurements
   mv=mema1%*%x+memaonep
   mvp=mv # ponderated by relative pools
   if (length(ipooled) > 0) {
      lapply(names(ipooled), function(nmpo) {
         if (nmpo=="ishort") {
            return(NULL)
         }
         po=ipooled[[nmpo]]
         mvp[po[1],] <<- colSums(mvp[po,,drop=F])
         return(NULL)
      })
      mvp=mvp[ipooled$ishort,,drop=F]
   }
   jx_f$usimcumom <- as.numeric(mvp)
#print(x)
   if (cjac) {
      # unreduced residuals derivated by scale params
      dur_dsc=matrix(0., nrow=nb_meas, ncol=nb_sc)
      jacobian=matrix(0., nb_meas, nb_ff+nb_sc+nb_poolf)
      # measurement vector before pool ponderation
      mx=measmat%*%x
      # scale part of jacobian
      if (nb_f$nb_sc > 0) {
         is2m=nb_f$is2m
         dur_dsc[is2m]=mvp[is2m[,1]]
      }
#browser()
      # free flux part of jacobian (and partially free pools if present in x_f)
      if (nb_ff+nb_fg > 0) {
         mffg=mema1%*%x_f
         if (length(ipooled) > 1) {
            lapply(names(ipooled), function(nmpo) {
               if (nmpo=="ishort") {
                  return(NULL)
               }
               po=ipooled[[nmpo]]
               mffg[po[1],] <<- colSums(mffg[po,,drop=F])
               return(NULL)
            })
            mffg=mffg[ipooled$ishort,,drop=F]
         }
      } else {
         mffg=Matrix(0., nrow=nb_meas, ncol=0)
      }
      # free pool part of jacobian
      if (nb_f$nb_poolf > 0) {
         if (length(dpwei) > 0) {
            # derivation of pool ponderation factor
            dur_dpf=as.numeric(mx)*fpw2m
            lapply(names(ipooled), function(nmpo) {
               if (nmpo=="ishort") {
                  return(NULL)
               }
               po=ipooled[[nmpo]]
               dur_dpf[po[1],] <<- colSums(dur_dpf[po,,drop=F])
               return(NULL)
            })
            mpf=dur_dpf[ipooled$ishort,,drop=F]
            # growth flux depending on free pools
            if (nb_fg > 0) {
               mpf=mpf+mffg[,nb_ff+1:nb_fg,drop=F]
               mff=mffg[,1:nb_ff]
            } else {
               mff=mffg
            }
#expp            # exponential part of poolf jacobian
#expp            mpf=as.matrix(mpf%mrv%pool[nm$poolf])
         }
      } else {
         mff=mffg
         mpf=matrix(0., nrow=nb_meas, ncol=0)
      }
      if (nb_sc > 0) {
         mff=vsc*mff
         mpf=vsc*mpf
      }
      # store usefull information in global list jx_f
      jacobian=cBind(mff, dur_dsc, mpf)
      jx_f$param <- param
      jx_f$x_f <- x_f
      jx_f$ujac <- jacobian
      jx_f$jacobian <- NULL # old jacobian is no more valid
   }
   return(append(list(x=x, x_f=x_f, jx_f=jx_f), jx_f$lf))
}

Tiso2cumo=function(len) {
   if (len<0) {
      return(FALSE)
   }
   if (len==0) {
      return(matrix(1,1,1))
   }
   # recursive call for len>1
   T=Tiso2cumo(len-1)
   return(rbind(cbind(T,T),cbind(diag(0,NROW(T)),T)))
}
Tcumo2iso=function(len) {
   if (len<0) {
      return(FALSE)
   }
   if (len==0) {
      return(matrix(1,1,1))
   }
   # recursive call for len>1
   T=Tcumo2iso(len-1)
   return(rbind(cbind(T,-T),cbind(diag(0,NROW(T)),T)))
}
Tiso2mass=function(len) {
   mass=matrix(0, len+1, 2**len)
   for (i in 0:(2**len-1)) {
      s=sumbit(i)
      mass[s+1,i+1]=1
   }
   return(mass)
}
Vcumo2iso0=function(len) {
   # coefficients of first row of matrix Tcumo2iso
   # giving the conversion to isotopomer of weight 0
   if (len<0) {
      return(FALSE)
   }
   if (len==0) {
      return(c(1))
   }
   # recursive call for len>1
   V=Vcumo2iso0(len-1)
   return(c(V,-V))
}
sumbit=function(i) {
   i=as.integer(i)
   res=0
   movi=1
   while (movi<=i) {
      res=res+(bitAnd(i,movi)>0)
      movi=movi*2
   }
   return(res)
}
cumo2mass=function(x, sep=":", emusep="+") {
   # convert cumomer or emu vector(s) to MID vector(s)
   # x may be multiple column matrix,
   # each of its column is then translated into MID column.
   # x names expected in format Metab:N, where N is an integer.
   # or Metab:N+m, where m is emu M+m weght.
   
   if (length(x)==0) {
      return(NULL)
   }
   # prepare x as matrix
   if (is.vector(x)) {
      nm_x=names(x)
      x=as.matrix(x)
   } else {
      nm_x=rownames(x)
   }
   
   # is it emu or cumomer vector
   emu=TRUE
   if (is.na(strsplit(nm_x[1], emusep, fixed=T)[[1]][2])) {
      emu=FALSE
   }
   if (emu) {
      # just take the longest fragments and return their MID
      spl=data.frame(t(vapply(strsplit(nm_x, "["%s+%sep%s+%emusep%s+%"]"), c, character(3))),
         stringsAsFactors=F)
      spl[,2]=as.integer(spl[,2])
      longest=tapply(spl[,2], list(spl[,1]), max)
      o=order(names(longest))
      longest=longest[o]
      nm_l=names(longest)
      # select the MIDs of the longest fragment
      res=x[unlist(sapply(nm_l, function(nm) which(spl[,1]==nm&spl[,2]==longest[nm]))),,drop=F]
      return(res)
   }

   # separate cumos by name and order by weight
   res=matrix(0., nrow=0, ncol=ncol(x))
   metabs=c(); # unique metab names
   spl=as.matrix(sapply(nm_x, function(s) {
      v=strsplit(s, sep, fixed=T)[[1]]
      if (length(v)==2) {
         return(v)
      } else {
         # badly formed cumomer name
         return(c(NA, NA))
      }
   }))
   i=!is.na(spl[2,])
   x=x[i,,drop=F]
   spl=spl[,i,drop=F]
   n=nrow(x)
   i=1:n
   icumo=as.integer(spl[2,])
   metabs=spl[1,]
   umetabs=union(metabs, NULL)
#cat("metabs:\n")
#print(metabs)
#cat("tbl:\n")
#print(tbl)
   # extract, order and convert each metab vector
   for (metab in umetabs) {
#      cat(paste(metab,":\n",sep=""))
      im=metabs==metab
#print(d)
      o=order(icumo[im])
      # ordered cumomer vector with #0==1 component
      vcumo=rbind(1,x[im,,drop=F][o,,drop=F])
      clen=log2(nrow(vcumo))
      # check that we have all components
      sb=sumbit(max(icumo[im]))
      if (!isTRUE(all.equal(sb, clen))) {
         next
      }
      # mass vector
      mass=as.matrix(Tiso2mass(clen)%*%(Tcumo2iso(clen)%*%vcumo))
      rownames(mass)=paste(metab, "+", 0:clen, sep="")
      res=rbind(res, mass)
   }
   return(res)
}
cumo2lab=function(x) {
   # converts cumomer vector to fraction of labeled isotopomer 1-i#0
   # separate cumos by name and order by weight
   n=length(x)
   nm_x=names(x)
   if (length(nm_x)!=n) {
      return()
   }
   metabs=c(); # unique metab names
   spl=unlist(strsplit(nm_x,":",fix=T))
   i=1:n
   icumo=as.integer(spl[2*i])
   metabs=spl[2*i-1]
   umetabs=union(metabs, NULL)
#cat("metabs:\n")
#print(metabs)
#cat("tbl:\n")
#print(tbl)
   # extract, order and convert each metab vector
   res=c()
   for (metab in umetabs) {
#      cat(paste(metab,":\n",sep=""))
      im=metabs==metab
#print(d)
      o=order(icumo[im])
      # ordered cumomer vector with #0==1 component
      vcumo=c(1,x[im][o])
      clen=log2(length(vcumo))
      # labeled fraction
      lab=1-Vcumo2iso0(clen)%*%vcumo
      names(lab)=metab
      res=c(res, lab)
   }
   return(res)
}
cumo_gradj=function(param, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, emu, pool, ipooled) {
   # calculate gradient of cost function for cumomer minimization probleme
   # method: mult jacobian by residual 2*jac*resid*invvar

   # grad=c(2*t(dr_df%*%df_dffp)*(measinvvar*rescumo)+
   #    t(drf_dff)*(invfmnvar*resfl), 2*(t(dr_dw)%*%rescumo))

   # gradient
   #grad=2*(jx_f$ures/c(measurements$dev$labeled, measurements$dev$flux, measurements$dev$pool)**2)%tmm%jx_f$udr_dp
   grad=2*as.numeric(crossprod(jx_f$res, jx_f$jacobian))
   return(grad)
}
# cost function for donlp2 solver
cumo_fn=function(p) {
   return(cumo_cost(p, jx_f, nb_f, nm, nb_rw, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, emu, pool, ipooled))
}
cumo_dfn=function(p) {
   return(cumo_gradj(p, jx_f, nb_f, nm, nb_rw, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb))
}
attr(cumo_fn, "gr")=cumo_dfn
#cumo_fn@gr=cumo_dfn
cumo_jacob=function(param, jx_f, nb_f, nm, nb_cumos,
   invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc) {
   # calculate jacobian dmeas_dparam and some annexe matrices
   # without applying invvar matrix
   # The result is stored in a global list jx_f.
   if (!all(param==jx_f$param)) {
      # Normally it must not be. The cost function must be already
      # called by this moment.
      # But for some strange reason it didn't happen.
      # So let recalculate the cost and some matricies
      #res=cumo_cost(param, jx_f, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, emu)
      stopifnot(!is.diff(param, jx_f$param))
   }
   ah=1.e-10; # a heavyside parameter to make it derivable in [-ah; ah]
   
   # store flux part of jacobian for sensitivity matrix
   jx_f$dfm_dff <- dfm_dff()
   nb_sc=nb_f$nb_sc
   nb_poolf=nb_f$nb_poolf
   nb_poolm=nb_f$nb_poolm
   if (nb_sc > 0 || nb_poolf > 0) {
      jx_f$dfm_dp <- cBind(jx_f$dfm_dff, Matrix(0, nrow=nrow(jx_f$dfm_dff), ncol=nb_sc+nb_poolf))
   } else {
      # scale and pool part is just zero
      jx_f$dfm_dp <- jx_f$dfm_dff
   }
   
   # pool part of measurements
   if (nb_poolm > 0) {
#expp      jx_f$dpm_dp <- cBind(Matrix(0., nb_poolm, nb_ff+nb_sc), measurements$mat$pool[,nm$poolf, drop=F]%mrv%exp(param[nm$poolf]))
      jx_f$dpm_dp <- cBind(Matrix(0., nb_poolm, nb_ff+nb_sc), measurements$mat$pool[,nm$poolf, drop=F])
   } else {
      jx_f$dpm_dp <- matrix(0, nrow=0, ncol=length(param))
   }
   rownames(jx_f$dpm_dp) <- rownames(measurements$mat$pool)
   
   #dmx_df=measmat%*%jx_f$x_f
   #drc_dff=c(1.,param)[ir2isc]*(dmx_df%*%jx_f$df_dffp)
   
   # store isotope part of jacobian for sensitivity matrix
   #jx_f$drc_dff <- drc_dff
   
   # add flux measures
   #jx_f$dr_dff <- rBind(drc_dff, jx_f$dfm_dff)
   
   # scale factor part
   #sm=jx_f$usimcumom; # just cumomer measure part
   #z=rep(0., NROW(sm))
   # each column is fullfilled with a part of residual vector
   # corresponding to a given scale parameter
   #if (nb_sc > 0) {
   #   jx_f$drc_sc <- apply(t(nb_ff+1+(1:nb_sc)), 2, function(isc){i=ir2isc==isc; v=z; v[i]=sm[i]; v;})
   #   jx_f$udr_dp <- cBind(jx_f$dr_dff, rBind(jx_f$drc_sc, matrix(0, nrow=nb_fmn, ncol=nb_sc)))
   #} else {
   #   jx_f$udr_dp <- jx_f$dr_dff
   #}
   #dimnames(jx_f$udr_dp)=list(names(jx_f$res), nm_par)
   jx_f$udr_dp <- rBind(jx_f$ujac, jx_f$dfm_dp, jx_f$dpm_dp)
   return(jx_f)
}

fwrv2Abr=function(fwrv, spAbr, incu, nm_rcumo, getb=T, emu=F) {
   # calculate sparse A and b (in A*x=b where x is cumomer vector)
   # from fields of the list spAbr
   # return a list(A, b)
   # 2012-02-21 sokol
   #
   # update :
   # added emu parameter.
   # In emu mode b_pre_emu, ind_fbe, ind_xe1 and ind_xe2 are lists
   # in spAbr one term per m+i weght.
   # incu is inemu vector
   # nm_rcumo is nm_emu
   # when returns b, it is a matrix with as many columns as cumomer weight.
   # emu+mi for i=0,N-1 N is the fragment weight can be calculated
   # from A*x=b.
   # emu+mN have to be calculated from 1-sum(lighter weights m+i)
   # 2012-07-16 sokol
   
   nb_c=spAbr$nb_c # cumomer or fragment number (when emu==T)
   w=spAbr$w # cumomer weight
   
   if (nb_c == 0) {
      a=Matrix(0., nb_c, nb_c)
      b=Matrix(0., nb_c, 1)
      return(list(A=a, b=b))
   }
   ind_a=spAbr$ind_a
   x=fwrv[ind_a[,"indf"]]
   i=which(ind_a[,"ir0"]==ind_a[,"ic0"])
   x[i]=-x[i] # diagonal terms are negative
   A=sparseMatrix(i=ind_a[,"ir0"]+1, j=ind_a[,"ic0"]+1, x=x, dims=c(nb_c, nb_c))
   dimnames(A)=list(head(nm_rcumo, nb_c), head(nm_rcumo, nb_c))
   
   # construct a complete b vector
   if (getb) {
      if (emu) {
         #b_pre=spAbr$b_pre_emu
         #iwc=length(b_pre)
         #for (iwe in 1:iwc) {
         #   b_pre[[iwe]]@x=fwrv[spAbr$ind_fbe[[iwe]]]*incu[spAbr$ind_xe1[[iwe]]]*incu[spAbr$ind_xe2[[iwe]]]
         #}
         #b=-vapply(b_pre, colSums, double(nb_c))
         ind_b=spAbr[["ind_b_emu"]]
         x=-fwrv[ind_b[,"indf"]]*incu[ind_b[,"indx1"]]*incu[ind_b[,"indx2"]]
         b=sparseMatrix(i=ind_b[,"irow"], j=ind_b[,"iwe"], x=x, dims=c(nb_c, w))
         dimnames(b)=list(head(nm_rcumo, nb_c), 1:w)
      } else {
         ind_b=spAbr[["ind_b"]]
         x=-fwrv[ind_b[,1]]*incu[ind_b[,2]]*incu[ind_b[,3]]
         b=sparseMatrix(i=ind_b[,"irow"], j=rep.int(1, nrow(ind_b)), x=x, dims=c(nb_c, 1))
         dimnames(b)=list(head(nm_rcumo, nb_c), 1)
      }
      return(list(A=A, b=b))
   } else {
      return(list(A=A, b=NULL))
   }
}

fx2jr=function(fwrv, spAbr, nb, incu, incup=NULL) {
   # calculate sparse j_rhs and b_x from fields of the list spAbr
   # according to conventions explained in comments
   # to python function netan2Abcumo_spr() generating spAbr
   # Return a list j_rhs, b_x
   # nb is a list of various numbers (cumomers, emus and so on)
   # 2012-02-22 sokol
   #
   # update: added emu approach
   # if emu then incu is inemu vector
   # 2012-07-18 sokol
   
   # we derivate a*x=b implicitly
   # a_f*x + a*x_f=b_f + b_xl*xl_f
   emu=is.matrix(spAbr$ind_b_emu)
   x0=c(0., incu)
   nb_c=spAbr$nb_c
   nb_fwrv=spAbr$nb_fwrv
   nb_cl=spAbr$nb_cl
   w=spAbr$w
   if (nb_c==0) {
      # no system at this weight
      return(list(j_rhsw=NULL, b_x=NULL, j_rhswp=NULL, b_xp=NULL))
   }
   
   # a_fx
   ind_a=spAbr$ind_a
   i=ind_a[,"ic0"]==ind_a[,"ir0"]
   if (emu) {
      x=tail(incu, (w+1)*nb_c)
      for (iwe in 1:w) {
         tmp=x[ind_a[,"ic0"]+(1+(iwe-1)*nb_c)]
         tmp[i]=-tmp[i]
         a_fx=sparseMatrix(i=ind_a[,"ir0"]+1, j=ind_a[,"indf"], x=tmp, dims=c(nb_c, nb_fwrv))
         if (iwe == 1) {
            res=a_fx
         } else {
            res=rBind(res, a_fx)
         }
      }
      a_fx=res
   } else {
      tmp=tail(incu, nb_c)[ind_a[,"ic0"]+1]
      tmp[i]=-tmp[i]
      a_fx=sparseMatrix(i=ind_a[,"ir0"]+1, j=ind_a[,"indf"], x=tmp, dims=c(nb_c, nb_fwrv))
   }
   
   # prepare b_f
   if (emu) {
      # NB: b is shorter than emuw (or xw) by M+N vector which is added to xw as (1-sum(lighter weights))
      ind_b=spAbr$ind_b_emu
      b_f=sparseMatrix(i=ind_b[,"irow"]+nb_c*(ind_b[,"iwe"]-1), j=ind_b[,"indf"],
         x=-incu[ind_b[,"indx1"]]*incu[ind_b[,"indx2"]],
         dims=c(nb_c*w, nb_fwrv)
      )
   } else {
      ind_b=spAbr$ind_b
      b_f=sparseMatrix(i=ind_b[,"irow"], j=ind_b[,"indf"],
         x=-incu[ind_b[,"indx1"]]*incu[ind_b[,"indx2"]],
         dims=c(nb_c, nb_fwrv)
      )
   }
   
   # prepare b_x
   if (all(dim(ind_bx) > 0)) {
      if (emu) {
         ind_bx=spAbr$ind_bx_emu
         b_x=sparseMatrix(
            i=ind_bx[,"irow"],
            j=ind_bx[,"ic1"],
            x=-fwrv[ind_bx[,"indf"]]*incu[ind_bx[,"indx"]],
            dims=c(nb_c*w, spAbr$nb_emul)
         )
      } else {
         ind_bx=spAbr$ind_bx
         b_x=sparseMatrix(
            i=ind_bx[,"irow"],
            j=ind_bx[,"ic1"],
            x=-fwrv[ind_bx[,"indf"]]*incu[ind_bx[,"indx"]],
            dims=c(nb_c, spAbr$nb_cl)
         )
      }
   } else {
      if (emu) {
         b_x=Matrix(0., nrow=nb_c, ncol=spAbr$nb_emul)
      } else {
         b_x=Matrix(0., nrow=nb_c, ncol=spAbr$nb_cl)
      }
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
      b_fp@x=incup[x2tb_f[,2]]*incu[x2tb_f[,3]]+incu[x2tb_f[,2]]*incup[x2tb_f[,3]]

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
   return(list(j_rhsw=b_f-a_fx, b_x=b_x, j_rhswp=j_rhswp, b_xp=b_xp))
}

put_inside=function(param, ui, ci) {
   # put param inside of feasible domain delimited by u%*%param >= ci
   mes=""
   ineq=as.numeric(ui%*%param)-ci
   if (all(ineq>1.e-10)) {
      # nothing to do, already inside and well inside
      return(param)
   }
   dp=ldp(as.matrix(ui), -ineq)
   if (!is.null(dp)) {
      # get new active inequalities
      ineqd=as.numeric(ui%*%(param+dp))-ci
      # check that we are at least at the border and not outside
      if (any(ineqd < -1.e-7)) {
         param=NA
         attr(param, "mes")="Inequality system is ill-conditionned. Failed to solve."
         attr(param, "err")=1
         return(param)
      }
      iact=ineqd<=1.e-10
#print(ineqd[iact])
      # solve an ldp pb to find non decreasing direction
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
      names(param)=nm_par
      if (nb_ff > 0) {
         i=abs(dpn[1:nb_ff])>=1.e-10
         if (sum(i) > 0) {
            tmp=cbind(param[1:nb_ff], param[1:nb_ff]+dpn[1:nb_ff], dpn[1:nb_ff])
            dimnames(tmp)=list(nm_par[1:nb_ff], c("outside", "inside", "delta"))
            obj2kvh(tmp[i,,drop=F], "Free parameters put inside of feasible domain")
         }
      }
      # move starting point slightly inside of feasible domain
      param=param+as.numeric(dpn)
   } else {
      param=NA
      mes="Infeasible inequalities."
      if (!is.null(rownames(ui))) {
         mes=join("\n", c(mes, rownames(ui)))
      }
      attr(param, "mes")=mes
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
   nb_c=spAbr$nb_c # cumomer or fragment number (when emu==T)
   if (gets) {
      ind_b=spAbr[["ind_b"]]
      s=fwrv[ind_b[,1]]*x[ind_b[,2]]*x[ind_b[,3]]
      s=sparseMatrix(i=ind_b[,"irow"], j=rep.int(1, nrow(ind_b)), x=s, dims=c(nb_c, 1))
   } else {
      s=NULL
   }
   # sp
   if (!is.null(xp)) {
      ind_b=spAbr[["ind_b"]]
      sp=fwrv[ind_b[,1]]*(xp[ind_b[,2]]*x[ind_b[,3]]+x[ind_b[,2]]*xp[ind_b[,3]])
      sp=sparseMatrix(i=ind_b[,"irow"], j=rep.int(1, nrow(ind_b)), x=sp, dims=c(nb_c, 1))
   } else {
      sp=NULL
   }
   return(list(s=s, sp=sp))
}
df_dffp=function(param, flnx, nb_f) {
   # derivation of fwrv by free_fluxes+poolf (and not growth fluxes neither log(poolf))
   ah=1.e-10; # a heavyside parameter to make it derivable in [-ah; ah]
   nb_fwrv=length(nm_fwrv)
   nm_par=names(param)
   i_fln=grep("d.n.", names(flnx), fixed=T)
   i_flx=grep("d.x.", names(flnx), fixed=T)
   i_ffn=grep("f.n.", nm_par, fixed=T)
   i_ffx=grep("f.x.", nm_par, fixed=T)
   if (nb_f$nb_fgr > 0) {
      i_fgn=grep("pf:", nm_par, fixed=T)
   } else {
      i_fgn=c()
   }
   i_fgx=c(); #grep("g.x.", nm_par, fixed=T) # must be always empty
   nb_ff=length(i_ffn)+length(i_ffx)
   nb_fgr=length(i_fgn)+length(i_fgx)
   df_dfl=Matrix(0., length(nm_fwrv), length(flnx))
   df_dffd=Matrix(0., length(nm_fwrv), nb_ff+nb_fgr)
   # derivation by dependent fluxes
   # net part
   net=flnx[i_fln]
   hnet=Heaviside(net)
   i=abs(net)<ah
   hnet[i]=net[i]/ah
   # xch part
   xch=flnx[i_flx]
   xch=1./(1.-xch)**2
   
   # forward fluxes
   df_dfl[cfw_fl]=c(hnet, xch)
   # reverse fluxes
   df_dfl[crv_fl]=c(hnet-1., xch)
   
   # derivation by free fluxes
   # net part
   net=param[i_ffn]
   hnet=Heaviside(net)
   i=abs(net)<ah
   hnet[i]=net[i]/ah
   # xch part
   xch=param[i_ffx]
   xch=1./(1.-xch)**2
   # forward fluxes
   df_dffd[cfw_ff]=c(hnet, xch)
   # reverse fluxes
   df_dffd[crv_ff]=c(hnet-1., xch)
   
   # derivation by growth fluxes
   # forward fluxes
   df_dffd[cfw_fg]=rep.int(1., length(i_fgn))
   # reverse fluxes
   df_dffd[crv_fg]=0.
   
   res=(df_dfl%*%dfl_dffg+df_dffd)%mrv%c(rep.int(1., nb_ff), rep(nb_f$mu, nb_fgr))
   dimnames(res)=list(nm_fwrv, names(param)[c(i_ffn, i_ffx, i_fgn, i_fgx)])
   return(res)
}
dfm_dff=function() {
   # measured fluxes derivation
   res=Matrix(0., length(nm_fmn), length(nm_ff))
   dimnames(res)=list(nm_fmn, nm_ff)
   # derivate free measured fluxes (trivial)
   i=grep("f.n.", nm_fmn, fixed=T)
   if (length(i) > 0) {
      res[i,nm_fmn[i]]=diag(length(i))
   }
   # derivate dependent measured fluxes
   i=grep("d.n.", nm_fmn, fixed=T, value=T)
   if (length(i) > 0) {
      res[i,]=dfl_dffg[i,1:length(nm_ff)]
   }
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
      inna=c()
   } else {
      tim=as.numeric(colnames(m))
      inna=which(!is.na(m))
   }
   plot(range(ti, tim), range(c(x,m[inna])), t="n", ylab="Labeling", xlab="Time", ...)
   matplot(ti, t(x), t="l", lty=1:nb_curve, col=1:nb_curve, lwd=2, add=T)
   legend("topright", legend=nm, lty=1:nb_curve, col=1:nb_curve, lwd=2, cex=0.75)
   if (!is.null(m)) {
      # plot measured points
      for (i in 1:nrow(m)) {
         inna=which(!is.na(m[i,]))
         points(tim[inna], m[i,inna], col=i, cex=0.5, t="b", lwd=0.5)
         if (nrow(m) == nb_curve) {
            # draw filled polygons between simulation and data
            polygon(c(ti,rev(tim[inna])), c(x[i,], rev(m[i,inna])), col=rgb(red=t(col2rgb(i)), alpha=31, max=255), bord=F)
         }
      }
   }
}
get_usm=function(f) {
   # return list of ti, usm from a _res.kvh file f

   # get skip and end number in the kvh
   d=kvh_get_matrix(f, c("simulated unscaled labeling measurements"))
   ti=as.numeric(colnames(d))
   o=order(rownames(d))
   return(list(ti=ti, usm=d[o,,drop=F]))
}
get_labcin=function(f, nm_meas=NULL) {
   # get labeling cinetic data form file f
   # with rows matching at the best nm_meas (if any)
   d=as.matrix(read.table(f, header=T, row.names=1, sep="\t", check=F, comment=""))
   # strip the last field (row id) in nm_meas and make it correspond
   # to rownames of d
   if (is.null(nm_meas)) {
      return(d)
   }
   nm_r=rownames(d)
   nm=sapply(nm_meas, function(s) {
      v=strsplit(s, ":", fix=T)[[1]]
      v[length(v)]=""
      pmatch(paste(v, sep=":", collapse=":"), nm_r)
   })
   nm=nm[!is.na(nm)]
   return(d[nm,])
}
get_hist=function(f, v) {
   # matrix from history field from a _res.kvh file f

   # get skip and end number in the kvh
   d=kvh_get_matrix(f, c("history", v))
   return(d)
}
spr2emu=function(spr, nm_incu, nm_inemu, nb) {
   # translate spAbr structure written for reduced cumo into
   # spemu structure for EMU
   # each term f_i*x_j*x_k in b is converted into a Cauchy product terms
   # f_i*Sum_p,q[(e_e(j)+Mp)*(e_e(k)+Mq)] where p+q=current emu weight (iwe)
   # emu weight (iwe) runs from 1 (+M0) to iwc (+M(iwc-1)). 
   # So for each cumoweight iwc there will be iwc columns in the final b vector
   # Input:
   # spr - cumomer sparse structures
   # nm_incu - vector of cumomer (i.e. fragment) names
   # of length 1+inp+cumo
   # nm_inemu - names of emu components ("one"+input+emu)
   # nb is list of various numbers (cumomers, emus and so on)
   
   # Returns spr structure for emu
   # 2012-07-12 sokol
   nw=length(spr)
   spemu=spr
   if (nw < 1) {
      return(spemu)
   }
   #x2tb_f=spr[[1]]$x2tb_f
   ind_b=spr[[1]]$ind_b
   ind_b_emu=cbind(ind_b, iwe=1)
   nme2iemu=1:length(nm_inemu)
   names(nme2iemu)=nm_inemu
   
   # cumo weights are allways start from 1 to weigth max
   # if there is only weight max, the lower weights must be present
   # but corresponding matrices and vectors are of dim=0

   # for iwc==1 emu+M1 are identical to cumomers
   # and the system A*x=b for emu+M0 does not change as
   # all fluxes in A and b sum to 0.
   for (iwc in 1:nw) {
      # iwc is the resulting fragment length
      sp=spr[[iwc]]
      spemu[[iwc]]=sp
      nb_c=sp$nb_c
      if (nb_c == 0) {
         next
      }
      ind_b=sp$ind_b
      nb_ind=nrow(ind_b)
      ba_e=1+nb$xiemu
      # prepare names
      nm_c1=nm_incu[ind_b[,"indx1"]]
      nm_c2=nm_incu[ind_b[,"indx2"]]
      # get fragment length for each ind_x which is product of two terms
      inot1=ind_b[,"indx2"]!=1
      flen1=iwc%rep%nb_ind
      flen1[inot1]=vapply(strsplit(nm_c1[inot1], ":"), function(v) sumbit(as.numeric(v[2])), 1)
      flen2=iwc-flen1
      onesind=(1%rep%nb_ind)
      
      for (iwe in 1:iwc) {
         # iwe runs from m+0 to m+(iwc-1) to form all masses but last for
         # the current fragment length.
         if (iwe==1) {
            # For m+0 (iwe=1) vector b is the same in cumo and emu
            ind_b_emu=ind_b
            ind_b_emu=cbind(ind_b_emu, iwe=1)
            ind_b_emu[,"indx2"]=nme2iemu[nm_incu[ind_b[,"indx2"]]%s+%"+0"]
            ind_b_emu[!inot1,"indx2"]=1
            ind_b_emu[,"indx1"]=nme2iemu[nm_incu[ind_b[,"indx1"]]%s+%"+0"]
         } else {
            onesiwe=(1%rep%iwe)
            # form cauchy product pairs to form the mass iwe-1
            # we start with m+0 for x1 and m+(iwe-1) for x2
            add1=(0:(iwe-1)) %o% onesind
            add2=((iwe-1):0) %o% onesind
            # and remove non physical indexes
            inph=add1>(onesiwe%o%flen1) | add2>(onesiwe%o%flen2)
            ind_xe1=nme2iemu[(outer(""%rep%iwe, nm_c1, paste, sep="")%s+%"+"%s+%add1)[!inph]]
            ind_xe2=nme2iemu[(outer(""%rep%iwe, nm_c2, paste, sep="")%s+%"+"%s+%add2)[!inph]]
            ind_xe2[is.na(ind_xe2)]=1 # ones stay ones
            ind_fbe=(onesiwe%o%ind_b[,"indf"])[!inph]
            ind_ice=(onesiwe%o%ind_b[,"irow"])[!inph]
            ind_b_emu=rbind(ind_b_emu, cbind(ind_fbe, ind_xe1, ind_xe2, ind_ice, iwe))
         }
      }
      spemu[[iwc]]$ind_b_emu=ind_b_emu
      # prepare b_x
      if (all(dim(spAbr[["ind_bx"]]) > 0)) {
         i=ind_b_emu[,"indx2"]!=1 # exclude from derivation plain input entries
         tmp=ind_b_emu[i,]

         # term of d/d_x1 ( is garanted to be internal, not input cumomer)
         # => indx remain in place in indx2, ind_store remain in column indx1

         # term of d/d_x2 (x2 can be an input cumomer => no derivation)
         # => indx is taken from indx1 and goes to indx2, while ind_store goes to indx1
         i=which(tmp[,"indx2"] > ba_e)
         if (length(i)) {
            ind_bx=rbind(tmp, tmp[i,c(1,3,2,4,5)])
         } else {
            ind_bx=tmp
         }
         colnames(ind_bx)=c("indf", "ic1", "indx", "irow", "iwe")
         ind_bx[,"ic1"]=ind_bx[,"ic1"]-ba_e
         ind_bx[,"irow"]=ind_bx[,"irow"]+(ind_bx[,"iwe"]-1)*nb_c
         spemu[[iwc]]$ind_bx_emu=ind_bx
      }
      spemu[[iwc]]$nb_emul=sum(head(nb$emus, iwc-1))
   }
   return(spemu)
}
opt_wrapper=function(measurements, jx_f, trace=1) {
   if (method == "BFGS") {
      control=list(maxit=500, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc, spa, emu, pool, ipooled)
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-08,
         jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc, spa, emu, pool, ipooled)
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-08,
         jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc, spa, emu, pool, ipooled)
   } else if (method == "nlsic") {
      control=list(trace=trace, btfrac=0.25, btdesc=0.75, maxit=50, errx=1.e-5,
         ci=list(report=F), history=FALSE, adaptbt=TRUE, sln=sln,
         maxstep=max(sqrt(norm2(param)), 1.)
      )
      control[names(control_ftbl)]=control_ftbl
      res=nlsic(param, cumo_resid, 
         ui, ci, control, e=NULL, eco=NULL, flsi=lsi_fun,
         jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi,
         measurements, ir2isc,
         spa, emu, pool, ipooled)
   } else if (method == "ipopt") {
      control=list(max_iter=500, print_level=trace*5)
      control[names(control_ftbl)]=control_ftbl
      tui=c(t(ui))
      eval_g=function(x, nb_f=nb_f, nm=nm_list, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, g2bfl=g2bfl, bp=bp, fc=fc, xi=xi,
         measurements=measurements, ir2isc=ir2isc, spAb=spa) {
         return(ui%*%x)
      }
      eval_jac_g=function(x, nb_f=nb_f, nm=nm_list, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, g2bfl=g2bfl, bp=bp, fc=fc, xi=xi,
         measurements=measurements, ir2isc=ir2isc, spAb=spa) {
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
         jx_f=jx_f, nb_f=nb_f, nm=nm_list, nb_cumos=nb_rcumos,
         invAfl=invAfl, p2bfl=p2bfl, g2bfl=g2bfl, bp=bp, fc=fc, xi=xi,
         measurements=measurements, ir2isc=ir2isc, spAb=spa, emu=emu, pool=pool, ipooled=ipooled)
      res$par=res$solution
      names(res$par)=nm_par
   } else {
      cat(paste("Unknown minimization method '", method, "'\\n", sep=""), file=fcerr)
      q("no", status=1)
   }
   if (is.null(res$err)) {
      res$err=0L
   }
   return(res)
}

# wrapper for Monte-Carlo simulations
mc_sim=function(i) {
   #set.seed(seeds[i])
   # random measurement generation
   if (nb_meas) {
      meas_mc=rnorm(nb_meas, simcumom, measurements$dev$labeled)
   } else {
      meas_mc=c()
   }
   if (nb_fmn) {
      fmn_mc=rnorm(nb_fmn, simfmn, measurements$dev$flux)
   } else {
      fmn_mc=c()
   }
   if (nb_poolm) {
      poolm_mc=rnorm(nb_poolm, simpool, measurements$dev$pool)
   } else {
      poolm_mc=c()
   }
   #cat("imc=", i, "\\n", sep="")
   # minimization
   measurements_mc=measurements
   measurements_mc$vec$labeled=meas_mc
   measurements_mc$vec$flux=fmn_mc
   measurements_mc$vec$pool=poolm_mc
   jx_f=list()
   #rres=cumo_resid(param, cjac=T, jx_f, nb_f, nm_list, nb_rcumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements_mc, ir2isc, spa, emu, pool, ipooled)
   #cat("mc_res=", sqrt(norm2(rres$res)), "\n", sep="")
   #jx_f=rres$jx_f
   res=opt_wrapper(measurements_mc, jx_f, trace=0)
   #save(res, file=sprintf("mc_%d.RData", i))
   if (!is.null(res$mes) && nchar(res$mes) > 0) {
      cat((if (res$err) "Error" else "Warning"), " in Monte-Carlo i=", i, ": ", res$mes, "\n", file=fcerr, sep="")
      if (res$err) {
         return(list(cost=NA, it=res$it, normp=res$normp, par=res$par, err=res$err))
      }
   }
   # return the solution
   iva=!is.na(res$res)
   vres=res$res[iva]
   return(list(cost=crossprod(vres)[1], it=res$it, normp=res$normp, par=res$par, err=res$err))
}
