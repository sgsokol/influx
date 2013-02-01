param2fl_usm_eul=function(param, cjac=TRUE, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, spAb, pool, ti, measurements, ipooled, tifull=ti) {
   # translate free params (fluxes+scales) to fluxes and
   # unscaled simulated measurements (usm) for labeling propagation.
   # tifull may be more fine grained than ti. All ti must be in tifull
   # only ti moments are reported in usm and jacobian
   
   # implicite euler scheme is used
   # jacobian is directly derived form dicrete scheme and not from ODE solving
   
   # branched from param2fl_usm() with excponential scheme.
   # 2012-06-01 sokol
   
   # recalculate or not the labeling?
   nb_w=length(spAb)
   calcx=is.diff(param, jx_f$param, tol=1.e-14) ||
      (length(jx_f$x)!=(spAb[[nb_w]]$nb_cl+spAb[[nb_w]]$nb_c+nb_xi+1))
   if (!calcx) {
      if (cjac) {
         if (!is.null(jx_f$ujaclab)) {
            return(list(usm=jx_f$usm, x=jx_f$x, ujaclab=jx_f$ujaclab, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx))
         } # else recalculate it here
      } else {
         # just x, fallnx, ... that are already calculated
         return(list(usm=jx_f$usm, x=jx_f$x, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv))
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
      cat("param2fl_usm_eul: recalculate labprop\n")
   }
   
   dt=diff(tifull)
   stopifnot(all(dt > 0.))
   
   # cumulated sum
   nbc_cumos=c(0, cumsum(nb_cumos))
   # calculate all fluxes from free fluxes
   fgr=numeric(nb_f$nb_fgr)
   names(fgr)=nm$nm_fgr
#expp   fgr[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*exp(param[nm$poolf])
   fgr[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   lf=param2fl(param, nb_f, nm, invAfl, p2bfl, g2bfl, bp, fc)
   jx_f$fallnx <<- lf$fallnx
   jx_f$fwrv <<- lf$fwrv
   jx_f$flnx <<- lf$flnx
   nb_fwrv=length(lf$fwrv)
   nb_xi=length(xi)
   nb_poolf=nb_f$nb_poolf
   nb_meas=length(ipooled$ishort)
   nb_sc=nb_f$nb_sc
   vsc=c(1.,param)[ir2isc]
   # fullfill pool with free pools
   if (nb_poolf > 0) {
#expp      pool[nm$poolf]=exp(param[nm$poolf])
      pool[nm$poolf]=param[nm$poolf]
   }
   
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
      pwei[[metabs]]=pool[ime]/vs # pool is assumed non negative, non zero vector
      # auxiliary matrix for jacobian part depending on ponderation by pools
      if (cjac && nb_poolf > 0) {
         in_pf=match(metabv, nm_metabf)
         imef=which(!is.na(in_pf))
         if (length(imef)==0) {
            next
         }
         in_pf=in_pf[!is.na(in_pf)]
         icoupl=cbind(imef, in_pf)
         vd=-vp%o%rep(1., nb_poolf)
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
   #mema1=measmatp[,-nb_mcol,drop=F]
   #memaone=measmatp[,nb_mcol]
   mema1=measmatp
   
   #irmeas_xi=irmeas+nb_xi+1
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
      usm=matrix(0., nrow=nb_meas, ncol=nb_ti-1)
      # construct the matrices invm*A in the systems pool*dx_dt=A*x+s from fluxes
      lwA=mclapply(1:nb_w, function(iw) {fwrv2Abr(lf$fwrv, spAb[[iw]], x1, nm$rcumo[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], getb=F)$A*invmw[[iw]]})
      inviadt=list()
      xsim=matrix(0., nrow=length(x1), ncol=nb_tifu)
      xsim[,1]=x1
   } else {
      xsim=jx_f$xsim
      inviadt=jx_f$inviadt
      lwA=jx_f$lwA
      usm=jx_f$usm
      
      x1=xsim[,1]
   }
   if (cjac) {
      cat("param2fl_usm_eul: recalculate jacobian\n")
      mdf_dffp=df_dffp(param, lf$flnx, nb_f)
      jx_f$df_dffp <<- mdf_dffp
      nb_ff=nb_f$nb_ff
      nb_poolf=nb_f$nb_poolf
      xff1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      xpf1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
      if (nb_poolf > 0) {
         dur_dpf=Matrix(0., nrow=nrow(measmat), ncol=nb_poolf) # we expect this matrix sparse
      }
      dur_dsc=matrix(0., nrow=nb_meas, ncol=nb_sc)
      jacobian=array(0., dim=c(nb_meas, nb_ff+nb_sc+nb_poolf, 0))
#browser()
   } else {
      jacobian=jx_f$ujaclab
   }

   if (cjac) {
      xff2=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      xpf2=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
   } else {
      xff2=NULL
      xpf2=NULL
   }
   # time run from time t1 to time t2
   iti=2
   while (iti <= nb_tifu) {
      dtr=as.character(round(dt[iti-1], 10))
      if (calcx) {
         x2=c(1., xi)
      } else {
         x2=xsim[,iti]
      }
      # prepare (I-a*dt)^-1
      if (calcx) {
         if (is.null(inviadt[[dtr]])) {
            # new inviadt
            inviadt[[dtr]] = mclapply(1:nb_w, function(iw) {
               solve(diag(nrow(lwA[[iw]]))-lwA[[iw]]*dt[iti-1])
            })
         }
      }
      subdiv=iti==2  # subdivide or no time step (when time run instability occurs)
      iw=1
      while (iw <= nb_w) {
         # prepare xw1
         icw=nbc_cumos[iw]+(1:nb_cumos[iw])
         ixw=1+nb_xi+icw
         xw1=x1[ixw]
         # prepare s
         if (calcx) {
            sw2=as.double(fwrv2sp(lf$fwrv, spAb[[iw]], x2)$s)*invmw[[iw]]
            # make a time step for xw
            xw2=as.double(inviadt[[dtr]][[iw]]%*%(xw1+dt[iti-1]*sw2))
            
            # subdivide or not time step ?
            if (!subdiv && (any(xw2 < -1.e-3) || any(xw2 > 1.001))) {
               #cat("subdivide at iti=", iti, " and iw=", iw, "\n", sep="")
               #cat("range xw1,xw2=", range(xw1), range(xw2), "\n", sep=" ")
               # first call for subdivision => restart this iw calculation
               subdiv=T
#browser()
               #break
            }
            
            # bring to interval [0; 1]
            #xw2[xw2<0.]=0.
            #xw2[xw2>1.]=1.
            x2=c(x2, xw2)
         } else {
            xw2=x2[ixw]
            if (cjac) {
               sw2=as.double(fwrv2sp(lf$fwrv, spAb[[iw]], x2)$s)*invmw[[iw]]
            }
         }
         if (cjac) {
            # prepare jacobian ff
            sj=fx2jr(lf$fwrv, spAb[[iw]], nb_f, x2)
            if (nb_ff > 0) {
               sj$jrhs=(-invmw[[iw]])*(sj$j_rhsw%*%mdf_dffp)
               if (iw > 1) {
                  # add lighter cumomers to jacobian source term
                  sj$jrhs=sj$jrhs-invmw[[iw]]*(sj$b_x%*%xff2[1:nbc_cumos[iw],])
               }
               xff2[icw,]=as.matrix(inviadt[[dtr]][[iw]]%*%(dt[iti-1]*sj$jrhs+xff1[icw,]))
            }
            # prepare jacobian poolf
            if (nb_poolf > 0) {
               spf=matrix(0., nrow=nb_cumos[iw], ncol=nb_poolf)
               i2x=nb_f$iparpf2ix[[iw]]
               spf[i2x]=-(lwA[[iw]]%*%xw2+sw2)[i2x[,1]]
               if (iw > 1) {
                  # add lighter cumomers to jacobian
                  spf=spf-sj$b_x%*%xpf2[1:nbc_cumos[iw],,drop=F]
               }
               spf=(dt[iti-1]*invmw[[iw]])*spf
               spf=xpf1[icw,]+spf
               xpf2[icw,]=as.matrix(inviadt[[dtr]][[iw]]%*%spf)
            }
         }
         iw=iw+1 # must be the last operator of while loop
      }
      subdiv=F
      if (iw <= nb_w) {
         # there were abnormal iw-loop exit via break
         # subdivide this time interval and redo the job
         # NB: potenial infinite loop
         tifull=append(tifull, tifull[iti-1]+dt[iti-1]/4., after=iti-1)
         nb_tifu=length(tifull)
         xsim=cbind(xsim, 0.)
         dt=diff(tifull)
         cat("added new time point. nb_tifu=", nb_tifu, ", dt=", dt[iti-1], ", dim(xsim)=", dim(xsim), "\n", sep=" ")
         next
      }
      it=which(tifull[iti] == ti)
      if (calcx) {
         xsim[,iti]=x2
         
         if (length(it) > 0) {
            m2=mema1%*%tail(x2, -nb_xi-1)+memaonep
            if (cjac && nb_f$nb_poolf > 0 && length(dpwei) > 0) {
               #mx=measmat%*%c(x2[irmeas_xi], 1)
               mx=measmat%*%tail(x2, -nb_xi-1)+memaone
            }
            if (length(ipooled) > 1) {
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  m2[po[1]] <<- sum(m2[po])
                  return(NULL)
               })
               m2=m2[ipooled$ishort]
            }
            usm[,it-1]=m2
         }
      } else {
         if (length(it) > 0) {
            m2=usm[,it-1]
            #mx=measmat%*%c(x2[irmeas_xi], 1.)
            mx=measmat%*%tail(x2, -nb_xi-1)+memaone
         }
      }
      if (cjac && length(it) > 0) {
#browser()
         # scale part of jacobian
         dur_dsc[]=0.
         if (nb_f$nb_sc > 0) {
            is2m=nb_f$is2m
            dur_dsc[is2m]=m2[is2m[,1]]
         }
         # free flux part of jacobian
         if (nb_ff > 0) {
            mff=mema1%*%xff2
            if (length(ipooled) > 1) {
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  mff[po[1],] <<- colSums(mff[po,,drop=F])
                  return(NULL)
               })
               mff=mff[ipooled$ishort,,drop=F]
            }
            if (nb_sc > 0) {
               mff=vsc*mff
            }
         } else {
            mff=matrix(0., nrow=nb_meas, ncol=0)
         }
         # free pool part of jacobian
         if (nb_f$nb_poolf > 0) {
            mpf=mema1%*%xpf2
            if (length(ipooled) > 1) {
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  mpf[po[1],] <<- colSums(mpf[po,,drop=F])
                  return(NULL)
               })
               mpf=mpf[ipooled$ishort,,drop=F]
            }
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
               mpf=mpf+as.matrix(dur_dpf[ipooled$ishort,,drop=F])
            }
            if (nb_sc > 0) {
               mpf=vsc*mpf
            }
         } else {
            mpf=matrix(0., nrow=nb_meas, ncol=0)
         }
         jacobian=array(c(jacobian, as.numeric(mff), dur_dsc, as.numeric(mpf)), dim=dim(jacobian)+c(0,0,1))
      }
      # prepare next step if any
      if (iti < nb_tifu) {
         x1=x2
         if (cjac) {
            xff1=xff2
            xpf1=xpf2
         }
      }
      iti=iti+1 # must be last operator of the while loop
   }
#browser()
   #cat("dim(xsim)=", dim(xsim), "len(tifull)=", length(tifull), "\n", sep=" ")
   dimnames(usm)=list(rownames(measmat)[ipooled$ishort], ti[-1])
   dimnames(xsim)=list(c("one", nm$xi, nm$rcumo), tifull)
#print(x)
   # store usefull information in global list jx_f
   jx_f$param<<-param
   jx_f$usm<<-usm
   jx_f$xsim<<-xsim
   jx_f$inviadt<<-inviadt
   jx_f$lwA<<-lwA
   if (cjac) {
      # unscaled and permuted jacobian
      # index run measures at time 1 then 2, ... for first free flux
      # then the same for second free flux and so on
      jacobian=matrix(aperm(jacobian, c(1,3,2)), ncol=length(param))
      # transform pool part from natural to log jacobian
      if (nb_poolf > 0) {
#browser()
         jacobian[,nb_ff+nb_sc+1:nb_poolf]=jacobian[,nb_ff+nb_sc+1:nb_poolf,drop=F]%mrv%pool[nm$poolf]
      }
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
   return(append(list(usm=usm, x=x, xff=xff2, xpf=xpf2, ujaclab=jx_f$ujaclab, tifull=tifull), lf))
}
icumo_resid=function(param, cjac=TRUE, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, pool, ti, tifull) {
   # claculates residual vector of labeling propagation corresponding to param
   #cat("icumo_resid: param=", param, ", cjac=", cjac, "\n")
   nb_w=length(spAb)
   sqm=measurements$dev$labeled
   sqf=measurements$dev$flux
   sqp=measurements$dev$pool
   nb_ti=length(ti)
   nb_sc=nb_f$nb_sc
   nb_poolf=nb_f$nb_poolf
   
   jacobian=NULL
   # find usimvec
   recalcx=is.diff(param, jx_f$param) ||
         (cjac && is.null(jx_f$ujaclab))
   recalcjac=cjac && (recalcx || is.diff(param, jx_f$param) ||
         is.null(jx_f$ujaclab) || is.null(jx_f$uujac))
   if (recalcx) {
      lres=param2fl_usm_eul(param, cjac, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, spAb, pool, ti, measurements, ipooled, tifull)
      if (!is.null(lres$err) && lres$err) {
         return(list(err=1, mes=lres$mes))
      }
      # scale simulated measurements scale*(usm)
      if (nb_sc) {
         simvec=jx_f$usm*c(1.,param)[ir2isc]
      } else {
         simvec=jx_f$usm
      }
      jx_f$simvec <<- simvec
      # diff between simulated and measured
      inna=which(!is.na(measvecti)) # for removing NA measurements
      if (is.null(measvecti)) {
         jx_f$res <<- NULL
      } else {
#expp         pool[nm$poolf]=exp(param[nm$poolf])
         pool[nm$poolf]=param[nm$poolf]
         jx_f$ureslab <<- simvec-measvecti
         jx_f$uresflu <<- jx_f$fallnx[nm$fmn]-measurements$vec$flux
         jx_f$urespool <<- as.numeric(measurements$mat$pool%*%pool)-measurements$vec$pool
         jx_f$reslab <<- jx_f$ureslab*sqm
         jx_f$resflu <<- jx_f$uresflu*sqf
         jx_f$respool <<- jx_f$urespool*sqp
         jx_f$res <<- c(jx_f$reslab[inna], jx_f$resflu, jx_f$respool)
         names(jx_f$res) <<- c(outer(rownames(jx_f$reslab), paste("ti=", colnames(jx_f$reslab), sep=""), paste)[inna], names(jx_f$resflu), names(jx_f$respool))
         jx_f$ures <<- c(jx_f$ureslab[inna], jx_f$uresflu, jx_f$urespool)
         names(jx_f$ures) <<- names(jx_f$res)
      }
   }
   if (recalcjac) {
      # add measured fluxes part of jacobian
      if (nb_f$nb_fmn) {
         mdfm_dff=dfm_dff()
      } else {
         mdfm_dff=Matrix(0, nrow=0, ncol=nb_ff)
      }
#browser()
      jx_f$uujac <<- as.matrix(rBind(jx_f$ujaclab[inna,,drop=F], cBind(mdfm_dff, Matrix(0., nrow=nrow(mdfm_dff), ncol=nb_sc+nb_poolf))))
      # scale jacobian by sd
      jacobian=(1./c(rep.int(sqm, nb_ti-1)[inna], sqf, sqp))*jx_f$uujac
      #colnames(jacobian)=nm$par
      if (nb_ff > 0) {
         jx_f$dr_dff <<- jx_f$uujac[,1:nb_ff,drop=F]
      } else {
         jx_f$dr_dff <<- jx_f$uujac[,0,drop=F]
      }
      jx_f$udr_dp <<- jx_f$uujac
      jx_f$jacobian <<- jacobian # must be in accordant dimention with res
   } # else we have everything in jx_f
   #if (cjac) {
   #   cat("icumo_resid: jx_f$param=", jx_f$param, "\n")
   #}
   return(list(res=jx_f$res, jacobian=jx_f$jacobian, ures=jx_f$ures, usm=jx_f$usm, simvec=jx_f$simvec, fallnx=jx_f$fallnx))
}

icumo_cost=function(param, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ipooled, ir2isc, spAb, pool, ti, tifull) {
   resl=icumo_resid(param, cjac=FALSE, nb_f, nm, nb_cumos, invAfl, p2bfl, g2bfl, bp, fc, xi, measurements, ir2isc, spAb, pool, ti, tifull)
   if (!is.null(resl$err) && resl$err) {
      return(NULL)
   }
   res=resl$res
   fn=sum(res*res)
   if (DEBUG) {
      write.matrix(fn, file="dbg_cost.txt", sep="\t")
   }
   return(fn)
}

param2fl_usm=function(param, cjac=TRUE, nb_f, nm, nb_cumos, invAfl, p2bfl,
   g2bfl, bp, fc, xi, spAb, pool, ti, measurements, ipooled, tifull=ti) {
   # translate free params (fluxes+scales) to fluxes and
   # unscaled simulated measurements (usm) for labeling propagation
   # tifull may be more fine grained than ti. All ti must be in tifull
   # only ti moments are reported in usm and jacobian
   
   # when time run instability occurs, a time interval is subdivided
   # in two: one third and two third. The first interval is redone
   # with linear scheme.
   # 2012-05-10 sokol
   
   #print(jx_f$param-param)
   
   # recalculate or not the labeling?
   nb_w=length(spAb)
   calcx=is.null(jx_f$param) ||
      is.diff(param, jx_f$param) ||
      (length(jx_f$x)!=dim(spAb[[nb_w]]$tb_x)[1]+nb_cumos[nb_w])
   if (!calcx) {
      if (cjac) {
         if (!is.null(jx_f$ujaclab)) {
            return(list(usm=jx_f$usm, x=jx_f$x, ujaclab=jx_f$ujaclab, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv, flnx=jx_f$flnx))
         } # else recalculate it here
      } else {
         # just x, fallnx, ... that are already calculated
         return(list(usm=jx_f$usm, x=jx_f$x, fallnx=jx_f$fallnx, fwrv=jx_f$fwrv))
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
   fgr=numeric(nb_f$nb_fgr)
   names(fgr)=nm$nm_fgr
#expp   fgr[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*exp(param[nm$poolf])
   fgr[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   lf=param2fl(param, nb_f, nm, invAfl, p2bfl, g2bfl, bp, fc)
   jx_f$fallnx <<- lf$fallnx
   jx_f$fwrv <<- lf$fwrv
   jx_f$flnx <<- lf$flnx
   nb_fwrv=length(lf$fwrv)
   nb_xi=length(xi)
   nb_poolf=nb_f$nb_poolf
   nb_meas=length(ipooled$ishort)
   nb_sc=nb_f$nb_sc
   vsc=c(1.,param)[ir2isc]
   # fullfill pool with free pools
   if (nb_poolf > 0) {
#expp      pool[nm$poolf]=exp(param[nm$poolf])
      pool[nm$poolf]=param[nm$poolf]
   }
   
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
      pwei[[metabs]]=pool[ime]/vs # pool is assumed non negative, non zero vector
      # auxiliary matrix for jacobian part depending on ponderation by pools
      if (cjac && nb_poolf > 0) {
         in_pf=match(metabv, nm_metabf)
         imef=which(!is.na(in_pf))
         if (length(imef)==0) {
            next
         }
         in_pf=in_pf[!is.na(in_pf)]
         icoupl=cbind(imef, in_pf)
         vd=-vp%o%rep(1., nb_poolf)
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
   #mema1=measmatp[,-nb_mcol,drop=F]
   mema1=measmatp
   #memaone=measmatp[,nb_mcol]
   
   #irmeas_xi=irmeas+nb_xi+1
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
      usm=matrix(0., nrow=nb_meas, ncol=nb_ti-1)
      # construct the matrices invm*A in the systems pool*dx_dt=A*x+s from fluxes
      lwA=lapply(1:nb_w, function(iw) {fwrv2Abr(lf$fwrv, spAb[[iw]], x1, nm$rcumo[(nbc_cumos[iw]+1):nbc_cumos[iw+1]], getb=F)$A*invmw[[iw]]})
      lwinva=lapply(1:nb_w, function(iw) {
         a=lwA[[iw]]
         ai=try(Matrix::solve(a))
         if (inherits(ai, "try-error")) {
            # find 0 rows if any
            izc=apply(as.matrix(a), 1, function(v)sum(abs(v))<=1.e-10)
            izf=names(which(abs(jx_f$fwrv)<1.e-7))
            if (sum(izc) || length(izf)) {
               mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
                  "Zero rows in cumomer matrix A at weight ", iw, ":\n",
                  paste(rownames(a)[izc], collapse="\n"), "\n",
                  "Zero fluxes are:\n",
                  paste(izf, collapse="\n"), "\n",
                  sep="")
            } else {
               mes="Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options option with small N, say 1.e-3."
            }
#browser()
            stop(mes)
         }
         return(ai)
      })
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
      mdf_dffp=df_dffp(param, lf$flnx)
      jx_f$df_dffp <<- mdf_dffp
      nb_ff=nb_f$nb_ff
      nb_poolf=nb_f$nb_poolf
      xff1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      xffp1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_ff)
      if (nb_ff > 0) {
         # jacobian ff rhs
         sj1=lapply(1:nb_w, function(iw) {
            sj=fx2jr(lf$fwrv, spAb[[iw]], nb_f, x1, xp1)
            sj$jrhs=(-invmw[[iw]])*(sj$j_rhsw%*%mdf_dffp)
            sj$jrhsp=(-invmw[[iw]])*(sj$j_rhswp%*%mdf_dffp)
            xffp1[nbc_cumos[iw]+1:nb_cumos[iw],] <<- as.matrix(sj$jrhs)
            return(sj)
         })
      }
      xpf1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
      xpfp1=matrix(0., nrow=nbc_cumos[nb_w+1], ncol=nb_poolf)
      if (nb_poolf > 0) {
         dur_dpf=Matrix(0., nrow=nrow(measmat), ncol=nb_poolf) # we expect this matrix sparse
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
   iti=2
   while (iti <= nb_tifu) {
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
         # see if matrix power can be used
         if (length(expadt)==0 || !is.null(expadt[[dtr]])) {
            usep=F
         } else {
            vdt=dt[iti-1]/as.numeric(names(expadt))
            ipow=which(abs(vdt-round(vdt)) < 1.e-7) # where is an integer power
            if (length(ipow) > 0) {
               idt=ipow[which.min(vdt[ipow])]
               mpow=vdt[idt]
               if (mpow > 8) {
                  usep=F
               } else {
                  usep=T
                  #cat("use matrix power idt=", idt, " mpow=", mpow, " iti=", iti, "\n", sep="")
               }
            } else {
               usep=F
            }
         }
         lapply(1:nb_w, function(iw) {
            if (is.null(expadt[[dtr]])) {
               expadt[[dtr]] <<- list()
            }
            if (length(expadt[[dtr]]) < iw) {
               # new exp(a*dt)
               if (usep) {
                  expadt[[dtr]][[iw]] <<- powmx(expadt[[idt]][[iw]], mpow)
               } else {
                  expadt[[dtr]][[iw]] <<- expm(as.matrix(lwA[[iw]]*dt[iti-1]))
               }
            }
         })
      }
      subdiv=iti==2  # subdivide or no time step (when time run instability occurs)
      iw=1
      while (iw <= nb_w) {
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
            if (iw == 1) {
               xw2=as.double(expm_const_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sw1, xw1))
            } else if (F && !subdiv) {
               # in the hope that the s curve is smooth enough
               xw2=as.double(expm_cub_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sw1, sw2, swp1, swp2, xw1))
            } else {
               # for a steep step subdivide and low order to preserve monotony
               xw2=as.double(expm_lin_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sw1, sw2, xw1))
            }
            # subdivide or not time step ?
            if (!subdiv && (any(xw2 < -1.e-3) || any(xw2 > 1.001))) {
               #cat("subdivide at iti=", iti, " and iw=", iw, "\n", sep="")
               #cat("range xw1,xw2=", range(xw1), range(xw2), "\n", sep=" ")
               # first call for subdivision => restart this iw calculation
               subdiv=T
#browser()
               #break
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
               sj=fx2jr(lf$fwrv, spAb[[iw]], nb_f, x2, xp2)
               sj$jrhs=(-invmw[[iw]])*(sj$j_rhsw%*%mdf_dffp)
               sj$jrhsp=(-invmw[[iw]])*(sj$j_rhswp%*%mdf_dffp)
               if (iw > 1) {
                  # add lighter cumomers to jacobian source term
                  sj$jrhs=sj$jrhs-invmw[[iw]]*(sj$b_x%*%xff2[1:nbc_cumos[iw],])
                  sj$jrhsp=sj$jrhsp-invmw[[iw]]*(sj$b_xp%*%xff2[1:nbc_cumos[iw],]+sj$b_x%*%xffp2[1:nbc_cumos[iw],])
               }
               sj2[[iw]]=sj
               # make a time step for jacobian xff
               if (F && !subdiv) {
                  xff2[icw,]=as.matrix(expm_cub_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sj1[[iw]]$jrhs, sj2[[iw]]$jrhs, sj1[[iw]]$jrhsp, sj2[[iw]]$jrhsp, xff1[icw,,drop=F]))
               } else {
                  xff2[icw,]=as.matrix(expm_lin_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], sj1[[iw]]$jrhs, sj2[[iw]]$jrhs, xff1[icw,,drop=F]))
               }
               xffp2[icw,]=as.matrix(lwA[[iw]]%*%xff2[icw,,drop=F])+as.matrix(sj2[[iw]]$jrhs)
            }
            # prepare source 2 for jacobian poolf
            if (nb_poolf > 0) {
               if (nb_ff==0 && iw > 1) {
                  sj=fx2jr(lf$fwrv, spAb[[iw]], nb_f, x2, xp2)
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
               if (F && !subdiv) {
                  xpf2[icw,]=as.matrix(expm_cub_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], spf1[[iw]]$s, spf2[[iw]]$s, spf1[[iw]]$sp, spf2[[iw]]$sp, xpf1[icw,,drop=F]))
               } else {
                  xpf2[icw,]=as.matrix(expm_lin_step(lwinva[[iw]], dt[iti-1], expadt[[dtr]][[iw]], spf1[[iw]]$s, spf2[[iw]]$s, xpf1[icw,,drop=F]))
               }
               xpfp2[icw,]=as.matrix(lwA[[iw]]%*%xpf2[icw,,drop=F])+as.matrix(spf)
            }
         }
         iw=iw+1 # must be the last operator of while loop
      }
      subdiv=F
      if (iw <= nb_w) {
         # there were abnormal iw-loop exit via break
         # subdivide this time interval and redo the job
         # NB: potenial infinite loop
         tifull=append(tifull, tifull[iti-1]+dt[iti-1]/4., after=iti-1)
         nb_tifu=length(tifull)
         xsim=cbind(xsim, 0.)
         xpsim=cbind(xpsim, 0.)
         dt=diff(tifull)
         cat("added new time point. nb_tifu=", nb_tifu, ", dt=", dt[iti-1], ", dim(xsim)=", dim(xsim), "\n", sep=" ")
         next
      }
      it=which(tifull[iti] == ti)
      if (calcx) {
         xsim[,iti]=x2
         xpsim[,iti]=xp2
         spx=append(spx, list(sp2))
         
         if (length(it) > 0) {
            m2=mema1%*%tail(x2, -nb_xi-1)+memaonep
            if (cjac && nb_f$nb_poolf > 0 && length(dpwei) > 0) {
               mx=measmat%*%tail(x2, -nb_xi-1)
            }
            if (length(ipooled) > 1) {
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  m2[po[1]] <<- sum(m2[po])
                  return(NULL)
               })
               m2=m2[ipooled$ishort]
            }
            usm[,it-1]=m2
         }
      } else {
         if (length(it) > 0) {
            m2=usm[,it-1]
            mx=measmat%*%tail(x2, -nb_xi-1)
         }
      }
      if (cjac && length(it) > 0) {
         # scale part of jacobian
         dur_dsc[]=0.
         if (nb_f$nb_sc > 0) {
            is2m=nb_f$is2m
            dur_dsc[is2m]=m2[is2m[,1]]
         }
         # free flux part of jacobian
         if (nb_ff > 0) {
            mff=mema1%*%xff2
            if (length(ipooled) > 1) {
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  mff[po[1],] <<- colSums(mff[po,,drop=F])
                  return(NULL)
               })
               mff=mff[ipooled$ishort,,drop=F]
            }
            if (nb_sc > 0) {
               mff=vsc*mff
            }
         } else {
            mff=matrix(0., nrow=nb_meas, ncol=0)
         }
         # free pool part of jacobian
         if (nb_f$nb_poolf > 0) {
            mpf=mema1%*%xpf2
            if (length(ipooled) > 1) {
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  mpf[po[1],] <<- colSums(mpf[po,,drop=F])
                  return(NULL)
               })
               mpf=mpf[ipooled$ishort,,drop=F]
            }
            if (length(dpwei) > 0) {
               # derivation of pool ponderation factor
               dur_dpf=c(mx)*fpw2m
               lapply(names(ipooled), function(nmpo) {
                  if (nmpo=="ishort") {
                     return(NULL)
                  }
                  po=ipooled[[nmpo]]
                  dur_dpf[po[1],] <<- colSums(dur_dpf[po,,drop=F])
                  return(NULL)
               })
               mpf=mpf+as.matrix(dur_dpf[ipooled$ishort,,drop=F])
            }
            if (nb_sc > 0) {
               mpf=vsc*mpf
            }
         } else {
            mpf=matrix(0., nrow=nb_meas, ncol=0)
         }
         jacobian=array(c(as.numeric(jacobian), as.numeric(mff), as.numeric(dur_dsc), as.numeric(mpf)), dim=dim(jacobian)+c(0,0,1))
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
      iti=iti+1 # must be last operator of the while loop
   }
#browser()
   #cat("dim(xsim)=", dim(xsim), "len(tifull)=", length(tifull), "\n", sep=" ")
   dimnames(usm)=list(rownames(measmat)[ipooled$ishort], ti[-1])
   dimnames(xsim)=list(c("one", nm$xi, nm$rcumo), tifull)
   dimnames(xpsim)=list(c("one", nm$xi, nm$rcumo), tifull)
#print(x)
   # store usefull information in global list jx_f
   jx_f$param<<-param
   jx_f$usm<<-usm
   jx_f$xsim<<-xsim
   jx_f$xpsim<<-xpsim
   jx_f$spx<<-spx
   jx_f$expadt<<-expadt
   jx_f$lwA<<-lwA
   jx_f$lwinva<<-lwinva
   if (cjac) {
      # unscaled and permuted jacobian
      # index run measures at time 1 then 2, ... for first free flux
      # then the same for second free flux and so on
      jacobian=matrix(aperm(jacobian, c(1,3,2)), ncol=length(param))
      # transform pool part from natural to log jacobian
      if (nb_poolf > 0) {
#browser()
         jacobian[,nb_ff+nb_sc+1:nb_poolf]=jacobian[,nb_ff+nb_sc+1:nb_poolf,drop=F]%mrv%pool[nm$poolf]
      }
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
   return(append(list(usm=usm, x=x, xff=xff2, xpf=xpf2, ujaclab=jx_f$ujaclab, tifull=tifull), lf))
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
expm_const_step=function(inva, dt, expadt, s1, x1) {
   # make a single time step for integration of ode xp=a*x+s
   # according to scheme with constant source term  and exponential matrix
   # return a vector x
   # 2012/05/10 sokol
   
   # increment
   dx=inva%*%(expadt%*%s1-s1)
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
