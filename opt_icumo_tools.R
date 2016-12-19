#suppressPackageStartupMessages(require(deSolve))
#suppressPackageStartupMessages(require(Matrix))

icumo_resid=function(param, cjac, labargs) {
   # claculates residual vector of labeling propagation corresponding to param
   #cat("icumo_resid: param=", param, ", cjac=", cjac, "\n")
#cjac=F # to remove in the final version
   # from labargs to local vars
   for (item in ls(labargs)) {
      assign(item, get(item, env=labargs))
   }
   
   nb_w=length(spa)
   sqm=measurements$dev$labeled
   sqf=measurements$dev$flux
   sqp=measurements$dev$pool
#browser()
   measvecti=measurements$vec$kin
   nb_ti=nb_f$ti
   nb_sc=nb_f$nb_sc
   nb_sc_tot=nb_f$nb_sc_tot
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   
   # find simulated cumomers
   lres=lab_sim(param, cjac, labargs)
   if (!is.null(lres$err) && lres$err) {
      return(list(err=1, mes=lres$mes))
   }
   nb_lab=nb_f$nb_meas*(nb_ti-1L)
   nb_lab_tot=sum(nb_lab)
   # list of indexes in residual vector by iexp
   ir=lapply(seq_len(nb_exp), function(iexp) sum(nb_lab[seq_len(iexp-1)])+seq_len(nb_lab[[iexp]]))
   if (is.null(jx_f$jacobian)) {
      # init variables in jx_f
      jx_f$jacobian=matrix(0., nrow=nb_lab_tot+nb_f$nb_fmn+nb_f$nb_poolm, ncol=nb_ff+nb_sc_tot+nb_poolf)
      dimnames(jx_f$jacobian)=list(nm$resid, nm$par)
      # constant part of jacobian
      jx_f$simlab=jx_f$ureslab=jx_f$reslab=vector("list", nb_exp)
      jx_f$udr_dp=jx_f$jacobian
   }
   for (iexp in seq_len(nb_exp)) {
      stopifnot(length(tifull[[iexp]])!=length(tifull2[[iexp]]))
      # scale simulated measurements scale*(usm)
      if (!noscale && nb_sc[[iexp]] > 0) {
         vsc=c(1.,param)[ir2isc[[iexp]]]
         jx_f$simlab[[iexp]]=jx_f$usm[[iexp]]*vsc
      } else {
         jx_f$simlab[[iexp]]=jx_f$usm[[iexp]]
      }
      rownames(jx_f$simlab[[iexp]])=rownames(jx_f$usm[[iexp]])=nm$meas[[iexp]]
      # diff between simulated and measured
      #inna=which(!is.na(measvecti)) # for removing NA measurements
      pool[nm$poolf]=param[nm$poolf]
      if (!is.null(measvecti[[iexp]])) {
         jx_f$ureslab[[iexp]]=jx_f$simlab[[iexp]]-measvecti[[iexp]] # unreduced labeled part
         jx_f$reslab[[iexp]]=jx_f$ureslab[[iexp]]/sqm[[iexp]]
      }
#browser()
      if (cjac) {
         # calculate jacobian
         # first, its labeling part (unscaled and unreduced)
         # next, flux part and measured pools part
         # for external usage
         jx_f$udr_dp[ir[[iexp]],]=jx_f$dux_dp[[iexp]]
         if (!noscale && nb_sc[[iexp]] > 0) {
            # scale it
            jx_f$jacobian[ir[[iexp]],]=vsc*jx_f$dux_dp[[iexp]]
            # add scale part of jacobian (the sparsity pattern doesn't change)
            jx_f$dr_dsc[[iexp]][nb_f$is2mti[[iexp]]]=jx_f$reslab[[iexp]][is2mti[[iexp]]][,1]
            jx_f$jacobian[ir[[iexp]], nb_ff+seq_len(nb_sc[[iexp]])]=jx_f$dr_dsc[[iexp]]
         } else {
            jx_f$jacobian[ir[[iexp]],]=jx_f$dux_dp[[iexp]]
         }
      }
   }
   if (cjac) {
      jx_f$jacobian[nb_lab_tot+seq_len(nb_f$nb_fmn),]=dufm_dp
      jx_f$jacobian[nb_lab_tot+nb_f$nb_fmn+seq_len(nb_f$nb_poolm),]=dupm_dp
      # reduce it
      jx_f$jacobian[]=with(measurements$dev, jx_f$jacobian/c(unlist(kin), flux, pool))
      # for later usage
      jx_f$dr_dff=jx_f$jacobian[,seq_len(nb_ff),drop=FALSE]
   }
   jx_f$simfmn=lres$lf$fallnx[nm$fmn]
   jx_f$uresflu=jx_f$simfmn-measurements$vec$flux
   jx_f$simpool=(measurements$mat$pool%*%pool)[,1]
   jx_f$urespool=jx_f$simpool - measurements$vec$pool
   jx_f$resflu=jx_f$uresflu/sqf
   jx_f$respool=jx_f$urespool/sqp
   jx_f$res=c(unlist(jx_f$reslab), jx_f$resflu, jx_f$respool)
   if (length(jx_f$res)) {
      names(jx_f$res)=nm$resid
   }
   jx_f$ures=c(unlist(jx_f$ureslab), jx_f$uresflu, jx_f$urespool)
   if (length(jx_f$ures)) {
      names(jx_f$ures)=nm$resid
   }
   return(list(res=jx_f$res, jacobian=if (cjac) jx_f$jacobian else NULL))
}

icumo_cost=function(param, labargs, resl=icumo_resid(param, cjac=FALSE, labargs)) {
   if (!is.null(resl$err) && resl$err) {
      return(NULL)
   }
   return(crossprod(resl$res))
}

fwrv2sp=function(fwrv, spAbr, incu, emu=FALSE) {
   # calculate s in A*x+s (where x is cumomer vector
   # according to conventions explained in comments to python function
   # netan2Abcumo_spr() generating spAbr
   # return s
   # 2012-03-07 sokol
   # 2014-04-11 sokol: added emu option
   # 2016-09-26 sokol: adapted for arbitrary long reactions; removed sp
   
   # construct the sources s
   # for this weight
#options(warn=2)
   nb_c=spAbr$nb_c # cumomer or fragment number (when emu==T)
   w=spAbr$w # cumomer weight
   emuw=ifelse(emu, w, 1L)
   incu=as.matrix(incu)
   nco=ncol(incu)
   
   if (nb_c == 0) {
      return(simple_triplet_zero_matrix(nb_c, 1))
   }
   l=spAbr
   nm_sx=as.character(nco)
   if (is.null(l$sx)) {
      # we need sx because of many different time sets (first and second order + possible parallel experiments))
      l$sx=list()
   }
   if (is.null(l$sx[[nm_sx]])) {
      li=list()
      li$nco=nco
      li$sxmat=simple_sparse_array(i=cbind(rep(l$bmat$i, nco), rep(l$bmat$j, nco), rep(seq_len(nco), each=length(l$bmat$i))),
         v=rep(l$bmat$v, nco), dim=c(l$bmat$nrow, l$bmat$ncol, nco))
      li$s=simple_triplet_matrix(i=rep(l$b$i, nco), j=l$b$j+emuw*rep(seq_len(nco)-1, each=length(l$b$i)), v=rep(l$b$v, nco), nrow=l$b$nrow, ncol=l$b$ncol*nco)
      l$sx[[nm_sx]]=li
   }
#browser()
   ind_b=if (emu) spAbr[["ind_b_emu"]] else spAbr[["ind_b"]]
   nprodx=ncol(ind_b)-2-emu
   prodx=incu[c(ind_b[,2+emu+seq_len(nprodx)]),]
   dim(prodx)=c(nrow(ind_b), nprodx, nco)
   l$sx[[nm_sx]]$sxmat$v=c(fwrv[ind_b[,"indf"]]*arrApply(prodx, 2, "prod"))
   l$sx[[nm_sx]]$s$v=c(arrApply(as.array(l$sx[[nm_sx]]$sxmat), 1, "sum"))
   s=l$sx[[nm_sx]]$s
   return(s)
}

param2fl_usm_eul2=function(param, cjac, labargs) {
   # translate free params (fluxes+scales+pools) to fluxes and
   # unscaled simulated measurements (usm) for label propagation.
   # tifull may be more fine grained than ti. All ti must be in tifull
   # only ti moments are reported in usm and jacobian
   
   # implicite euler scheme is used on all time points in a given weight.
   # => no possibility to add a time point during a run.
   # jacobian is directly derived form discrete scheme and not from ODE solving
   
   # branched from param2fl_usm_eul().
   # 2014-07-09 sokol
   
#browser()   
   # from labargs to local vars
   for (item in ls(labargs)) {
      assign(item, get(item, env=labargs))
   }
   nb_w=length(spa)
   nb_ti=nb_f$ti
   nb_tifu=nb_f$tifu
   # cumulated sum
   nb_rcumos=nb_f$rcumos
   nbc_cumos=c(0L, cumsum(nb_rcumos))
   # calculate all fluxes from free fluxes
   fgr=numeric(nb_f$nb_fgr)
   names(fgr)=nm$nm_fgr
   fgr[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   lf=param2fl(param, labargs)
   fwrv=lf$fwrv
   nb_fwrv=length(lf$fwrv)
   nb_xi=nb_f$xi
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   nb_fgr=nb_f$nb_fgr
   nb_meas=nb_f$nb_meas
   nb_sc=nb_f$nb_sc
   # fullfill pool with free pools
   if (nb_poolf > 0) {
      pool[nm$poolf]=param[nm$poolf]
   }
   # prepare pool vectors
   # vm has the same length as the full label vector
   vm=pool[nb_f$ip2ircumo]
   if (cjac) {
      #cat("param2fl_usm_eul2: recalculate jacobian\n")
      jx_f$df_dffp=df_dffp(param, lf$flnx, nb_f, nm_list)
   }
   Alit=lapply(seq_len(nb_w), function(iw) -fwrv2Abr(fwrv, spa[[iw]], x1, nm$x[nbc_x[iw]+seq_len(nb_x[iw])], getb=F,  emu=emu)$A$triplet())
   # prepare place for (diag(vm)/dt-a)^-1
   ali_w=list()
   for (iexp in seq_len(nb_exp)) {
      if (nb_ti[iexp] < 2) {
         return(list(err=1, mes="Number of time points is less than 2"))
      }
      if (!all(ti[[iexp]] %in% tifull[[iexp]])) {
         return(list(err=1, mes="Not all time moments in ti are present in tifull vector"))
      }
      
      # prepare vectors at t1=0 with zero labeling
      # incu, xi is supposed to be in [0; 1]
      x1=c(1., xi[,iexp], rep(0., nbc_x[nb_w+1])) # later set m+0 to 1 in x1
      names(x1)=c("one", nm$inp, nm$x)
      # prepare time vectors
      dt=diff(tifull[[iexp]])
      stopifnot(all(dt > 0.))
      ntico=nb_tifu[[iexp]]-1L # number of time columns
      idt=seq_len(ntico)
      itifu=seq_len(nb_tifu[[iexp]])
      invdt=1./dt
      
      nb_mcol=ncol(measmat[[iexp]])
      # prepare ponderation with actual metab pools
      pwe[[iexp]][ipwe[[iexp]]]=pool[ip2ipwe[[iexp]]]
      spwe=tapply(pwe[[iexp]], pool_factor[[iexp]], sum)
      spwe=1./as.numeric(spwe[nm$measmat[[iexp]]])
      pwe[[iexp]]=pwe[[iexp]]*spwe
      
#browser()
      xsim=matrix(x1, nrow=length(x1), ncol=ntico)
      xsim[1+seq_len(nb_xi),]=xi[,iexp] # set input label profile
      
      dimnames(xsim)=list(names(x1), tifull[[iexp]][-1L])
      if (cjac) {
         #cat("param2fl_usm_eul2: recalculate jacobian\n")
         xpf=double(nbc_x[nb_w+1L]*(nb_ff+nb_poolf)*ntico)
         dim(xpf)=c(nbc_x[nb_w+1L], nb_ff+nb_poolf, ntico)
         if (length(ijpwef[[iexp]])) {
            dpwe=-pwe[[iexp]]*spwe
            dpwe[-ipwe[[iexp]]]=0.
            dpwe=(dpwe+dp_ones[[iexp]]*spwe)
         }
      }
      # prepare data for iadt array
      dtr=as.character(round(dt, 6L))
      dtru=unique(dtr)
      # just update already inversed matricies from previous iexp
      for (iw in seq_len(nb_w)) {
         emuw=ifelse(emu, iw, 1L)
         nb_c=spa[[iw]]$nb_c
#browser()
         ixw=nbc_x[iw]+seq_len(nb_x[iw])
         inxw=(1L+nb_xi)+ixw
         nb_row=nb_c*emuw
         inrow=(1L+nb_xi+nbc_x[iw])+seq_len(nb_row)
         imw=nbc_x[iw]+seq_len(nb_row) # mass index in x (all but last)
         vmw=vm[nbc_cumos[iw]+seq_len(nb_c)]
         vmw=rep(vmw, emuw)
         dim(vmw)=c(nb_c, emuw)
         # prepare (diag(vm)/dt-a)^-1
         am=Alit[[iw]]
         if (iexp == 1) {
            ali_w[[iw]]=lapply(invdt[pmatch(dtru, dtr)], function(dti) {
               am$v[spa[[iw]]$iadiag]=vmw[,1L]*dti+am$v[spa[[iw]]$iadiag]
               asp=Rmumps$new(am)
               return(asp@.xData[[".pointer"]])
            })
            names(ali_w[[iw]])=dtru
         } else {
            dt_add=setdiff(dtru, names(ali_w[[iw]]))
            dt_rm=setdiff(names(ali_w[[iw]]), dtru)
            if (length(dt_rm)) {
               ali_w[[iw]][dt_rm]=NULL
            }
            if (length(dt_add)) {
               nm_tmp=names(ali_w[[iw]])
               ali_w[[iw]]=append(ali_w[[iw]], lapply(invdt[pmatch(dt_add, dtr)], function(dti) {
                  am$v[spa[[iw]]$iadiag]=vmw[,1L]*dti+am$v[spa[[iw]]$iadiag]
                  asp=Rmumps$new(am)
                  return(asp@.xData[[".pointer"]])
               }))
               names(ali_w[[iw]])=c(nm_tmp, dt_add)
            }
         }
         ilua=pmatch(dtr, dtru, dup=T)
         if (emu) {
            # for the first time point, set m+0 to 1 in x1
            x1[(1L+nb_xi+nbc_x[iw])+seq_len(nb_c)]=1.
            xsim[inxw,1L]=x1[inxw]
            imwl=nbc_x[iw]+nb_row+seq_len(nb_c) # the last mass index in x
         }
         # source terms
         st=as.matrix(fwrv2sp(fwrv, spa[[iw]], xsim, emu=emu))
         # calculate labeling for all time points
         xw1=x1[inrow]
         dim(xw1)=c(nb_c, emuw)
         #xsim[inxw,]=vapply(idt, function(idtr) {
         #   xw2=lusolve(iadt[[dtr]], (xw1+st[,idtr]))
         #   xw1[] <<- xw2
         #   if (emu) {
         #      return(c(xw2, 1.-rowSums(xw2)))
         #   } else {
         #      return(xw2)
         #   }
         #}, x1[inxw])
         #xw2=vapply(idt, function(idtr) {
         #   xw2=solve(ali[[ilua[idtr]]], vmw*xw1/dt[idtr]+st[,idtr])
         #   xw1[] <<- xw2
         #}, xw1)
         dim(st)=c(nb_c, emuw, ntico)
         solve_ieu(invdt, xw1, vmw, ali_w[[iw]], st, ilua)
         #dim(xw2)=c(nb_c, emuw, ntico)
         xsim[inrow,]=st #xw2
         if (emu)
            xsim[1L+nb_xi+imwl,]=1.-arrApply(st, 2, "sum")
         if (cjac) {
#browser()
            # prepare jacobian ff, pf
            # rhs for all time points on this weight
            # parts before b_x%*%...
            #Rprof(file="fx2jr.Rprof", append=TRUE)
            xpfw=double(nb_row*(nb_ff+nb_poolf)*ntico)
            dim(xpfw)=c(nb_row, nb_ff+nb_poolf, ntico)
            xpf1=double(nb_c*emuw*(nb_ff+nb_poolf))
            dim(xpf1)=c(nb_c, emuw*(nb_ff+nb_poolf))
            sfpw=double(nb_row*ntico*nb_poolf)
            dim(sfpw)=c(nb_c, emuw, ntico, nb_poolf)
            sj=fx2jr(fwrv, spa[[iw]], nb_f, xsim)
            #Rprof(NULL)
            # ff part
            if (nb_ff+nb_fgr > 0) {
               tmp=sj$j_rhsw%stm%jx_f$df_dffp
               dim(tmp)=c(nb_row, ntico, nb_ff+nb_fgr)
               xpfw[,seq_len(nb_ff+nb_fgr),]=-aperm(tmp, c(1L, 3L, 2L))
            }
            # poolf part
            if (nb_poolf > 0) {
               i2x=nb_f$ipf2ircumo[[iexp]][[iw]]
               xw2=(cbind(x1[inrow], xsim[inrow, -ntico, drop=FALSE])-xsim[inrow, , drop=FALSE])%mrv%invdt
               dim(xw2)=c(nb_c, emuw, ntico)
               #dim(xw2)=c(nb_c, emuw*ntico)
               #tmp=at%stm%xw2+c(st)
               #dim(tmp)=c(nb_c, emuw, ntico)
               #sfpw[i2x]=tmp[i2x[,-4L]]
               sfpw[i2x]=xw2[i2x[,-4L]]
               xpfw[,nb_ff+seq_len(nb_poolf),]=(if (nb_fgr > 0) xpfw[,nb_ff+seq_len(nb_poolf,)] else 0.) + aperm(sfpw, c(1L, 2L, 4L, 3L))
            }
            # add lighter xpf (b_x%*%...)
            if (iw > 1L) {
#browser()
               #ir=seq_len(nb_row)
               #ic=seq_len(nbc_x[iw])
               #pti=proc.time()
               #for (idtr in idt) {
               #   xpfw[,,idtr]=xpfw[,,idtr]-sj$b_x[ir+(idtr-1L)*nb_row,]%stm%xpf[ic,,idtr]
               #}
               mult_bxxc(xpfw, sj$b_x, xpf)
               #cat("iw=", iw, "; len(x)=", length(sj$b_x[ir,]@x), "\n", sep="")
               #print(proc.time()-pti)
            }
            # solve the system at each time point
#browser()
            #xpfw=vapply(idt, function(idtr) {
            #   #dtr=as.character(round(dt[idtr], 10L))
            #   xpf1[] <<- lusolve(ali[[pmatch(dtr[idtr], dtru)]], c(xpfw[,,idtr])+xpf1) # +xpf1 makes that it is a matrix of a suitable size            xpf1[]=xpfw[,idtr,]
            #}, xpf1)
            #xpfw[]=vapply(idt, function(idtr) {
            #   xpf1[] <<- solve(ali[[ilua[idtr]]], c(xpfw[,,idtr])+c(vmw)*xpf1/dt[idtr]) # +xpf1 makes that it is a matrix of a suitable size            xpf1[]=xpfw[,idtr,]
            #}, xpf1)
            #dim(xpfw)=c(nb_c, emuw, nb_ff+nb_poolf, ntico)
            vmw=rep(vmw, nb_ff+nb_poolf)
            dim(vmw)=c(nb_c, emuw*(nb_ff+nb_poolf))
            dim(xpfw)=c(nb_c, emuw*(nb_ff+nb_poolf), ntico)
            solve_ieu(invdt, xpf1, vmw, ali_w[[iw]], xpfw, ilua)
            dim(xpfw)=c(nb_c, emuw, nb_ff+nb_poolf, ntico)
            #xpfw[]=aperm(xpfw, c(1L, 2L, 4L, 3L))
            xpf[imw,,]=xpfw
            if (emu) {
               # treat the last weight
               xpf[imwl,,]=-arrApply(xpfw, 2, "sum")
            }
         }
      } # iw loop
      #jx_f$xsim=xsim # for debugging only
      # get ti moments corresponding to measurements
      isel=itifu[tifull[[iexp]] %in% ti[[iexp]]][-1L]-1L
      ntise=length(isel)
      xsimf=xsim[-seq_len(1L+nb_xi),, drop=FALSE] # full simulation (in time)
      xsim=xsimf[,isel,drop=FALSE]
      # usm
      mx=measmat[[iexp]]%stm%(if (nrow(xsim) == nb_mcol) xsimf else xsimf[nm$rcumo_in_cumo,,drop=FALSE])+memaone[[iexp]]
      if (length(ipooled[[iexp]]) > 1L) {
         usmf=as.matrix(meas2sum[[iexp]]%stm%(pwe[[iexp]]*mx)) # full simulated measurements (in time)
      } else {
         usmf=mx
      }
      usm=usmf[,isel,drop=FALSE]
      if (cjac) {
         #jx_f$xpf=xpf # for debugging only
         xpf=aperm(xpf[,,isel,drop=FALSE], c(1L, 3L, 2L))
         dim(xpf)=c(nbc_x[nb_w+1L], ntise*(nb_ff+nb_poolf))
         # scale part of jacobian
         if (nb_f$nb_sc_tot != 0L) {
            stop_mes("nb_sc != 0L is not implemented as meaningless for dynamic labeling. Use --noscale option.", fcerr)
         }
         dux_dp=measmat[[iexp]]%stm%xpf
         dim(xpf)=c(nbc_x[nb_w+1L], ntise, (nb_ff+nb_poolf))
         dimnames(xpf)=list(nm$x, ti[[iexp]][-1], nm$par)
         if (length(ipooled[[iexp]]) > 1L) {
            dux_dp=meas2sum[[iexp]]%stm%(pwe[[iexp]]*dux_dp) # resize
         }
         dim(dux_dp)=c(nb_meas[iexp], ntise, nb_ff+nb_poolf)
         if (length(ijpwef[[iexp]]) > 0L) {
            # derivative of pool ponderation factor
            dpw_dpf=double(nrow(measmat[[iexp]])*ntise*nb_poolf)
            dim(dpw_dpf)=c(nrow(measmat[[iexp]]), ntise, nb_poolf)
            tmp=c(mx)*dpwe
            dim(tmp)=dim(dpw_dpf)
            dpw_dpf[ijpwef[[iexp]]]=tmp[ijpwef[[iexp]]]
            dim(dpw_dpf)= c(nrow(measmat[[iexp]]), ntise*nb_poolf)
            dux_dp[,,nb_ff+seq_len(nb_poolf)]=dux_dp[,,nb_ff+seq_len(nb_poolf)]+c(meas2sum[[iexp]]%stm%dpw_dpf)
         }
         dimnames(dux_dp)=list(nm$meas[[iexp]], ti[[iexp]][-1L], nm$par)
         jx_f$dux_dp[[iexp]]=dux_dp
         jx_f$xpf[[iexp]]=xpf
      }
      dimnames(usm)=list(nm$meas[[iexp]], ti[[iexp]][-1L])
      # store usefull information in global list jx_f
      jx_f$param=param
      jx_f$usm[[iexp]]=usm
      jx_f$usmf[[iexp]]=usmf
      jx_f$xsim[[iexp]]=xsim
      jx_f$xsimf[[iexp]]=xsimf
   }
   names(jx_f$usm)=names(jx_f$usmf)=names(jx_f$xsim)=names(jx_f$xsimf)=nm$nm_exp
   return(list(usm=jx_f$usm, x=jx_f$xsim, lf=lf))
}

param2fl_usm_rich=function(param, cjac, labargs) {
   # Richardson extrapolation to get order 2 in ODE solve
   res1=param2fl_usm_eul2(param, cjac, labargs)
   if (labargs$time_order=="2") {
      # find usimvec with step h/2
      if (is.null(labargs$jx_f2)) {
         labargs$jx_f2=new.env()
      }
      nm1=c("tifull", "tifull2", "jx_f", "jx_f2", "nb_f", "nb_exp") # names to save for orders 1 and 2
      for (item in nm1) {
         assign(item, get(item, env=labargs))
      }
      # solve with step=h/2
      labargs$tifull=tifull2
      labargs$jx_f=jx_f2
      labargs$nb_f$ipf2ircumo=nb_f$ipf2ircumo2
      labargs$nb_f$tifu=nb_f$tifu2
      res2=param2fl_usm_eul2(param, cjac, labargs)
      # restore labargs for order=1
      labargs$tifull=tifull
      labargs$nb_f=nb_f
      labargs$jx_f2=labargs$jx_f
      labargs$jx_f=jx_f
      labargs$nb_f$ipf2ircumo=nb_f$ipf2ircumo
      labargs$nb_f$tifu=nb_f$tifu
      if (!is.null(res2$err) && res2$err) {
         return(list(err=1, mes=res2$mes))
      }
      # Richardson interpolation
      for (iexp in seq_len(nb_exp)) {
         jx_f$usm[[iexp]]=2*jx_f2$usm[[iexp]]-jx_f$usm[[iexp]]
         jx_f$xsim[[iexp]]=2*jx_f2$xsim[[iexp]]-jx_f$xsim[[iexp]]
         if (cjac) {
            jx_f$dux_dp[[iexp]]=2*jx_f2$dux_dp[[iexp]]-jx_f$dux_dp[[iexp]]
            #jx_f$dux_dp=labargs$jx_f$dux_dp
         }
         res1$usm=jx_f$usm
         res1$x=jx_f$xsim
      }
   }
   return(res1)
}
