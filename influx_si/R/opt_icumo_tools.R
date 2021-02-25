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
   nbc_lab=c(0L, cumsum(nb_lab))
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
         #jx_f$udr_dp[ir[[iexp]],]=jx_f$dux_dp[[iexp]]
         bop(jx_f$udr_dp, c(1, nbc_lab[[iexp]], nb_lab[[iexp]]), "=", jx_f$dux_dp[[iexp]])
         if (!noscale && nb_sc[[iexp]] > 0) {
            # scale it
            #jx_f$jacobian[ir[[iexp]],]=vsc*jx_f$dux_dp[[iexp]]
            bop(jx_f$udr_dp, c(1, nbc_lab[[iexp]], nb_lab[[iexp]]), "*=", vsc)
            # add scale part of jacobian (the sparsity pattern doesn't change)
            jx_f$dr_dsc[[iexp]][nb_f$is2mti[[iexp]]]=jx_f$reslab[[iexp]][is2mti[[iexp]]][,1]
            jx_f$jacobian[ir[[iexp]], nb_ff+seq_len(nb_sc[[iexp]])]=jx_f$dr_dsc[[iexp]]
         } else {
            #jx_f$jacobian[ir[[iexp]],]=jx_f$dux_dp[[iexp]]
            bop(jx_f$jacobian, c(1, nbc_lab[[iexp]], nb_lab[[iexp]]), "=", jx_f$dux_dp[[iexp]])
         }
      }
   }
   if (cjac) {
      #jx_f$jacobian[nb_lab_tot+seq_len(nb_f$nb_fmn),]=dufm_dp
      bop(jx_f$jacobian, c(1, nb_lab_tot, nb_f$nb_fmn), "=", dufm_dp)
      #jx_f$jacobian[nb_lab_tot+nb_f$nb_fmn+seq_len(nb_f$nb_poolm),]=dupm_dp
      bop(jx_f$jacobian, c(1, nb_lab_tot+nb_f$nb_fmn, nb_f$nb_poolm), "=", dupm_dp)
      # reduce it
      #jx_f$jacobian[]=with(measurements$dev, jx_f$jacobian/c(unlist(kin), flux, pool))
      bop(jx_f$jacobian, 1, "/=", with(measurements$dev, c(unlist(kin), flux, pool)))
      # for later use
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
   if (any(ina <- is.na(resl$res))) {
      return(crossprod(resl$res[!ina]))
   } else {
      return(crossprod(resl$res))
   }
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
      #li$sxmat=simple_sparse_array(i=cbind(rep(l$bmat$i, nco), rep(l$bmat$j, nco), rep(seq_len(nco), each=length(l$bmat$i))),
      #   v=rep(l$bmat$v, nco), dim=c(l$bmat$nrow, l$bmat$ncol, nco))
      li$sxmat=simple_triplet_matrix(i=rep(l$bmat$i, nco), j=l$bmat$j+rep(seq(0, nco-1), each=length(l$bmat$i))*l$bmat$ncol, v=rep(l$bmat$v, nco), nrow=l$bmat$nrow, ncol=l$bmat$ncol*nco)
      li$s=simple_triplet_matrix(i=rep(l$b$i, nco), j=l$b$j+emuw*rep(seq_len(nco)-1, each=length(l$b$i)), v=rep(l$b$v, nco), nrow=l$b$nrow, ncol=l$b$ncol*nco)
      l$sx[[nm_sx]]=li
   }
#browser()
   ind_b=if (emu) spAbr[["ind_b_emu"]] else spAbr[["ind_b"]]
   nprodx=ncol(ind_b)-2-emu
   prodx=incu[c(ind_b[,2+emu+seq_len(nprodx)]),]
   dim(prodx)=c(nrow(ind_b), nprodx, nco)
   l$sx[[nm_sx]]$sxmat$v[]=fwrv[ind_b[,"indf"]]*arrApply(prodx, 2, "prod")
   #l$sx[[nm_sx]]$s$v[]=arrApply(as.array(l$sx[[nm_sx]]$sxmat), 1, "sum")
   l$sx[[nm_sx]]$s$v[]=slam::col_sums(l$sx[[nm_sx]]$sxmat)
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
   #cat("labargs=", format(labargs), "\n")
   for (item in ls(labargs)) {
      assign(item, get(item, env=labargs))
   }
   if (is.null(labargs$getx))
      getx=FALSE
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
   Alit=lapply(seq_len(nb_w), function(iw) -fwrv2Abr(fwrv, spa[[iw]], NULL, nm$x[nbc_x[iw]+seq_len(nb_x[iw])], getb=FALSE,  emu=emu)$A$triplet())
   dtru=unique(unlist(lapply(seq_len(nb_exp), function(iexp) as.character(round(diff(tifull[[iexp]]), 6)))))
   # prepare place for (diag(vm)/dt-a)^-1 common to all iexp
   if (is.null(labargs$ali_w)) {
      ali_w=list()
      for (iw in seq_len(nb_w)) {
         nb_c=spa[[iw]]$nb_c
         emuw=ifelse(emu, iw, 1L)
         vmw=vm[nbc_cumos[iw]+seq_len(nb_c)]
         vmw=rep(vmw, emuw)
         redim(vmw, c(nb_c, emuw))
         ali_w[[iw]]=lapply(dtru, function(dtu) {
            dti=1./as.double(dtu)
            a=Alit[[iw]]
            #a$v[spa[[iw]]$iadiag]=vmw[,1L]*dti+a$v[spa[[iw]]$iadiag]
            asp=Rmumps$new(a) # just a place-holder
            asp$set_icntl(3, 7) # 2=amf, 3=scotch, 4=pord, 5=metis
            #asp$set_icntl(400, 14) # increase by 400% working space
            #if (packageVersion("rmumps") >= "5.1.1-1")
            #   asp$set_keep(40, 1) # undocumented feature of mumps (cf. their mail on mumps-user group from 12/04/2017)
            return(asp)
         })
         names(ali_w[[iw]])=dtru
      }
   }
   ntico_max=max(sapply(seq_len(nb_exp), function(iexp) nb_tifu[[iexp]]-1L))
   nb_row_max=max(sapply(seq_len(nb_w), function(iw) {emuw=ifelse(emu, iw, 1L); spa[[iw]]$nb_c*emuw}))
   xpf=double(nbc_x[nb_w+1L]*(nb_ff+nb_poolf)*ntico_max)
   #sfpw=double(nb_row_max*ntico_max*nb_poolf)
   xpfw=double(nb_row_max*(nb_ff+nb_poolf)*ntico_max)
   #xpf1=double(nb_row_max*(nb_ff+nb_poolf))
   for (iexp in seq_len(nb_exp)) {
      if (nb_ti[iexp] < 2) {
         return(list(err=1, mes="Number of time points is less than 2"))
      }
      if (!all(ti[[iexp]] %in% tifull[[iexp]])) {
         return(list(err=1, mes="Not all time moments in ti are present in tifull vector"))
      }
      
      # prepare vectors at t1=0 with zero labeling
      # incu, xi is supposed to be in [0; 1]
      x1=c(1., xi[[iexp]][,1L], rep(0., nbc_x[nb_w+1])) # later set m+0 to 1 in x1
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
      bop(xsim, c(1, 1, nb_xi), "=", xi[[iexp]][,-1])
      
      dimnames(xsim)=list(names(x1), tifull[[iexp]][-1L])
      if (cjac) {
         #cat("param2fl_usm_eul2: recalculate jacobian\n")
         #xpf=double(nbc_x[nb_w+1L]*(nb_ff+nb_poolf)*ntico)
         #redim(xpf, c(nbc_x[nb_w+1L], nb_ff+nb_poolf, ntico))
         resize(xpf, c(nbc_x[nb_w+1L], nb_ff+nb_poolf, ntico))
         if (length(ijpwef[[iexp]])) {
            dpwe=-pwe[[iexp]]*spwe
            dpwe[-ipwe[[iexp]]]=0.
            dpwe=(dpwe+dp_ones[[iexp]]*spwe)
         }
      }
      # prepare data for iadt array
      dtr=as.character(round(dt, 6L))
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
         redim(vmw, c(nb_c, emuw))
         # prepare (diag(vm)/dt-a)^-1
         if (iexp == 1) {
            ali_w[[iw]][]=lapply(names(ali_w[[iw]]), function(dtu) {
               dti=1./as.double(dtu)
               asp=ali_w[[iw]][[dtu]]
               if (!inherits(asp, "Rcpp_Rmumps"))
                  asp=attr(asp, "asp")
               av=Alit[[iw]]$v
               av[spa[[iw]]$iadiag]=vmw[,1L]*dti+av[spa[[iw]]$iadiag]
               asp$set_mat_data(av)
               if (cjac && nb_c < 65) {
                  m=asp$inv()
                  attr(m, "asp")=asp
                  return(m)
               } else {
                  return(asp)
               }
            })
         }
         ilua=pmatch(dtr, names(ali_w[[iw]]), dup=T)
         if (emu) {
            # for the first time point, set m+0 to 1 in x1
            x1[(1L+nb_xi+nbc_x[iw])+seq_len(nb_c)]=1.
            xsim[inxw,1L]=x1[inxw]
            imwl=nbc_x[iw]+nb_row+seq_len(nb_c) # the last mass index in x
         }
         # source terms
         st=fwrv2sp(fwrv, spa[[iw]], xsim, emu=emu)
         st=as.matrix(st)
         # calculate labeling for all time points
         xw1=x1[inrow]
         redim(xw1, c(nb_c, emuw))
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
         redim(st, c(nb_c, emuw, ntico))
         solve_ieu(invdt, xw1, vmw, ali_w[[iw]], st, ilua)
         #       dim(xw2)=c(nb_c, emuw, ntico)
         #xsim[inrow,]=st #xw2
         bop(xsim, c(1, 1L+nb_xi+nbc_x[iw], nb_row), "=", st)
         if (emu) {
            #xsim[1L+nb_xi+imwl,]=1.-arrApply(st, 2, "sum")
            bop(xsim, c(1, 1L+nb_xi+nbc_x[iw]+nb_row, nb_c), "=", 1.-arrApply(st, 2, "sum"))
         }
         if (cjac) {
#browser()
            # prepare jacobian ff, pf
            # rhs for all time points on this weight
            # parts before b_x%*%...
            #Rprof(file="fx2jr.Rprof", append=TRUE)
            #xpfw=double(nb_row*(nb_ff+nb_poolf)*ntico)
            resize(xpfw, c(nb_c*emuw, nb_ff+nb_poolf, ntico))
            bop(xpfw, 1, "=", 0.)
            #xpf1=double(nb_c*emuw*(nb_ff+nb_poolf))
            #resize(xpf1, c(nb_c, emuw*(nb_ff+nb_poolf)))
            #bop(xpf1, 1, "=", 0.) # xpf1 is always 0
            #sfpw=double(nb_row*ntico*nb_poolf)
            #resize(sfpw, c(nb_c, emuw, ntico, nb_poolf))
            #bop(sfpw, 1, "=", 0.)
            sj=fx2jr(fwrv, spa[[iw]], nb_f, xsim)
            #Rprof(NULL)
            # ff part
            if (nb_ff+nb_fgr > 0) {
               #tmp=sj$j_rhsw%stm%(-jx_f$df_dffp)
               #redim(tmp, c(nb_row, ntico, nb_ff+nb_fgr))
               ##xpfw[,seq_len(nb_ff+nb_fgr),]=-aperm(tmp, c(1L, 3L, 2L))
               #bop(xpfw, c(3, 0, nb_ff+nb_fgr), "=", aperm(tmp, c(1L, 3L, 2L)))
               jrhs_ff(sj$j_rhs, jx_f$df_dffp, xpfw)
               #if (sum(xpfw[,seq(nb_ff+nb_fgr),]-aperm(tmp, c(1L, 3L, 2L))) > 1.e-14)
               #   browser()
            }
            # poolf part
            if (nb_poolf > 0) {
#browser()
               i2x=nb_f$ipf2ircumo[[iexp]][[iw]]
               i2xf=nb_f$ipf2ircumo[[iexp]][[iw]]
               #bop(i2xf, c(2,2,1), "+=", nb_ff)
               i2xf[,3]=i2xf[,3]+nb_ff
               xw2=(cbind(x1[inrow], xsim[inrow, -ntico, drop=FALSE])-xsim[inrow, , drop=FALSE])%mrv%invdt
               redim(xw2, c(nb_c, emuw, ntico))
               #dim(xw2)=c(nb_c, emuw*ntico)
               #tmp=at%stm%xw2+c(st)
               #dim(tmp)=c(nb_c, emuw, ntico)
               #sfpw[i2x]=tmp[i2x[,-4L]]
               #bop(sfpw, i2x, "=", xw2[i2x[,-4L,drop=FALSE]])
               #xpfw[,nb_ff+seq_len(nb_poolf),]=(if (nb_fgr > 0) xpfw[,nb_ff+seq_len(nb_poolf,)] else 0.) + aperm(sfpw, c(1L, 2L, 4L, 3L))
               redim(xpfw, c(nb_c, emuw, nb_ff+nb_poolf, ntico))
               if (nb_fgr > 0) {
                  #bop(xpfw, c(2, nb_ff, nb_poolf), "+=", aperm(sfpw, c(1L, 2L, 4L, 3L)))
                  bop(xpfw, i2xf, "+=", xw2[i2x[,-3L,drop=FALSE]])
               } else {
                  #bop(xpfw, c(2, nb_ff, nb_poolf), "=", aperm(sfpw, c(1L, 2L, 4L, 3L)))
                  bop(xpfw, i2xf, "=", xw2[i2x[,-3L,drop=FALSE]])
                  redim(xpfw, c(nb_row, nb_ff+nb_poolf, ntico))
               }
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
            redim(vmw, c(nb_c, emuw*(nb_ff+nb_poolf)))
            redim(xpfw, c(nb_c, emuw*(nb_ff+nb_poolf), ntico))
            solve_ieu(invdt, NULL, vmw, ali_w[[iw]], xpfw, ilua)
            redim(xpfw, c(nb_c, emuw, nb_ff+nb_poolf, ntico))
            #xpfw[]=aperm(xpfw, c(1L, 2L, 4L, 3L))
#browser()
            #xpf[imw,,]=xpfw
            bop(xpf, c(1, nbc_x[iw], nb_row), "=", xpfw)
            if (emu) {
               # treat the last weight
               #xpf[imwl,,]=-arrApply(xpfw, 2, "sum")
               bop(xpf, c(1, nbc_x[iw]+nb_row, nb_c), "=", -arrApply(xpfw, 2, "sum"))
            }
         }
      } # iw loop
      #jx_f$xsim=xsim # for debugging only
      # get ti moments corresponding to measurements
      isel=itifu[tifull[[iexp]] %in% ti[[iexp]]][-1L]-1L
      ntise=length(isel)
      xsim=xsim[-seq_len(1L+nb_xi),, drop=FALSE]
      if (getx) {
         xsimf=xsim # full simulation (in time)
         xsim=xsim[,isel,drop=FALSE]
         # usm
         mx=measmat[[iexp]]%stm%(if (nrow(xsim) == nb_mcol) xsimf else xsimf[nm$rcumo_in_cumo,,drop=FALSE])+memaone[[iexp]]
         if (length(ipooled[[iexp]]) > 1L) {
            usmf=as.matrix(meas2sum[[iexp]]%stm%(pwe[[iexp]]*mx)) # full simulated measurements (in time)
         } else {
            usmf=mx
         }
         usm=usmf[,isel,drop=FALSE]
      } else {
         # usm
         mx=measmat[[iexp]]%stm%(if (nrow(xsim) == nb_mcol) xsim[,isel,drop=FALSE] else xsim[nm$rcumo_in_cumo,isel,drop=FALSE])+memaone[[iexp]]
         if (length(ipooled[[iexp]]) > 1L) {
            usm=as.matrix(meas2sum[[iexp]]%stm%(pwe[[iexp]]*mx))
         } else {
            usm=mx
         }
      }
      
      if (cjac) {
         #jx_f$xpf=xpf # for debugging only
         #xpf=aperm(xpf[,,isel,drop=FALSE], c(1L, 3L, 2L))
         #redim(xpf, c(nbc_x[nb_w+1L], ntise*(nb_ff+nb_poolf)))
         # scale part of jacobian
         if (nb_f$nb_sc_tot != 0L) {
            stop_mes("nb_sc != 0L is not implemented as meaningless for dynamic labeling. Use --noscale option.", file=fcerr)
         }
         #dux_dp=measmat[[iexp]]%stm%xpf
         dux_dp=mm_xpf(measmat[[iexp]], xpf, isel)
         #redim(xpf, c(nbc_x[nb_w+1L], ntise, (nb_ff+nb_poolf)))
         #dimnames(xpf)=list(nm$x, ti[[iexp]][-1], nm$par)
         nr=nrow(measmat[[iexp]])
         if (length(ipooled[[iexp]]) > 1L) {
            redim(dux_dp, c(dim(dux_dp)[1], ntise*(nb_ff+nb_poolf)))
            dux_dp=meas2sum[[iexp]]%stm%(pwe[[iexp]]*dux_dp) # resize
            nr=nrow(meas2sum[[iexp]])
         }
         redim(dux_dp, c(nr, ntise, nb_ff+nb_poolf))
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
         #jx_f$xpf[[iexp]]=xpf
      }
      dimnames(usm)=list(nm$meas[[iexp]], ti[[iexp]][-1L])
      # store usefull information in global list jx_f
      jx_f$param=param
      jx_f$usm[[iexp]]=usm
      if (getx) {
         jx_f$usmf[[iexp]]=usmf
         jx_f$xsim[[iexp]]=xsim
         jx_f$xsimf[[iexp]]=xsimf
      }
   }
   if (is.null(labargs$ali_w))
      labargs$ali_w=ali_w
   if (getx) {
      names(jx_f$usm)=names(jx_f$usmf)=names(jx_f$xsim)=names(jx_f$xsimf)=nm$nm_exp
      res=list(usm=jx_f$usm, usmf=jx_f$usmf, x=jx_f$xsim, xf=jx_f$xsimf, dux_dp=jx_f$dux_dp, lf=lf, df_dffp=jx_f$df_dffp)
   } else {
      names(jx_f$usm)=nm$nm_exp
      res=list(usm=jx_f$usm, dux_dp=jx_f$dux_dp, lf=lf, df_dffp=jx_f$df_dffp)
   }
   return(res)
}

param2fl_usm_rich=function(param, cjac, labargs) {
   # Richardson extrapolation to get order 2 in ODE solve
   if (labargs$time_order=="2") {
      getx=FALSE
      if (!is.null(labargs$getx))
         getx=labargs$getx
      labargs$labargs2$getx=getx
      # solve with step=h/2
      if (!is.null(labargs$cl) && is.null(.GlobalEnv$mc_iter)) {
         # do in parallel only out off MC iterations
         clusterExport(labargs$cl, c("param", "cjac"), envir=environment())
#print(labargs$cl[[1]])
#lila=parLapply(labargs$cl, c("labargs", "labargs"), function(nm) format(labargs))
#cat("lila=\n")
#print(lila)
#clusterEvalQ(labargs$cl, {cat("usm idth=", idth, "\n"); print(labargs)})
         if (getx)
            clusterEvalQ(labargs$cl, {labargs$labargs2$getx=labargs$getx=TRUE})
         else
            clusterEvalQ(labargs$cl, {labargs$labargs2$getx=labargs$getx=FALSE})
         res=clusterEvalQ(labargs$cl, if (idw < 3) param2fl_usm_eul2(param, cjac, if (idw == 1) labargs else labargs$labargs2))
         #res=parLapply(labargs$cl, seq(2), function(ith) {cat ("parlap ith=", ith, "\n"); print (labargs); cl_worker(funth=param2fl_usm_eul2, argth=list(param=param, cjac=cjac, labargs=if (ith == 1) labargs else labargs$labargs2))})
      } else {
         # do sequentially
         res=lapply(seq(2), function(i) param2fl_usm_eul2(param, cjac, if (i == 1) labargs else labargs$labargs2))
      }
      res1=res[[1]]
      if (!is.null(res1$err) && res1$err) {
         return(list(err=1, mes=res1$mes))
      }
      res2=res[[2]]
      if (!is.null(res2$err) && res2$err) {
         return(list(err=1, mes=res2$mes))
      }
      # Richardson interpolation
      jx_f=labargs$jx_f
      for (iexp in seq_len(labargs$nb_exp)) {
         jx_f$usm[[iexp]]=2*res2$usm[[iexp]]-res1$usm[[iexp]]
         if (getx) {
            jx_f$usmf[[iexp]]=2*res2$usmf[[iexp]][,seq.int(2, ncol(res2$usmf[[iexp]]), 2)]-res1$usmf[[iexp]]
            jx_f$xsim[[iexp]]=2*res2$x[[iexp]]-res1$x[[iexp]]
         }
         if (cjac) {
            jx_f$dux_dp[[iexp]]=2*res2$dux_dp[[iexp]]-res1$dux_dp[[iexp]]
            #jx_f$dux_dp=labargs$jx_f$dux_dp
         }
      }
      res1$usm[]=jx_f$usm
      res1$dux_dp[]=jx_f$dux_dp
      jx_f$df_dffp=res1$df_dffp
      if (getx) {
         res1$usmf[]=jx_f$usmf
         res1$x[]=jx_f$xsim
      }
   } else {
      res1=param2fl_usm_eul2(param, cjac, labargs)
   }
   return(res1)
}
# transform a funlab list 'li' into matrix #nm_inp x #tp having input labeling in time points 'tp'
# li entries: "metab name" -> "int iso mask" -> R-code depending on scalar 't'
# nm_inp cumo: "Glc:63"
# nm_inp emu: "Glc:63+0"
 
funlab=function(tp, nm_inp, li, env, emu, fname, fcerr, tol=.Machine$double.eps*2**7) {
   # lit is nested a list: met => str(isoint) => vector of legth #tp
   lit=lapply(structure(seq_along(li), names=names(li)),
   function(i) {
      m=li[[i]]
      met=names(li)[[i]]
      lapply(structure(names(m), names=names(m)),
      function(n) {
         vapply(tp,
         function(t) {
            env$t=t
#browser()
            v=try(eval(m[[n]], env), silent=TRUE) # time dependent isotopomers, i.e. functions applied on t
            if (inherits(v, "try-error")) {
               stop_mes("Error in R code '", format(m[[n]][[1L]]), "' for input label '", met, "#", n, "' from '", fname, "':\n", v, file=fcerr)
            }
            if (any(ibad <- is.na(suppressWarnings(as.double(v))))) {
               ibad=which(ibad)[1]
               stop_mes("Input label '", met, "#", n, "' from '", fname, "' produced a non numeric value at t=", tp[ibad], ": '", v[ibad], "'.", file=fcerr)
            }
            v[v < 0. && v >= -tol]=0.
            v[v > 1. && v <= 1+tol]=1.
            if (any(ibad <- v < 0.)) {
               ibad=which(ibad)[1]
               stop_mes("Input label '", met, "#", n, "' from '", fname, "' produced a negative value at t=", tp[ibad], ": '", v[ibad], "'.", file=fcerr)
            }
            if (any(ibad <- v > 1.)) {
               ibad=which(ibad)[1]
               stop_mes("Input label '", met, "#", n, "' from '", fname, "' produced a value > 1 at t=", tp[ibad], ": '", v[ibad], "'.", file=fcerr)
            }
            v
         }, double(1L)
         )
      }
      )
   }
   )
   # complete to 1 if necessary
   for (met in names(lit)) {
      m=lit[[met]]
      # sanity check, sum > 1
      su=Reduce("+", m)
      if (any(ibad <- su > 1+tol)) {
         ibad=which(ibad)[1L]
         stop_mes("Input labeled metabolite '", met, "' from '", fname, "' sums up to a value greater than 1 at t=", tp[ibad], ": '", su[ibad], "'.", file=fcerr)
      }
      if (!"0" %in% names(m)) {
         # "the rest is unlabeled"
         lit[[met]][["0"]]=1.-su
      } else if (length(m) == 2**clen[met]-1) {
         # "guess the lacking one"
         nm_all=as.character(seq(2**len[met])-1)
         lack=nm_all[which(!nm_all %in% names(m))]
         lit[[met]][[lack]]=1.-su
      } else if (any(ibad <- su < 1-tol)) {
         # sanity check, sum < 1
         ibad=which(ibad)[1L]
         stop_mes("Input labeled metabolite '", met, "' from '", fname, "' sums up to a value less than 1 at t=", tp[ibad], ": '", su[ibad], "'.", file=fcerr)
      }
   }
   sp=matrix(unlist(strsplit(nm_inp, ":", fixed=TRUE)), nrow=2)
   cres=matrix(0., nrow=length(nm_inp), ncol=length(tp))
   if (emu) {
      for (j in seq(ncol(sp))) {
         # find what isotopomers in lit contributes to this nm_inp
         met=sp[1L, j]
         iemu=as.integer(strsplit(sp[2L, j], "+", fixed=TRUE)[[1L]])
         iso=as.integer(names(lit[[met]]))
         i=sapply(iso, function(ii) sum(as.integer(intToBits(bitwAnd(ii, iemu[1])))) == iemu[2])
         res=double(length(tp))
         if (any(i))
            res=Reduce("+", lit[[met]][i])
         cres[j,]=res
      }
   } else {
      for (j in seq(ncol(sp))) {
         # find what isotopomers in lit contributes to this nm_inp
         met=sp[1L, j]
         icu=as.integer(sp[2L, j])
         iso=as.integer(names(lit[[met]]))
         i=sapply(iso, function(ii) bitwAnd(ii, icu) == icu)
         res=double(length(tp))
         if (any(i))
            res=Reduce("+", lit[[met]][i])
         cres[j,]=res
      }
   }
   return(cres)
}
