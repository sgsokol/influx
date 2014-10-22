#suppressPackageStartupMessages(require(deSolve))
suppressPackageStartupMessages(require(Matrix))

icumo_resid=function(param, cjac, labargs) {
   # claculates residual vector of labeling propagation corresponding to param
   #cat("icumo_resid: param=", param, ", cjac=", cjac, "\n")
#cjac=F # to remove in the final version
   # from labargs to local vars
   for (item in ls(labargs)) {
      assign(item, get(item, env=labargs))
   }
   stopifnot(length(tifull)!=length(tifull2))
   
   nb_w=length(spAb)
#browser()
   sqm=measurements$dev$labeled
   sqf=measurements$dev$flux
   sqp=measurements$dev$pool
   measvecti=measurements$vec$kin
   nb_ti=length(ti)
   nb_sc=nb_f$nb_sc
   nb_poolf=nb_f$nb_poolf
   
   jacobian=NULL
   # find usimvec with step h
   lres=lab_sim(param, cjac, labargs)
   if (!is.null(lres$err) && lres$err) {
      return(list(err=1, mes=lres$mes))
   }
   if (time_order == 2) {
      # find usimvec with step h/2
      labargs$tifull=tifull2
      labargs$jx_f=new.env()
      labargs$nb_f$ipf2ircumo=labargs$nb_f$ipf2ircumo2
      lres2=lab_sim(param, cjac, labargs)
      labargs$tifull=tifull
      jx_f2=labargs$jx_f
      labargs$jx_f=jx_f
      labargs$nb_f$ipf2ircumo=nb_f$ipf2ircumo
      if (!is.null(lres2$err) && lres2$err) {
         return(list(err=1, mes=lres2$mes))
      }
      # Richardson interpolation
      jx_f$usm=2*jx_f2$usm-jx_f$usm
      if (cjac) {
         jx_f$dux_dp=2*jx_f2$dux_dp-jx_f$dux_dp
         #jx_f$dux_dp=labargs$jx_f$dux_dp
      }
   }
   
   # scale simulated measurements scale*(usm)
   simlab=jx_f$usm
   if (nb_sc > 0) {
      vsc=c(1.,param)[ir2isc]
      simlab=simlab*vsc
   }
   jx_f$simlab=simlab
   # diff between simulated and measured
   #inna=which(!is.na(measvecti)) # for removing NA measurements
   if (is.null(measvecti)) {
      jx_f$res=NULL
   } else {
      pool[nm$poolf]=param[nm$poolf]
      jx_f$ureslab=simlab-measvecti
      jx_f$simfmn=jx_f$lf$fallnx[nm$fmn]
      jx_f$uresflu=jx_f$simfmn-measurements$vec$flux
      jx_f$simpool=as.matrix(measurements$mat$pool%*%pool)
      jx_f$urespool=jx_f$simpool - measurements$vec$pool
      jx_f$reslab=jx_f$ureslab/sqm
      jx_f$resflu=jx_f$uresflu/sqf
      jx_f$respool=jx_f$urespool/sqp
      jx_f$res=c(jx_f$reslab, jx_f$resflu, jx_f$respool)
      names(jx_f$res)=nm$resid
      jx_f$ures=c(jx_f$ureslab, jx_f$uresflu, jx_f$urespool)
      names(jx_f$ures)=nm$resid
   }
#browser()
   if (cjac) {
      # calculate jacobian
      # first, its labeling part (unscaled and unreduced)
      #dusm=dusm_dffpf(param, labargs)
      #dux_dp=cBind(dusm[,iseq(nb_ff),drop=F], Matrix(0., nrow=nrow(dusm), ncol=nb_sc), if (nb_ff > 0) dusm[,-iseq(nb_ff),drop=F] else dusm)
      # next, flux part and measured pools part
      nb_lab=nb_f$nb_meas*(nb_ti-1L)
      if (is.null(jx_f$jacobian)) {
         jacobian=matrix(0., nrow=nb_lab+nb_f$nb_fmn+nb_f$nb_poolm, ncol=nb_f$nb_ff+nb_sc+nb_poolf)
         dimnames(jacobian)=list(nm$resid, nm$par)
         # constant part of jacobian
         jacobian[nb_lab+iseq(nb_f$nb_fmn),]=dufm_dp
         jacobian[nb_lab+nb_f$nb_fmn+iseq(nb_f$nb_poolm),]=dupm_dp
      } else {
         jacobian=jx_f$udr_dp
      }
      jacobian[iseq(nb_lab),]=jx_f$dux_dp
      #jacobian[nb_lab+iseq(nb_f$nb_fmn),]=dufm_dp
      #jacobian[nb_lab+nb_f$nb_fmn+iseq(nb_f$nb_poolm),]=dupm_dp
      #as.matrix(rBind(jx_f$dux_dp, dufm_dp, dupm_dp))
      jx_f$udr_dp=jacobian
      jx_f$df_dffp=df_dffp(param, lres$lf$flnx, nb_f, nm_list)
      # reduce it
      jacobian[]=with(measurements$dev, jacobian/c(kin, flux, pool))
      if (nb_sc > 0) {
         # scale it
         i=iseq(length(measvecti))
         jacobian[i,]=vsc*jacobian[i,]
         # add scale part of jacobian (the sparsity pattern doesn't change)
         jx_f$dr_dsc[nb_f$is2mti]=jx_f$res[is2mti[,1]]
         # combine three parts: ff, scale, poolf
         jacobian[]=as.matrix(cBind(jacobian[,iseq(nb_ff)], jx_f$dr_dsc, if (nb_ff > 0) jacobian[,-iseq(nb_ff)] else jacobian))
      }
      
      jx_f$jacobian=jacobian
      jx_f$dr_dff=jacobian[,iseq(nb_ff),drop=F]
   }
   
   return(list(res=jx_f$res, jacobian=if (cjac) jx_f$jacobian else NULL))
}

icumo_cost=function(param, labargs) {
   resl=icumo_resid(param, cjac=FALSE, labargs)
   if (!is.null(resl$err) && resl$err) {
      return(NULL)
   }
   res=resl$res
   #res=labargs$jx_f$res
   fn=sum(res*res)
   return(fn)
}

fwrv2sp=function(fwrv, spAbr, incu, incup=NULL, gets=TRUE, emu=FALSE) {
   # calculate s and ds/dt in A*x+s (where x is cumomer vector
   # and xp its first derivative in time)
   # from fields of the spAbr
   # according to conventions explained in comments to python function
   # netan2Abcumo_spr() generating spAbr
   # return a list with s and sp
   # 2012-03-07 sokol
   # 2014-04-11 sokol, added emu option
   
   # construct the sources s and its derivatives (if xp not null)
   # for this weight
   nb_c=spAbr$nb_c # cumomer or fragment number (when emu==T)
   w=spAbr$w # cumomer weight
   incu=as.matrix(incu)
   nco=ncol(incu)
   
   if (gets) {
      if (nb_c == 0) {
         s=Matrix(0., nb_c, 1)
         return(list(s=s), sp=if (is.null(incup)) NULL else s)
      }
      if (emu) {
         ind_b=spAbr[["ind_b_emu"]]
         s=c(fwrv[ind_b[,"indf"]]*incu[ind_b[,"indx1"],]*incu[ind_b[,"indx2"],])
         i=rep(ind_b[,"irow"]+(ind_b[,"iwe"]-1L)*nb_c, nco)
         j=rep(iseq(nco), each=nrow(ind_b))
         s=sparseMatrix(i=i, j=j, x=s, dims=c(nb_c*w, nco))
      } else {
         ind_b=spAbr[["ind_b"]]
         s=c(fwrv[ind_b[,"indf"]]*incu[ind_b[,"indx1"],]*incu[ind_b[,"indx2"],])
         s=sparseMatrix(i=rep(ind_b[,"irow"], nco), j=rep(iseq(nco), each=nrow(ind_b)), x=s, dims=c(nb_c, nco))
      }
   } else {
      s=NULL
   }
   # sp
   if (!is.null(incup)) {
      if (emu) {
         ind_b=spAbr[["ind_b_emu"]]
         sp=fwrv[ind_b[,"indf"]]*(incup[ind_b[,"indx1"]]*incu[ind_b[,"indx2"]] + incu[ind_b[,"indx1"]]*incup[ind_b[,"indx2"]])
      } else {
         ind_b=spAbr[["ind_b"]]
         sp=fwrv[ind_b[,"indf"]]*(xp[ind_b[,"indx1"]]*x[ind_b[,"indx2"]] + x[ind_b[,"indx1"]]*xp[ind_b[,"indx2"]])
         sp=sparseMatrix(i=ind_b[,"irow"], j=rep.int(1, nrow(ind_b)), x=sp, dims=c(nb_c, 1))
      }
   } else {
      sp=NULL
   }
   return(list(s=s, sp=sp))
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
   
   
   # from labargs to local vars
   for (item in ls(labargs)) {
      assign(item, get(item, env=labargs))
   }
   nb_w=length(spAb)
   nb_ti=length(ti)
   nb_tifu=length(tifull)
   if (nb_ti < 2) {
      return(list(err=1, mes="Number of time points is less than 2"))
   }
   if (!all(ti %in% tifull)) {
      return(list(err=1, mes="Not all time moments in ti are present in tifull vector"))
   }
   
   dt=diff(tifull)
   stopifnot(all(dt > 0.))
   ntico=nb_tifu-1L # number of time columns
   idt=iseq(ntico)
   itifu=iseq(nb_tifu)
   
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
   nb_xi=length(xi)
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   nb_fgr=nb_f$nb_fgr
   nb_meas=nb_f$nb_meas
   nb_sc=nb_f$nb_sc
   # fullfill pool with free pools
   if (nb_poolf > 0) {
      pool[nm$poolf]=param[nm$poolf]
   }
   
   nb_mcol=ncol(measmat)
   # prepare ponderation with actual metab pools
   pwe[ipwe]=pool[ip2ipwe]
   spwe=tapply(pwe, pool_factor, sum)
   spwe=1./as.numeric(spwe[nm$measmat])
   pwe=pwe*spwe
   
   # prepare inverse of pool vectors
   invpool=1./pool
   # invm has the same length as the full label vector
   invm=invpool[nb_f$ip2ircumo]
#browser()
   # prepare vectors at t1=0 with zero labeling
   # incu, xi is supposed to be in [0; 1]
   x1=c(1., xi, rep(0., nbc_x[nb_w+1])) # later set m+0 to 1 in x1 later
   names(x1)=c("one", nm$inp, nm$x)

   xsim=matrix(x1, nrow=length(x1), ncol=ntico)
   dimnames(xsim)=list(names(x1), tifull[-1L])
   if (cjac) {
      #cat("param2fl_usm_eul2: recalculate jacobian\n")
      mdf_dffp=df_dffp(param, lf$flnx, nb_f, nm_list)
      jx_f$df_dffp=mdf_dffp
      xpf=double(nbc_x[nb_w+1L]*ntico*(nb_ff+nb_poolf))
      dim(xpf)=c(nbc_x[nb_w+1L], ntico, nb_ff+nb_poolf)
      if (length(ijpwef)) {
         dpwe=-pwe*spwe
         dpwe[-ipwe]=0.
         dpwe=(dpwe+dp_ones*spwe)
      }
   }
   # prepare data for iadt array
   dtr=as.character(round(dt, 10L))
   dtru=unique(dtr)
   ilua=pmatch(dtr, dtru, dup=T)
#browser()
   for (iw in iseq(nb_w)) {
      emuw=ifelse(emu, iw, 1L)
      nb_c=spAb[[iw]]$nb_c
      ixw=nbc_x[iw]+iseq(nb_x[iw])
      inxw=(1L+nb_xi)+ixw
      nb_row=nb_c*emuw
      inrow=(1L+nb_xi+nbc_x[iw])+iseq(nb_row)
      imw=nbc_x[iw]+iseq(nb_row) # mass index in x (all but last)
      if (cjac) {
         xpfw=double(nb_row*(nb_ff+nb_poolf)*ntico)
         dim(xpfw)=c(nb_row, nb_ff+nb_poolf, ntico)
         #xpf1=matrix(0., nb_c, emuw*(nb_ff+nb_poolf))
         sfpw=double(nb_row*ntico*nb_poolf)
         dim(sfpw)=c(nb_c, emuw, ntico, nb_poolf)
      }
      invmw=invm[nbc_cumos[iw]+iseq(nb_c)]
      Aw=fwrv2Abr(fwrv, spAb[[iw]], x1, nm$x[nbc_x[iw]+iseq(nb_x[iw])], getb=F,  emu=emu)$A*invmw
      # prepare (I-a*dt)^-1
      eye=diag(nb_c)
      ali=mclapply(dt[pmatch(dtru, dtr)], function(dti) lu(eye-as.matrix(Aw*dti)))
      lua=vapply(ali, function(item) item@x, eye)
      dim(lua)=c(nb_c, nb_c, length(dtru))
      perm=vapply(ali, function(item) item@perm, integer(nb_c))
      if (emu) {
         # for the first time point, set m+0 to 1 in x1
         x1[(1L+nb_xi+nbc_x[iw])+iseq(nb_c)]=1.
         xsim[inxw,1L]=x1[inxw]
         imwl=nbc_x[iw]+nb_row+iseq(nb_c) # the last mass index in x
      }
      # source terms
      st=as.matrix(fwrv2sp(fwrv, spAb[[iw]], xsim, emu=emu)$s*invmw)
      stt=st%mrv%dt
      # calculate labeling for all time points
      #xw1=matrix(x1[inrow], nb_c, emuw)
      stt[,1L]=stt[,1L]+x1[inrow]
      #xsim[inxw,]=vapply(idt, function(idtr) {
      #   xw2=lusolve(iadt[[dtr]], (xw1+st[,idtr]))
      #   xw1[] <<- xw2
      #   if (emu) {
      #      return(c(xw2, 1.-rowSums(xw2)))
      #   } else {
      #      return(xw2)
      #   }
      #}, x1[inxw])
      dim(stt)=c(nb_c, emuw, ntico)
      solve_lut(lua, perm, stt, ilua, dirx)
      dim(stt)=c(nb_c, emuw, ntico)
      xsim[inrow,]=stt
      if (emu) {
         stt=aperm(stt, c(2L, 1L, 3))
         xsim[(1L+nb_xi)+imwl,]=1.-colSums(stt)
      }
      if (cjac) {
         # prepare jacobian ff, pf
         # rhs for all time points on this weight
         # parts before b_x%*%...
         sj=fx2jr(fwrv, spAb[[iw]], nb_f, xsim)
         # ff part
         if (nb_ff+nb_fgr > 0) {
            tmp=as.numeric(sj$j_rhsw%*%mdf_dffp)
            dim(tmp)=c(nb_row, ntico, nb_ff+nb_fgr)
            xpfw[,iseq(nb_ff+nb_fgr),]=aperm(tmp, c(1L, 3L, 2L))
         }
         # poolf part
         if (nb_poolf > 0) {
            i2x=nb_f$ipf2ircumo[[iw]]
            xw2=xsim[inrow, , drop=F]
            dim(xw2)=c(nb_c, emuw*ntico)
            tmp=as.matrix(Aw%*%xw2)+c(st)
            dim(tmp)=c(nb_c, emuw, ntico)
            sfpw[i2x]=tmp[i2x[,-4L]]
            xpfw[,nb_ff+iseq(nb_poolf),]=if (nb_fgr > 0) xpfw[,nb_ff+iseq(nb_poolf,)] else 0. + aperm(sfpw, c(1L, 2L, 4L, 3L))
         }
         # add lighter xpf (b_x%*%...)
         if (iw > 1L) {
#browser()
            ir=iseq(nb_row)
            ic=iseq(nbc_x[iw])
            #pti=proc.time()
            #for (idtr in idt) {
            #   xpfw[,,idtr]=xpfw[,,idtr]+ as.double(sj$b_x[ir+(idtr-1L)*nb_row,]%*%xpf[ic,idtr,])
            mult_bxx(xpfw, sj$b_x, xpf, ntico, dirx)
            #}
            #cat("iw=", iw, "; len(x)=", length(sj$b_x[ir,]@x), "\n", sep="")
            #print(proc.time()-pti)
         }
         # solve the system at each time point
#browser()
         xpfw[]=(rep(-invmw, emuw)%o%rep(1., nb_ff+nb_poolf)%o%dt)*xpfw
         #xpfw=vapply(idt, function(idtr) {
         #   #dtr=as.character(round(dt[idtr], 10L))
         #   xpf1[] <<- lusolve(ali[[pmatch(dtr[idtr], dtru)]], c(xpfw[,,idtr])+xpf1) # +xpf1 makes that it is a matrix of a suitable size            xpf1[]=xpfw[,idtr,]
         #}, xpf1)
         dim(xpfw)=c(nb_c, emuw*(nb_ff+nb_poolf), ntico)
         solve_lut(lua, perm, xpfw, ilua, dirx)
         dim(xpfw)=c(nb_c, emuw, nb_ff+nb_poolf, ntico)
         xpfw=aperm(xpfw, c(1L, 2L, 4L, 3L))
         xpf[imw,,]=xpfw
         if (emu) {
            # treat the last weight
            #xpf[imwl,,]=-apply(xpfw, c(1L, 3L, 4L), sum)
            xpfw=aperm(xpfw, c(2L, 1L, 3L, 4L))
            xpf[imwl,,]=-colSums(xpfw)
         }
      }
   } # iw loop
   #jx_f$xsim=xsim # for debugging only
   # get ti moments corresponding to measurements
   isel=itifu[tifull %in% ti][-1L]-1L
   ntise=length(isel)
   xsim=xsim[-iseq(1L+nb_xi), isel, drop=F]
   # usm
   mx=measmat%*%(if (nrow(xsim) == nb_mcol) xsim else xsim[nm$rcumo_in_cumo,,drop=F])+memaone
   if (length(ipooled) > 1L) {
      usm=as.matrix(meas2sum%*%(pwe*mx))
   } else {
      usm=mx
   }
   if (cjac) {
      #jx_f$xpf=xpf # for debugging only
      xpf=xpf[,isel,,drop=F]
      dim(xpf)=c(nbc_x[nb_w+1L], ntise*(nb_ff+nb_poolf))
      # scale part of jacobian
      if (nb_f$nb_sc != 0L) {
         stop_mes("nb_sc != 0L is not implemented as meaningless for dynamic labeling", fcerr)
      }
      dux_dp=measmat%*%xpf
      dim(xpf)=c(nbc_x[nb_w+1L], ntise, (nb_ff+nb_poolf))
      dimnames(xpf)=list(nm$x, ti[-1], nm$par)
      if (length(ipooled) > 1L) {
         dux_dp=meas2sum%*%(pwe*dux_dp) # resize
      }
      dim(dux_dp)=c(nb_meas, ntise, nb_ff+nb_poolf)
      if (length(ijpwef) > 0L) {
         # derivation of pool ponderation factor
         dpw_dpf=double(nrow(measmat)*ntise*nb_poolf)
         dim(dpw_dpf)=c(nrow(measmat), ntise, nb_poolf)
         tmp=c(mx)*dpwe
         dim(tmp)=dim(dpw_dpf)
         dpw_dpf[ijpwef]=tmp[ijpwef]
         dim(dpw_dpf)= c(nrow(measmat), ntise*nb_poolf)
         dux_dp[,,nb_ff+iseq(nb_poolf)]=dux_dp[,,nb_ff+iseq(nb_poolf)]+c(meas2sum%*%dpw_dpf)
      }
      dimnames(dux_dp)=list(nm$meas, ti[-1L], nm$par)
      jx_f$dux_dp=dux_dp
      jx_f$xpf=xpf
   }
   dimnames(usm)=list(nm$meas, ti[-1L])
   # store usefull information in global list jx_f
   jx_f$param=param
   jx_f$usm=usm
   jx_f$xsim=xsim
   
   return(list(usm=usm, x=xsim, lf=lf))
}
mult_bxx=function(a, bx, c, ntico, dirx) {
   # R wrapper for a fortran call mult_bxt()
   if (!is.loaded("mult_bxt")) {
      dyn.load(sprintf("%s/mult_bxx%s", dirx, .Platform$dynlib.ext))
   }
   stopifnot(class(bx)=="dgCMatrix")
   .Fortran("mult_bxt", bx@x, bx@i, bx@p, as.integer(nrow(bx)/ntico), as.integer(ntico), ncol(bx), c, dim(c)[1L], dim(c)[3L], a,
      NAOK=T, DUP=F
   )
   return(NULL) # the matrix a is modified in place.
}
solve_lut=function(lua, pivot, b, ilua, dirx) {
   # call lapack dgters() for solving a series of linea systems a_i%*%(b_{i-1}+b_i)=b_i
   # The result is stored inplace in b.
   # The sizes:
   #  - LU matrices of a : (nr_a, nr_a, nlua)
   #  - pivot : (nr_a, nlua)
   #  - b : (nr_a, nc_c, ntico)
   #  - ilua : (ntico). It is a index vector. ilua[i] indicate which lua corresponds to the i-th time point
   if (!is.loaded("dgetrs")) {
      lapack.path <- file.path(R.home(), ifelse(.Platform$OS.type == "windows",
         file.path("bin", "Rlapack"), file.path("lib", "libRlapack")))
      dyn.load(paste(lapack.path,.Platform$dynlib.ext, sep=""))
   }
   if (!is.loaded("mult_bxt")) {
      dyn.load(sprintf("%s/mult_bxx%s", dirx, .Platform$dynlib.ext))
   }
   dlu=dim(lua)
   db=dim(b)
   .Fortran("solve_lut", lua, dlu[1L], dlu[3L], pivot, b, db[2L], db[3L], ilua,
      NAOK=F, DUP=F)
   return(NULL)
}
