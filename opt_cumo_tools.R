#TIMEIT=0; # 1 to enable time printing at some stages
if (length(find("TIMEIT")) && TIMEIT && length(find("fclog"))) {
   cat("load    : ", format(Sys.time()), " cpu=", proc.time()[1], "\n", sep="", file=fclog)
}
# get compiled code
#browser()
so=.Platform$dynlib.ext
fso=paste("mult_bxxc", so, sep="")
fcpp="mult_bxxc.cpp"
if (!file.exists(file.path(dirx, fso)) ||
   file.mtime(file.path(dirx, fso)) < file.mtime(file.path(dirx, fcpp))) {
   # freshly compile the code (==first use or .so is outdated)
   Sys.setenv(PKG_LIBS=file.path(system.file(package="rmumps"), "libs", paste("rmumps", .Platform$dynlib.ext, sep="")))
   Sys.setenv(PKG_CXXFLAGS="-std=c++11")
   tes=capture.output(sourceCpp(file.path(dirx, "mult_bxxc.cpp"), verbose=TRUE))
   ftmp=sub(".*'(.*)'.*", "\\1", grep("dyn.load", tes, value=TRUE))
   file.copy(ftmp, file.path(dirx, fso), copy.date=TRUE)
}
# define R functions from mult_bxxc.so
multdll=dyn.load(file.path(dirx, fso))
mult_bxxc <- Rcpp:::sourceCppFunction(function(a, b, c) {}, TRUE, multdll, 'sourceCpp_0_mult_bxxc')
solve_ieu <- Rcpp:::sourceCppFunction(function(invdt, x0, M, ali, s, ilua) {}, TRUE, multdll, 'sourceCpp_0_solve_ieu')
match_ij <- Rcpp:::sourceCppFunction(function(ix, jx, ti, tj) {}, TRUE, multdll, 'sourceCpp_0_match_ij')
rm(multdll)

dfcg2fallnx=function(nb_f, flnx, param, fc, fg) {
   # produce complete flux (net,xch)*(dep,free,constr,growth) vector
   # from dep,free,constr,growth
   f=c(flnx[seq_len(nb_f$nb_fln)], param[seq_len(nb_f$nb_ffn)], fc[seq_len(nb_f$nb_fcn)], fg,
      flnx[nb_f$nb_fln+seq_len(nb_f$nb_flx)], param[nb_f$nb_ffn+seq_len(nb_f$nb_ffx)], fc[nb_f$nb_fcn+seq_len(nb_f$nb_fcx)], numeric(nb_f$nb_fgr))
   return(f)
}

cumo_resid=function(param, cjac=TRUE, labargs) {
   nm_local=c("jx_f", "pool", "measurements", "nm", "ir2isc", "jacobian", "nb_f")
   for (item in nm_local) {
      assign(item, labargs[[item]])
   }
   # find x for all weights
   lres=lab_sim(param, cjac, labargs)
   if (!is.null(lres$err) && lres$err) {
      return(list(err=1, mes=lres$mes))
   }

   # find simulated scaled measure vector scale*(measmat*x)
   jx_f$simlab=jx_f$usimcumom*c(1.,param)[ir2isc]

   # diff between simulated and measured
   pool[nm$poolf]=param[nm$poolf]
   jx_f$simfmn=jx_f$lf$fallnx[nm$fmn]
   jx_f$simpool=as.numeric(measurements$mat$pool%*%pool)

   jx_f$ureslab=(jx_f$simlab-measurements$vec$labeled)
   jx_f$reslab=jx_f$ureslab/measurements$dev$labeled
   jx_f$uresflu=jx_f$simfmn-measurements$vec$flux
   jx_f$resflu=jx_f$uresflu/measurements$dev$flux
   jx_f$urespool=jx_f$simpool-measurements$vec$pool
   jx_f$respool=jx_f$urespool/measurements$dev$pool
   jx_f$res=c(jx_f$reslab, jx_f$resflu, jx_f$respool)
   jx_f$res[measurements$outlier]=NA
   jx_f$ures=c(jx_f$ureslab, jx_f$uresflu, jx_f$urespool)
   if (cjac) {
#browser()
      # jacobian
      cumo_jacob(param, labargs)
      if (is.null(labargs$jacobian)) {
         jacobian=matrix(0., length(nm$resid), length(param))
         dimnames(jacobian)=list(nm$resid, nm$par)
         labargs$jacobian=jacobian
      }
#require(numDeriv)
#r=function(p) cumo_resid(p, F, labargs)$res
#jacobian=jacobian(r, param)
      jacobian[]=as.numeric(with(measurements$dev, jx_f$udr_dp/c(labeled, flux, pool)))
      jx_f$jacobian=jacobian
      jx_f$dr_dff=jacobian[,seq_len(nb_f$nb_ff),drop=F]
      labargs$jacobian=jacobian
      return(list(res=jx_f$res, jacobian=jacobian))
   } else {
      return(list(res=jx_f$res))
   }
}

cumo_cost=function(param, labargs, resl=lab_resid(param, cjac=FALSE, labargs)) {
   if (!is.null(resl$err) && resl$err) {
      return(NULL)
   }
   res=resl$res
   if (is.null(res))
      return(NA)
   iva=!is.na(res)
   vres=res[iva]
   fn=crossprod(vres)[1L]
   return(fn)
}

param2fl=function(param, labargs) {
   # claculate all fluxes from free fluxes
   
   # local variabl assignments form labargs
   nm_local=c("nb_f", "nm", "invAfl", "p2bfl", "bp", "g2bfl", "fc")
   for (item in nm_local) {
      assign(item, labargs[[item]])
   }
   
   fg=numeric(nb_f$nb_fgr)
   names(fg)=nm$fgr
   if (nb_f$nb_fgr > 0) {
      fg[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   }
   flnx=as.numeric(invAfl%*%(p2bfl%*%param[seq_len(nb_f$nb_ff)]+c(bp)+g2bfl%*%fg))
   names(flnx)=nm$flnx
   fallnx=c(dfcg2fallnx(nb_f, flnx, param, fc, fg))
   names(fallnx)=nm$fallnx
   fwrv=c(fallnx2fwrv(fallnx, nb_f))
   names(fwrv)=nm$fwrv
   lf=list(fallnx=fallnx, fwrv=fwrv, flnx=flnx)
   labargs$jx_f$lf=lf
   return(lf)
}

param2fl_x=function(param, cjac=TRUE, labargs) {
   # translate free params (fluxes+scales) to fluxes and cumomers
   # or emus
   # local variabl assignments form labargs
   nm_local=ls(labargs)
   for (item in nm_local) {
      assign(item, labargs[[item]])
   }
   
   nb_xi=length(xi)
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   nb_fgr=nb_f$nb_fgr
   nb_meas=nb_f$nb_meas
   nb_sc=nb_f$nb_sc
   fg=nb_f$mu*param[nm$poolf] # the same alphabetical order
   names(fg)=nm$fgr
   pool[nm$poolf]=param[nm$poolf] # inject variable pools to pool vector

   # calculate all fluxes from free fluxes
   lf=param2fl(param, labargs)
   # prepare measurement pooling operations
   pwe[ipwe]=pool[ip2ipwe]
   spwe=tapply(pwe, pool_factor, sum)
   spwe=1./as.numeric(spwe[nm$measmat])
   pwe=pwe*spwe
   # construct the system A*x=b from fluxes
   # and find x for every weight
   # if fj_rhs is not NULL, calculate jacobian x_f
   if (is.null(labargs$incu)) {
      labargs$incu=incu=c(1, xi, double(nbc_x[nb_w+1L]))
   }
   if (cjac) {
      # derivation of fwrv fluxes by free parameters: free fluxes+concentrations
      mdf_dffp=df_dffp(param, jx_f$lf$flnx, nb_f, nm)
      jx_f$df_dffp=mdf_dffp
      if (is.null(labargs$x_f)) {
         labargs$x_f=x_f=matrix(0., nrow=sum(nb_x), ncol=nb_ff+nb_fgr)
      }
      
      if (length(ijpwef)) {
         dpwe=-pwe*spwe
         # dpwe is shortened to non zero entries in dpw_dpf
         dpwe=ifelse(ipf_in_ppw[ijpwef[,1L]]==ijpwef[,2L], (spwe+dpwe)[ijpwef[,1L]], dpwe[ijpwef[,1L]])
      }
   } else {
      x_f=NULL
   }
   
   # simulate labeling weight by weight
   ba_x=0
   for (iw in seq_len(nb_w)) {
      nb_c=spAb[[iw]]$nb_c
      emuw=ifelse(emu, iw, 1L)
      if (nb_c == 0) {
         next
      }
      ixw=(nbc_x[iw]+1L):nbc_x[iw+1]
      incuw=(1L+nb_xi)+ixw
      if (emu) {
         lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], incu, nm$emu[ixw], emu=emu)
      } else {
         lAb=fwrv2Abr(lf$fwrv, spAb[[iw]], incu, nm$rcumo[ixw], emu=emu)
      }
      #if (use_mumps) {
#browser()
      xw=try(solve(lAb$A, lAb$b), silent=TRUE)
      if (inherits(xw, "try-error")) {
         # find 0 rows if any
         l=spAb[[iw]]
         ag=aggregate(abs(lf$fwrv[l$ind_a[,"indf"]]), list(l$ind_a[,"ir0"]), sum)
         izc=ag$Group.1[ag$x <= 1.e-10]
         izf=names(which(abs(lf$fwrv)<1.e-7))
         mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
            "Zero rows in cumomer matrix A at weight ", iw, ":\n",
            paste(nm$rcumo[ixw][izc+1], collapse="\n"), "\n",
            "Zero fluxes are:\n",
            paste(izf, collapse="\n"), "\n",
            sep="")
         #izc=apply(lAb$A, 1L, function(v)sum(abs(v))<=1.e-10)
         #izf=names(which(abs(lf$fwrv)<1.e-7))
         #if (sum(izc) && length(izf)) {
         #   mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
         #      "Zero rows in cumomer matrix A at weight ", iw, ":\n",
         #      paste(rownames(lAb$A)[izc], collapse="\n"), "\n",
         #      "Zero fluxes are:\n",
         #      paste(izf, collapse="\n"), "\n",
         #      sep="")
         #} else {
         #   mes="Cumomer matrix is singular.\n"
         #}
         return(list(x=NULL, fA=lAb$A, err=1L, mes=mes))
      }
      #} else {
      #   wa=options(warn=2) # to cath singular matrix as error and not just a warning
      #   lua=try(if (use_magma) magma::lu(magma(as.matrix(lAb$A))) else lu(as.matrix(lAb$A), errSing=T), silent=T)
      #   options(wa) # restore warning situation
      #   if (inherits(lua, "try-error")) {
      #      # find 0 rows if any
      #      izc=apply(lAb$A, 1L, function(v)sum(abs(v))<=1.e-10)
      #      izf=names(which(abs(lf$fwrv)<1.e-7))
      #      if (sum(izc) && length(izf)) {
      #         mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
      #            "Zero rows in cumomer matrix A at weight ", iw, ":\n",
      #            paste(rownames(lAb$A)[izc], collapse="\n"), "\n",
      #            "Zero fluxes are:\n",
      #            paste(izf, collapse="\n"), "\n",
      #            sep="")
      #      } else {
      #         mes="Cumomer matrix is singular.\n"
      #      }
#browser()
      #      return(list(x=NULL, fA=lAb$A, err=1L, mes=mes))
      #   }
      #   b=as.matrix(lAb$b); # may have several columns if emu is TRUE
      #   #solve the system A*x=b
      #   #lsolv=trisparse_solv(lAb$A, lAb$b, iw, lf$fwrv, method="sparse")
      #   xw=lusolve(lua, b)
      #}
      if (emu) {
         xw=c(xw, 1.-rowSums(xw))
      }
      incu[incuw]=xw
      if (cjac) {
         # calculate jacobian x_f
         # first, calculate right hand side for jacobian calculation
         # j_rhsw, b_x from sparse matrices
         # bind cumomer vector
         j_b_x=fx2jr(jx_f$lf$fwrv, spAb[[iw]], nb_f, incu)
         j_rhsw=j_b_x$j_rhsw%stm%mdf_dffp
         b_x=j_b_x$b_x
         if (iw > 1) {
            if (ba_x > 0) {
               j_rhsw=j_rhsw+b_x%stm%x_f[1L:ba_x,,drop=F]
            }
         }
         #j_rhsw=as.double(j_rhsw)
         dim(j_rhsw)=c(nb_c, emuw*(nb_ff+nb_fgr))
         #if (use_mumps) {
         tmp=try(solve(lAb$A, j_rhsw))
         if (inherits(tmp, "try-error")) {
            #browser()
            mes="Some obscure problem with label matrix.\n"
            return(list(x=NULL, fA=lAb$A, err=1L, mes=mes))
         } else {
            j_rhsw=tmp
         }
         #} else {
         #   j_rhsw=lusolve(lua, j_rhsw)
         #}
         if (emu) {
            dim(j_rhsw)=c(nb_c, iw, nb_ff+nb_fgr)
            x_f[ba_x+seq_len(iw*nb_c),]=j_rhsw
            # m+N component
            #x_f[ba_x+iw*nb_c+seq_len(nb_c),]= -apply(j_rhsw, c(1L,3L), sum)
            x_f[ba_x+iw*nb_c+seq_len(nb_c),]=-arrApply(j_rhsw, 2, "sum")
         } else {
            x_f[ba_x+seq_len(nb_c),]=j_rhsw
         }
      }
      ba_x=ba_x+nb_x[iw]
   }
   names(incu)=c("one", nm$inp, nm$x)
   x=tail(incu, -nb_xi-1L)
   jx_f$x=x
   
   # calculate unreduced and unscaled measurements
   if (length(x) == ncol(measmat)) {
      mx=measmat%*%x+memaone
   } else {
      mx=measmat%*%x[nm$rcumo_in_cumo]+memaone
   }
   if (length(ipooled) > 1L) {
      mv=meas2sum%*%(pwe*mx)
   } else {
      mv=mx
   }
   jx_f$usimcumom=as.numeric(mv)
   names(jx_f$usimcumom)=nm$meas
   if (cjac) {
      # unreduced residuals derivated by scale params
      if (is.null(labargs$dur_dsc)) {
         labargs$dur_dsc=dur_dsc=matrix(0., nrow=nb_meas, ncol=nb_sc)
      }
      if (is.null(labargs$dux_dp)) {
         labargs$dux_dp=dux_dp=matrix(0., nb_meas, nb_ff+nb_sc+nb_poolf)
      }
      # measurement vector before pool ponderation
      # scale part of jacobian
      if (nb_sc > 0) {
         is2m=nb_f$is2m
         dur_dsc[is2m]=mv[is2m[,1]]
      }
#browser()
      # free flux part of jacobian (and partially free pools if present in x_f)
      if (nb_ff+nb_fgr > 0) {
         mffg=measmat%*%x_f
         if (length(ipooled) > 1L) {
            mffg=meas2sum%*%(pwe*mffg)
         }
      } else {
         mffg=matrix(0., nrow=nb_meas, ncol=0L)
      }
      # free pool part of jacobian
      mff=mffg
      mpf=matrix(0., nrow=nb_meas, ncol=nb_f$nb_poolf)
      if (nb_f$nb_poolf > 0L) {
         if (length(ijpwef) > 0L) {
            # derivation of pool weights
            dpw_dpf$v=as.double(mx[ijpwef[,1L]]*dpwe)
            mpf[]=meas2sum%stm%dpw_dpf
            # growth flux depending on free pools
            if (nb_fgr > 0L) {
               mpf=mpf+as.matrix(mffg[,nb_ff+seq_len(nb_fgr),drop=F])
               mff=mffg[,seq_len(nb_ff)]
            } else {
               mff=mffg
            }
         }
      }
      mff=as.matrix(mff)
      if (nb_sc > 0) {
         vsc=c(1., param)[ir2isc]
         mff=vsc*mff
         mpf=vsc*mpf
      }
      # store usefull information in global list jx_f
      dux_dp[, seq_len(nb_ff)]=mff
      dux_dp[, nb_ff+seq_len(nb_sc)]=as.matrix(dur_dsc)
      dux_dp[, nb_ff+nb_sc+seq_len(nb_fgr)]=mpf
      jx_f$param=param
      jx_f$x_f=x_f
      jx_f$dux_dp=dux_dp
   }
   return(list(x=x, lf=jx_f$lf))
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
   # return the sum of bits in every component of the integer vector i
   # The result has the length=length(i)
   return(colSums(outer(2**(0:30), as.integer(i), bitwAnd) > 0))
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

cumo_gradj=function(param, labargs) {
   # calculate gradient of cost function for cumomer minimization probleme
   grad=2*as.numeric(crossprod(jx_f$res, jx_f$jacobian))
   return(grad)
}

# cost function for donlp2 solver
cumo_fn=function(p) {
   return(cumo_cost(p, labargs))
}

cumo_dfn=function(p) {
   return(cumo_gradj(p, labargs))
}

attr(cumo_fn, "gr")=cumo_dfn
#cumo_fn@gr=cumo_dfn
cumo_jacob=function(param, labargs) {
   # calculate jacobian dmeas_dparam and some annexe matrices
   # without applying invvar matrix
   # The result is in a returned list jx_f.
   jx_f=labargs$jx_f
   if (is.null(jx_f$udr_dp)) {
      jx_f$udr_dp=rbind(jx_f$dux_dp, labargs$dufm_dp, labargs$dupm_dp)
   } else {
      jx_f$udr_dp[1L:nrow(jx_f$dux_dp),]=jx_f$dux_dp
   }
   dimnames(jx_f$udr_dp)=list(labargs$nm$resid, labargs$nm$par)
   return(NULL)
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
   
   #nb_c=spAbr$nb_c # cumomer or fragment number (when emu==T)
   #w=spAbr$w # cumomer weight
   nb_c=spAbr$nb_c
   if (nb_c == 0) {
      A=simple_triplet_zero_matrix(nb_c, nb_c)
      b=simple_triplet_zero_matrix(nb_c, 1)
      return(list(A=A, b=b))
   }
   ind_a=spAbr$ind_a
   if (!is.matrix(ind_a))
      ind_a=t(ind_a)
   #x=fwrv[ind_a[,"indf"]]
   #i=which(ind_a[,"ir0"]==ind_a[,"ic0"])
   #x[i]=-x[i] # diagonal terms are negative
   #A=sparseMatrix(i=ind_a[,"ir0"]+1, j=ind_a[,"ic0"]+1, x=x, dims=c(nb_c, nb_c), giveCsparse=FALSE)
   #dimnames(A)=list(head(nm_rcumo, nb_c), head(nm_rcumo, nb_c))
   #A=Rmumps$new(ind_a[,"ir0"], ind_a[,"ic0"], x, nb_c)
   spAbr$xmat$v <- fwrv[ind_a[,"indf"]]
   x <- col_sums(spAbr$xmat)
   x[spAbr$iadiag] <- -x[spAbr$iadiag]
   spAbr$a$set_mat_data(x)
   
   # construct a complete b vector
   if (getb) {
      ind_b=if (emu) spAbr[["ind_b_emu"]] else spAbr[["ind_b"]]
      if (!is.matrix(ind_b))
         ind_b=t(ind_b)
      spAbr$bmat$v=fwrv[ind_b[,"indf"]]*incu[ind_b[,"indx1"]]*incu[ind_b[,"indx2"]]
      spAbr$b$v=-col_sums(spAbr$bmat)
      return(list(A=spAbr$a, b=spAbr$b))
      #if (emu) {
         #b_pre=spAbr$b_pre_emu
         #iwc=length(b_pre)
         #for (iwe in 1:iwc) {
         #   b_pre[[iwe]]@x=fwrv[spAbr$ind_fbe[[iwe]]]*incu[spAbr$ind_xe1[[iwe]]]*incu[spAbr$ind_xe2[[iwe]]]
         #}
         #b=-vapply(b_pre, colSums, numeric(nb_c))
      #   ind_b=spAbr[["ind_b_emu"]]
      #   x=-fwrv[ind_b[,"indf"]]*incu[ind_b[,"indx1"]]*incu[ind_b[,"indx2"]]
      #   b=sparseMatrix(i=ind_b[,"irow"], j=ind_b[,"iwe"], x=x, dims=c(nb_c, w))
      #   dimnames(b)=list(head(nm_rcumo, nb_c), 1:w)
      #} else {
      #   ind_b=spAbr[["ind_b"]]
      #   x=-fwrv[ind_b[,1]]*incu[ind_b[,2]]*incu[ind_b[,3]]
      #   b=sparseMatrix(i=ind_b[,"irow"], j=rep.int(1, nrow(ind_b)), x=x, dims=c(nb_c, 1))
      #   dimnames(b)=list(head(nm_rcumo, nb_c), 1)
      #}
      #return(list(A=spAbr$a, b=b))
   } else {
      return(list(A=spAbr$a, b=NULL))
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
   #
   # update: added vectorization for multiple incu vectors considered
   # as matrix columns (there are nco of them) . output is a matrix with dim
   # (nb_x*nco, nb_fwrv)
   # 2014-07-15 sokol
   
   # we derivate a*x=b implicitly
   # a_f*x + a*x_f=b_f + b_xl*xl_f
#browser()
   nb_c=spAbr$nb_c
   if (nb_c==0) {
      # no system at this weight
      return(list(j_rhsw=NULL, b_x=NULL, j_rhswp=NULL, b_xp=NULL))
   }
   emu=is.matrix(spAbr$ind_b_emu)
   nb_fwrv=spAbr$nb_fwrv
   nb_cl=spAbr$nb_cl
   w=spAbr$w
   
   # a_fx
   ind_a=spAbr$ind_a
   if (!is.matrix(ind_a))
      ind_a=t(ind_a)
   i=ind_a[,"ic0"]==ind_a[,"ir0"]
   incu=as.matrix(incu)
   nco=ncol(incu)
   nro=nrow(ind_a)
   emuw=ifelse(emu, w, 1L)
   nb_xw=nb_c*emuw
   # first a_fx creation for futher updates (a_fx can be of two different sizes
   # depending on the time scheme for ODE)
   first_a_fx=FALSE
   if (is.null(spAbr$a_fx1)) {
      # create auxiliary data for a_fx1
      nm_a_fx="a_fx1"
      first_a_fx=TRUE
   } else if (spAbr$a_fx1$nco != nco && is.null(spAbr$a_fx2)) {
      nm_a_fx="a_fx2"
      first_a_fx=TRUE
   } else if (spAbr$a_fx1$nco==nco) {
      nm_a_fx="a_fx1"
   } else if (spAbr$a_fx2$nco==nco) {
      nm_a_fx="a_fx2"
   } else {
      stop("Must not happen: not a_fx1 neither a_fx2 fits nco")
   }
   if (first_a_fx) {
      l=list()
      l$nco=nco
      if (emu) {
         ia=(ind_a[,"ir0"]+1L)+rep((seq_len(w)-1L)*nb_c, each=nro)
         ja=rep(ind_a[,"indf"], w)
      } else {
         ia=ind_a[,"ir0"]+1L
         ja=ind_a[,"indf"]
      }
      l$ineg=which(ind_a[,"ic0"]==ind_a[,"ir0"])
      ia=ia+rep((seq_len(nco)-1L)*(nb_c*emuw), each=nro*emuw)
      ja=rep(ja, nco)
      nar=nb_xw*nco # row number in a_fx, b_f, b_x
      iv1=ia+(ja-1)*nar
      o=order(iv1)
      l$oa_x=o
      iv1=iv1[o]
      lrep=lrepx=rle(iv1)
      lrepx$values=seq_along(lrep$values)
      l$xmat=simple_triplet_matrix(i=unlist(lapply(lrep$lengths, seq_len)),
         j=inverse.rle(lrepx), v=rep(1, length(iv1)))
      iu1=lrep$values
      i=as.integer((iu1-1)%%nar)+1L
      j=as.integer((iu1-1)%/%nar)+1L
      l$a_fx=simple_triplet_matrix(i=i, j=j, v=rep(1, length(i)), nrow=nar, ncol=nb_fwrv)
      
      # prepare b_f auxiliaries
      ind_b=if(emu) spAbr$ind_b_emu else spAbr$ind_b
      if (!is.matrix(ind_b))
         ind_b=t(ind_b)
      nro=nrow(ind_b)
      ib=ind_b[,"irow"]+if(emu) nb_c*(ind_b[,"iwe"]-1L) else 0L
      jb=ind_b[,"indf"]
      ib=ib+rep((seq_len(nco)-1L)*nb_xw, each=nro)
      jb=rep(jb, nco)
      iv1=ib+(jb-1)*nar
      o=order(iv1)
      l$ob_f=o
      iv1=iv1[o]
      lrep=lrepx=rle(iv1)
      lrepx$values=seq_along(lrep$values)
      l$b_fmat=simple_triplet_matrix(i=unlist(lapply(lrep$lengths, seq_len)),
         j=inverse.rle(lrepx), v=rep(1, length(iv1)))
      iu1=lrep$values
      i=as.integer((iu1-1)%%nar)+1L
      j=as.integer((iu1-1)%/%nar)+1L
      l$b_f=simple_triplet_matrix(i=i, j=j, v=rep(1, length(i)), nrow=nar, ncol=nb_fwrv)
      
      # prepare b_x auxiliaries
      if (all(dim(spAbr$ind_bx) > 0)) {
         ind_bx=if (emu) spAbr$ind_bx_emu else spAbr$ind_bx
         nro=nrow(ind_bx)
         ib=ind_bx[,"irow"]
         jb=ind_bx[,"ic1"]
         ib=ib+rep((seq_len(nco)-1L)*nb_xw, each=nro)
         jb=rep(jb, nco)
         iv1=ib+(jb-1)*nar
         o=order(iv1)
         l$ob_x=o
         iv1=iv1[o]
         lrep=lrepx=rle(iv1)
         lrepx$values=seq_along(lrep$values)
         l$b_xmat=simple_triplet_matrix(i=unlist(lapply(lrep$lengths, seq_len)),
            j=inverse.rle(lrepx), v=rep(1, length(iv1)))
         iu1=lrep$values
         i=as.integer((iu1-1)%%nar)+1L
         j=as.integer((iu1-1)%/%nar)+1L
         l$b_x=simple_triplet_matrix(i=i, j=j, v=rep(1, length(i)), nrow=nar, ncol=nb$nbc_x[w])
      } else {
         l$b_x=simple_triplet_zero_matrix(nrow=nar, ncol=nb$nbc_x[w])
      }
      spAbr[[nm_a_fx]]=l
   }
   x=incu[(1L+nb$xi+nb$nbc_x[w])+seq_len(nb_xw),,drop=F]
   if (emu) {
      dim(x)=c(nb_c, w, nco)
      tmp=x[ind_a[,"ic0"]+1L,,,drop=F]
      tmp[spAbr[[nm_a_fx]]$ineg,,]=-tmp[spAbr[[nm_a_fx]]$ineg,,]
      #ia=(ind_a[,"ir0"]+1L)+rep((seq_len(w)-1L)*nb_c, each=nro)
      #ja=rep(ind_a[,"indf"], w)
   } else {
      tmp=x[ind_a[,"ic0"]+1L,,drop=F]
      tmp[spAbr[[nm_a_fx]]$ineg,]=-tmp[spAbr[[nm_a_fx]]$ineg,]
      #ia=ind_a[,"ir0"]+1L
      #ja=ind_a[,"indf"]
   }
   #ia=ia+rep((seq_len(nco)-1L)*(nb_c*emuw), each=nro*emuw)
   #ja=rep(ja, nco)
   #a_fx=sparseMatrix(i=ia, j=ja, x=as.double(tmp), dims=c(nb_xw*nco, nb_fwrv))
   spAbr[[nm_a_fx]]$xmat$v=tmp[spAbr[[nm_a_fx]]$oa_x]
   spAbr[[nm_a_fx]]$a_fx$v=slam::col_sums(spAbr[[nm_a_fx]]$xmat)
   a_fx=spAbr[[nm_a_fx]]$a_fx
   
   # prepare b_f
   # NB emu: b is shorter than xw by the last M+N vector which is added as (1-sum(lighter weights))
   ind_b=if(emu) spAbr$ind_b_emu else spAbr$ind_b
   if (!is.matrix(ind_b))
      ind_b=t(ind_b)
   ##nro=nrow(ind_b)
   #ib=ind_b[,"irow"]+if(emu) nb_c*(ind_b[,"iwe"]-1L) else 0L
   #jb=ind_b[,"indf"]
   #ib=ib+rep((seq_len(nco)-1L)*nb_xw, each=nro)
   #jb=rep(jb, nco)
   #b_f=sparseMatrix(i=ib, j=jb,
   #   x=as.double(-incu[ind_b[,"indx1"],]*incu[ind_b[,"indx2"],]),
   #   dims=c(nb_xw*nco, nb_fwrv)
   spAbr[[nm_a_fx]]$b_fmat$v=(-incu[ind_b[,"indx1"],]*incu[ind_b[,"indx2"],])[spAbr[[nm_a_fx]]$ob_f]
   spAbr[[nm_a_fx]]$b_f$v=slam::col_sums(spAbr[[nm_a_fx]]$b_fmat)
   b_f=spAbr[[nm_a_fx]]$b_f
   
   # prepare b_x
   if (all(dim(spAbr$ind_bx) > 0)) {
      ind_bx=if (emu) spAbr$ind_bx_emu else spAbr$ind_bx
      #nro=nrow(ind_bx)
      #ib=ind_bx[,"irow"]
      #jb=ind_bx[,"ic1"]
      #ib=ib+rep((seq_len(nco)-1L)*nb_xw, each=nro)
      #jb=rep(jb, nco)
      #b_x=sparseMatrix(
      #   i=ib,
      #   j=jb,
      #   x=as.double(-fwrv[ind_bx[,"indf"]]*incu[ind_bx[,"indx"],]),
      #   dims=c(nb_xw*nco, nb$nbc_x[w])
      #)
      spAbr[[nm_a_fx]]$b_xmat$v=(-fwrv[ind_bx[,"indf"]]*incu[ind_bx[,"indx"],])[spAbr[[nm_a_fx]]$ob_x]
      spAbr[[nm_a_fx]]$b_x$v=slam::col_sums(spAbr[[nm_a_fx]]$b_xmat)
      b_x=spAbr[[nm_a_fx]]$b_x
   } else {
      b_x=simple_triplet_zero_matrix(nrow=nb_xw*nco, ncol=nb$nbc_x[w])
   }
   if (!is.null(incup)) {
   #   # calculate first derivative in time
   #   xp0=c(0., incup)
   #   # form a_fx_pre
   #   a_fx_pre@x=xp0[x2ta_fx[,2]]-xp0[x2ta_fx[,3]]
   #   # calculate @x slot of a_fxp
   #   a_fxp=a_fx
   #   a_fxp@x=colSums(as.matrix(a_fx_pre))

   #   # prepare b_fp
   #   b_fp=spAbr$tb_f
   #   b_fp@x=incup[x2tb_f[,2]]*incu[x2tb_f[,3]]+incu[x2tb_f[,2]]*incup[x2tb_f[,3]]

   #   # prepare b_xp
   #   b_xp=spAbr$tb_x
   #   if (nrow(b_x) > 0 && all(dim(fx2tb_x) > 0)) {
   #      # form b_x_pre
   #      b_x_pre@x=fwrv[fx2tb_x[,2]]*incup[fx2tb_x[,3]]
   #      # calculate @x slot of b_x
   #      b_xp@x=Matrix::colSums(b_x_pre)
   #   }
   #   j_rhswp=-t(b_fp+a_fxp)
   #   b_xp=-t(b_xp)
   } else {
      j_rhswp=NULL
      b_xp=NULL
   }
#browser()
   return(list(j_rhsw=stm_pm(b_f, a_fx, "-"), b_f=b_f, a_fx=a_fx, b_x=b_x, j_rhswp=j_rhswp, b_xp=b_xp))
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

df_dffp=function(param, flnx, nb_f, nm_list) {
   # derivation of fwrv by free_fluxes+poolf (and not growth fluxes neither log(poolf))
   ah=1.e-10; # a heavyside parameter to make it derivable in [-ah; ah]
   nb_fwrv=nb_f$nb_fwrv
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
   #df_dfl=Matrix(0., nb_fwrv, length(flnx))
   #df_dffd=Matrix(0., nb_fwrv, nb_ff+nb_fgr)
   df_dfl=simple_triplet_zero_matrix(nb_fwrv, length(flnx))
   df_dffd=simple_triplet_zero_matrix(nb_fwrv, nb_ff+nb_fgr)
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
   df_dfl[nb_f$cfw_fl]=c(hnet, xch)
   # reverse fluxes
   df_dfl[nb_f$crv_fl]=c(hnet-1., xch)
   
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
   if (length(nb_f$cfw_ff) > 0)
      df_dffd[nb_f$cfw_ff]=c(hnet, xch)
   # reverse fluxes
   if (length(nb_f$crv_ff) > 0)
      df_dffd[nb_f$crv_ff]=c(hnet-1., xch)
   
   # derivation by growth fluxes
   # forward fluxes
   if (length(nb_f$cfw_fg) > 0)
      df_dffd[nb_f$cfw_fg]=rep.int(1., length(i_fgn))
   # reverse fluxes
   if (length(nb_f$crv_fg) > 0)
      df_dffd[nb_f$crv_fg]=0.
   
   res=(df_dfl%stm%nb_f$dfl_dffg+df_dffd)%mrv%c(rep.int(1., nb_ff), rep(nb_f$mu, nb_fgr))
   dimnames(res)=list(nm_list$fwrv, names(param)[c(i_ffn, i_ffx, i_fgn, i_fgx)])
   return(res)
}

dufm_dff=function(nb_f, nm_list) {
   # measured fluxes derivation (non reduced by SD)
   res=matrix(0., length(nm_list$fmn), length(nm_list$ff))
   dimnames(res)=list(nm_list$fmn, nm_list$ff)
   # derivate free measured fluxes (trivial)
   i=grep("f.n.", nm_list$fmn, fixed=T)
   if (length(i) > 0) {
      res[i,nm_list$fmn[i]]=diag(length(i))
   }
   # derivate dependent measured fluxes
   i=grep("d.n.", nm_list$fmn, fixed=T, value=T)
   if (length(i) > 0) {
      res[i,]=as.matrix(nb_f$dfl_dffg[i,1:length(nm_list$ff)])
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

opt_wrapper=function(param, measurements, jx_f, trace=1) {
   oldmeas=labargs$measurements
   labargs$measurements=measurements
   labargs$jx_f=jx_f
   if (method == "BFGS") {
      control=list(maxit=500, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         labargs)
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-08,
         labargs)
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl)]=control_ftbl
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-08,
         labargs)
   } else if (method == "nlsic") {
      control=list(trace=trace, btfrac=0.25, btdesc=0.1, maxit=50, errx=1.e-5,
         ci=list(report=F), history=FALSE, adaptbt=TRUE, sln=sln,
         maxstep=max(10.*sqrt(norm2(param)), 1.)
      )
      control[names(control_ftbl)]=control_ftbl
      res=nlsic(param, lab_resid, 
         ui, ci, control, e=ep, eco=cp, flsi=lsi_fun,
         labargs)
   } else if (method == "ipopt") {
      control=list(max_iter=500, print_level=trace*5)
      control[names(control_ftbl)]=control_ftbl
      tui=c(t(ui))
      eval_g=function(x, labargs) {
         return(ui%*%x)
      }
      eval_jac_g=function(x, labargs) {
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
         labargs)
      res$par=res$solution
      names(res$par)=nm_par
   } else {
      stop_mes(paste("Unknown minimization method '", method, "'\\n", sep=""), file=fcerr)
   }
   if (is.null(res$err)) {
      res$err=0L
   }
   labargs$measurements=oldmeas
   return(res)
}

# wrapper for Monte-Carlo simulations
mc_sim=function(i, refsim, labargs=NULL) {
   #set.seed(seeds[i])
   #cat(sort(ls(pos=1)), sep="\n", file=sprintf("tmp_%d.log", i))
   for (item in c("nb_f", "measurements", "case_i")) {
      assign(item, labargs[[item]])
   }
   # random measurement generation
#cat("simlab=\n")
#print(head(as.double(refsim$simlab)))
   if (nb_f$nb_meas) {
      if (case_i) {
         meas_mc=refsim$usm+rnorm(n=length(refsim$usm))*measurements$dev$labeled
      } else {
         meas_mc=refsim$simlab+rnorm(n=length(refsim$simlab))*measurements$dev$labeled
      }
   } else {
      meas_mc=c()
   }
#cat("meas=\n")
#print(head(meas_mc))
   if (nb_f$nb_fmn) {
      fmn_mc=refsim$simfmn+rnorm(n=length(refsim$simfmn))*measurements$dev$flux
   } else {
      fmn_mc=c()
   }
   if (nb_f$nb_poolm) {
      poolm_mc=refsim$simpool+rnorm(n=length(refsim$simpool))*measurements$dev$pool
   } else {
      poolm_mc=c()
   }
#cat("imc=", i, "\\n", sep="")
#browser()
   # minimization
   measurements_mc=measurements
   if (case_i) {
      measurements_mc$vec$kin=meas_mc
   } else {
      measurements_mc$vec$labeled=meas_mc
   }
   measurements_mc$vec$flux=fmn_mc
   measurements_mc$vec$pool=poolm_mc
   loc_jx_f=new.env()
   res=opt_wrapper(param, measurements_mc, loc_jx_f, trace=0)
   #save(res, file=sprintf("mc_%d.RData", i))
   if (!is.null(res$mes) && nchar(res$mes) > 0) {
      cat((if (res$err) "Error" else "Warning"), " in Monte-Carlo i=", i, ": ", res$mes, "\n", file=fcerr, sep="")
      if (res$err) {
         res=list(cost=NA, it=res$it, normp=res$normp, par=res$par, err=res$err)
         rm(loc_jx_f)
         gc()
         return(list(cost=NA, it=res$it, normp=res$normp, par=res$par, err=res$err))
      }
   }
   # return the solution
   iva=!is.na(res$res)
   vres=res$res[iva]
   res=list(cost=crossprod(vres)[1], it=res$it, normp=res$normp, par=res$par, err=res$err)
   rm(loc_jx_f)
   gc() # for big problems we run easily out of memory
   return(res)
}

fallnx2fwrv=function(fallnx, nb_f) {
   n=length(fallnx)
   # extract and reorder in fwrv order
   net=fallnx[nb_f$inet2ifwrv]
   xch=fallnx[nb_f$ixch2ifwrv]
   # expansion 0;1 -> 0;+inf of xch (second half of fallnx)
   xch=xch/(1-xch)
   # fw=xch-min(-net,0)
   # rv=xch-min(net,0)
   fwrv=c(xch-pmin(-net,0),xch-pmin(net,0))
   return(fwrv)
}
stm_pm=function(e1, e2, pm=c("+", "-")) {
   if (pm == "-") {
      return(stm_pm(e1, -e2, "+"))
   }
   stopifnot(pm == "+")
   #ipos=match(e2$i+e2$j*e2$nrow, e1$i+e1$j*e1$nrow, nomatch=0L)
   pos=match_ij(e2$i, e2$j, e1$i, e1$j)
   #if (any(pos!=ipos)) {
   #   browser()
   #}
   ind=which(pos == 0L)
   e1$v[pos] = e1$v[pos] + e2$v[pos > 0L]
   e1$i = c(e1$i, e2$i[ind])
   e1$j = c(e1$j, e2$j[ind])
   e1$v = c(e1$v, e2$v[ind])
   return(e1)
}
