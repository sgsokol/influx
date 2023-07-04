#TIMEIT=0; # 1 to enable time printing at some stages
tmp=try(flush(fclog), silent=TRUE)
if (length(find("TIMEIT")) && TIMEIT && !inherits(tmp, "try-error")) {
   cat("load    : ", format(Sys.time()), " cpu=", proc.time()[1], "\n", sep="", file=fclog)
}
build_mult_bxxc=function(dirx) {
   fcpp="mult_bxxc.cpp"
   fso=paste("mult_bxxc", .Platform$dynlib.ext, sep="")
   if (!file.exists(file.path(dirr, "mult_bxxc.txt")) || !file.exists(file.path(dirr, fso)) ||
      file.mtime(file.path(dirr, fso)) < file.mtime(file.path(dirr, fcpp))) {
      # freshly compile the code (==first use or .so is outdated)
      frmu=file.path(system.file(package="rmumps"), "libs", .Platform$r_arch, paste("rmumps", .Platform$dynlib.ext, sep=""))
      Sys.setenv(PKG_LIBS=sprintf('"%s"', frmu))
      Sys.setenv(PKG_CXXFLAGS="-std=c++11 -fopenmp") # ‑mveclibabi=svml
      tes=capture.output(sourceCpp(file.path(dirr, "mult_bxxc.cpp"), verbose=TRUE))
      dl_str=grep("dyn.load", tes, value=TRUE)
      ftmp=sub(".*'(.*)'.*", "\\1", dl_str)
      dl_inf=sub("^(.*) <- dyn.load\\(.*$", "\\1", dl_str)
      fu_li=sub("// ", "", grep("// ", tes, value=TRUE))
      file.copy(ftmp, file.path(dirr, fso), overwrite = TRUE, copy.date=TRUE)
      sy=sapply(fu_li, function(it) {s=grep(paste(it, " <- ", sep=""), tes, value=TRUE, fixed=TRUE); sub(dl_inf, "multdll", s)})
      write.table(sy, file=file.path(dirr, "mult_bxxc.txt"), col.names=FALSE, row.names=FALSE)
   }
}
# build compiled code
#build_mult_bxxc(dirx)
#browser()
so=.Platform$dynlib.ext
#fso=paste("mult_bxxc", so, sep="")
## define R functions from mult_bxxc.so
#multdll=dyn.load(file.path(dirr, fso))
#sy=as.matrix(read.table(file=file.path(dirr, "mult_bxxc.txt"), header=FALSE))[,1]
#for (rsy in sy) {
#   eval(parse(text=rsy))
#}
#rm(multdll)

dfcg2fallnx=function(nb_f, flnx, param, fc, fg) {
   # produce complete flux (net,xch)*(dep,free,constr,growth) vector
   # from dep,free,constr,growth
   f=c(flnx[seq_len(nb_f$nb_fln)], param[seq_len(nb_f$nb_ffn)], fc[seq_len(nb_f$nb_fcn)], fg,
      flnx[nb_f$nb_fln+seq_len(nb_f$nb_flx)], param[nb_f$nb_ffn+seq_len(nb_f$nb_ffx)], fc[nb_f$nb_fcn+seq_len(nb_f$nb_fcx)], numeric(nb_f$nb_fgr))
   return(f)
}

cumo_resid=function(param, cjac=TRUE, labargs) {
   nm_local=c("jx_f", "pool", "measurements", "nm", "ir2isc", "jacobian", "nb_f", "nb_exp", "noscale")
   for (item in nm_local) {
      assign(item, labargs[[item]])
   }
   # find x for all weights and experiments
   lres=lab_sim(param, cjac, labargs)
   if (!is.null(lres$err) && lres$err) {
      return(list(err=1, mes=lres$mes))
   }

   # find simulated scaled measure vector scale*(measmat*x)
   if (is.null(jx_f$simlab)) {
      jx_f$simlab=jx_f$ureslab=jx_f$reslab=vector("list", nb_exp)
      names(jx_f$simlab)=names(jx_f$ureslab)=names(jx_f$reslab)=nm$nm_exp
   }
   for (iexp in seq_len(nb_exp)) {
      if (noscale) {
         jx_f$simlab[[iexp]]=jx_f$usimlab[[iexp]]
      } else {
         jx_f$simlab[[iexp]]=jx_f$usimlab[[iexp]]*c(1.,param)[ir2isc[[iexp]]]
      }
      names(jx_f$simlab[[iexp]])=names(jx_f$usimlab[[iexp]])=names(measurements$vec$labeled[[iexp]])
      # diff between simulated and measured
      jx_f$ureslab[[iexp]]=(jx_f$simlab[[iexp]]-measurements$vec$labeled[[iexp]])
      jx_f$reslab[[iexp]]=jx_f$ureslab[[iexp]]/measurements$dev$labeled[[iexp]]
   }

   # diff between simulated and measured
   pool[nm$poolf]=param[nm$poolf]
   jx_f$simfmn=jx_f$lf$fallnx[nm$fmn]
   jx_f$simpool=(measurements$mat$pool%*%pool)[,1]
   jx_f$uresflu=jx_f$simfmn-measurements$vec$flux
   jx_f$resflu=jx_f$uresflu/measurements$dev$flux
   jx_f$urespool=jx_f$simpool-measurements$vec$pool
   jx_f$respool=jx_f$urespool/measurements$dev$pool
   jx_f$res=c(unlist(jx_f$reslab), jx_f$resflu, jx_f$respool)
   jx_f$res[measurements$outlier]=NA
   jx_f$ures=c(unlist(jx_f$ureslab), jx_f$uresflu, jx_f$urespool)

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
      jacobian[]=as.numeric(with(measurements$dev, jx_f$udr_dp/c(unlist(labeled), flux, pool)))
      jx_f$jacobian=jacobian
      jx_f$dr_dff=jacobian[,seq_len(nb_f$nb_ff),drop=FALSE]
      labargs$jacobian=jacobian
      return(list(res=jx_f$res, jacobian=jacobian))
   } else {
      return(list(res=jx_f$res))
   }
}

cumo_cost=function(param, labargs, resl=lab_resid(param, cjac=FALSE, labargs)) {
   if (!is.null(resl$err) && resl$err) {
      return(NA)
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
   
   # local variabl assignments from labargs
   nm_local=c("nb_f", "nm", "invAfl", "p2bfl", "bp", "g2bfl", "fc")
   for (item in nm_local) {
      assign(item, labargs[[item]])
   }
   
   fg=numeric(nb_f$nb_fgr)
   names(fg)=nm$fgr
   if (nb_f$nb_fgr > 0) {
      fg[paste("g.n.", substring(nm$poolf, 4), "_gr", sep="")]=nb_f$mu*param[nm$poolf]
   }
   flnx=as.numeric(invAfl%*%(p2bfl%stm%param[seq_len(nb_f$nb_ff)]+c(bp)+g2bfl%stm%fg))
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
   
   nb_xi=nb_f$xi
   nb_ff=nb_f$nb_ff
   nb_poolf=nb_f$nb_poolf
   nb_fgr=nb_f$nb_fgr
   nb_meas=nb_f$nb_meas
   nm_exp=nm_list$nm_exp
   nm_x=nm_list$x
   nb_sc=nb_f$nb_sc
   nb_sc_tot=nb_f$nb_sc_tot
   fg=nb_f$mu*param[nm$poolf] # the same alphabetical order
   names(fg)=nm$fgr
   pool[nm$poolf]=param[nm$poolf] # inject variable pools to pool vector

   # calculate all fluxes from free fluxes
   lf=param2fl(param, labargs)
   if (is.null(jx_f$x)) {
      jx_f$x=matrix(0, nrow=sum(nb_f$x), nb_exp)
      dimnames(jx_f$x)=list(nm_x, nm_exp)
      jx_f$usimlab=vector("list", nb_exp)
      names(jx_f$usimlab)=nm_exp
   }
   x=jx_f$x
   usimlab=jx_f$usimlab

   if (is.null(labargs$incu)) {
      labargs$incu=incu=lapply(seq_len(nb_exp), function(iexp) c(1, xi[,iexp], double(nbc_x[nb_w+1L])))
      # unreduced residuals derivated by scale params
      labargs$dur_dsc=dur_dsc=lapply(seq_len(nb_exp), function(iexp) matrix(0., nrow=nb_meas[[iexp]], ncol=nb_sc_tot))
      labargs$dux_dp=dux_dp=lapply(seq_len(nb_exp), function(iexp) matrix(0., nb_meas[[iexp]], nb_ff+nb_sc_tot+nb_poolf))
   }
   Ali=lapply(seq_len(nb_w), function(iw) fwrv2Abr(lf$fwrv, spa[[iw]], NULL, NULL, getb=FALSE, emu=emu)$A)
   if (cjac) {
      # derivation of fwrv fluxes by free parameters: free fluxes+concentrations
      mdf_dffp=df_dffp(param, jx_f$lf$flnx, nb_f, nm)
      jx_f$df_dffp=mdf_dffp
      if (is.null(labargs$x_f)) {
         labargs$x_f=matrix(0., nrow=sum(nb_x), ncol=nb_ff+nb_fgr)
      }
      x_f=labargs$x_f
   } else {
      x_f=NULL
   }
   for (iexp in seq_len(nb_exp)) {
      # prepare measurement pooling operations
      pwe[[iexp]][ipwe[[iexp]]]=pool[ip2ipwe[[iexp]]]
      spwe=tapply(pwe[[iexp]], pool_factor[[iexp]], sum)
      spwe=1./spwe[nm$measmat[[iexp]]]
      pwe[[iexp]]=c(pwe[[iexp]]*spwe)
      # construct the system A*x=b from fluxes
      # and find x for every weight
      # if fj_rhs is not NULL, calculate jacobian x_f
      if (cjac) {
         if (length(ijpwef[[iexp]])) {
            dpwe=-pwe[[iexp]]*spwe
            # dpwe is shortened to non zero entries in dpw_dpf
            dpwe=ifelse(ipf_in_ppw[[iexp]][ijpwef[[iexp]][,1L]]==ijpwef[[iexp]][,2L], (spwe+dpwe)[ijpwef[[iexp]][,1L]], dpwe[ijpwef[[iexp]][,1L]])
         }
      }
      
      # simulate labeling weight by weight
      ba_x=0
      for (iw in seq_len(nb_w)) {
         nb_c=spa[[iw]]$nb_c
         emuw=ifelse(emu, iw, 1L)
         if (nb_c == 0) {
            next
         }
         ixw=nbc_x[iw]+seq_len(nb_x[iw])
         incuw=(1L+nb_xi)+ixw
#cat("iw=", iw, "\there 1\n")
#if (iw==1) {
#   print(labargs)
#   print(spa[[iw]])
#}
         if (emu) {
            b=fwrv2Abr(lf$fwrv, spa[[iw]], incu[[iexp]], nm$emu[[iexp]][ixw], getA=FALSE, emu=emu)$b
         } else {
            b=fwrv2Abr(lf$fwrv, spa[[iw]], incu[[iexp]], nm$rcumo[[iexp]][ixw], getA=FALSE, emu=emu)$b
         }
#browser()
         xw=try(solve(Ali[[iw]], b), silent=TRUE)
#print(head(xw))
         if (inherits(xw, "try-error")) {
#browser()
            rerr=attr(xw, "condition")
            if (length(grep("rmumps:.*info\\[1\\]=-10,", rerr$message, fixed=FALSE))) {
               # find 0 rows if any
               l=spa[[iw]]
               ag=aggregate(abs(lf$fwrv[l$ind_a[,"indf"]]), list(l$ind_a[,"ir0"]), sum)
               izc=ag$Group.1[ag$x <= 1.e-10]
               izf=names(which(abs(lf$fwrv)<1.e-7))
               mes=paste("Cumomer matrix is singular. Try '--clownr N' or/and '--zc N' options with small N, say 1.e-3\nor constrain some of the fluxes listed below to be non zero\n",
                  "Zero rows in cumomer matrix A at weight ", iw, ":\n",
                  paste(nm$rcumo[ixw][izc+1], collapse="\n"), "\n",
                  "Zero fluxes are:\n",
                  paste(izf, collapse="\n"), "\n",
                  sep="")
            } else {
               mes=as.character(rerr$message)
            }
            return(list(x=NULL, fA=Ali[[iw]], err=1L, mes=mes))
         }
         if (emu) {
            xw=c(xw, 1.-rowSums(xw))
         }
         incu[[iexp]][incuw]=xw
         if (cjac && nb_ff+nb_fgr > 0) {
            # calculate jacobian x_f
            # first, calculate right hand side for jacobian calculation
            # j_rhsw, b_x from sparse matrices
            # bind cumomer vector
            j_b_x=fx2jr(jx_f$lf$fwrv, spa[[iw]], nb_f, incu[[iexp]])
            j_rhsw=j_b_x$j_rhsw%stm%mdf_dffp
            b_x=j_b_x$b_x
            if (iw > 1) {
               if (ba_x > 0) {
                  bop(j_rhsw, 1, "+=", b_x%stm%x_f[seq_len(ba_x),,drop=FALSE])
               }
            }
            redim(j_rhsw, c(nb_c, emuw*(nb_ff+nb_fgr)))
            tmp=try(solve(Ali[[iw]], j_rhsw))
            if (inherits(tmp, "try-error")) {
               #browser()
               mes="Some obscure problem with label matrix.\n"
               return(list(x=NULL, fA=Ali[[iw]], err=1L, mes=mes))
            } else {
               bop(j_rhsw, 1, "=", tmp)
            }
            if (emu) {
               redim(j_rhsw, c(nb_c, iw, nb_ff+nb_fgr))
               bop(x_f, c(1, ba_x, iw*nb_c), "=", j_rhsw)
               # m+N component
               #x_f[ba_x+iw*nb_c+seq_len(nb_c),]= -apply(j_rhsw, c(1L,3L), sum)
               bop(x_f, c(1, ba_x+iw*nb_c, nb_c), "=", -arrApply(j_rhsw, 2, "sum"))
            } else {
               bop(x_f, c(1, ba_x, nb_c), "=", j_rhsw)
            }
         }
         ba_x=ba_x+nb_x[iw]
      }
   
      #rownames(incu)=c("one", nm$inp, nm$x)
      x[, iexp]=tail(incu[[iexp]], -nb_xi-1L)
      
      # calculate unreduced and unscaled measurements
      if (nrow(x) == ncol(measmat[[iexp]])) {
         mx=(measmat[[iexp]]%stm%x[, iexp])[,1]+memaone[[iexp]]
      } else {
         mx=(measmat[[iexp]]%stm%x[nm$rcumo_in_cumo, iexp])[,1]+memaone[[iexp]]
      }
      # measurement vector before pool ponderation
#browser()
      if (length(ipooled[[iexp]]) > 1L) {
         mv=meas2sum[[iexp]]%stm%(pwe[[iexp]]*mx)
      } else {
         mv=mx
      }
      jx_f$usimlab[[iexp]]=as.numeric(mv)
      if (cjac) {
         # scale part of jacobian
         if (!noscale && nb_sc[[iexp]] > 0) {
            is2m=nb_f$is2m
            dur_dsc[[iexp]][is2m[[iexp]]]=mv[is2m[[iexp]][,1]]
         }
   #browser()
         # free flux part of jacobian (and partially free pools if present in x_f)
         if (nb_ff+nb_fgr > 0) {
            mffg=measmat[[iexp]]%stm%x_f
            if (length(ipooled[[iexp]]) > 1L) {
               mffg=meas2sum[[iexp]]%stm%(pwe[[iexp]]*mffg)
            }
         } else {
            mffg=matrix(0., nrow=nb_meas[[iexp]], ncol=0L)
         }
         # free pool part of jacobian
         mff=mffg
         mpf=matrix(0., nrow=nb_meas[[iexp]], ncol=nb_f$nb_poolf)
         if (nb_f$nb_poolf > 0L) {
            if (length(ijpwef[[iexp]]) > 0L) {
               # derivation of pool weights
               dpw_dpf[[iexp]]$v=as.double(mx[ijpwef[[iexp]][,1L]]*dpwe)
               mpf[]=meas2sum[[iexp]]%stm%dpw_dpf[[iexp]]
               # growth flux depending on free pools
               if (nb_fgr > 0L) {
                  bop(mpf, 1, "+=", as.matrix(mffg[,nb_ff+seq_len(nb_fgr),drop=FALSE]))
                  mff=mffg[,seq_len(nb_ff)]
               } else {
                  mff=mffg
               }
            }
         }
         mff=as.matrix(mff)
         if (!noscale) {
            vsc=c(1., param)[ir2isc[[iexp]]]
            mff=vsc*mff
            mpf=vsc*mpf
         }
         # store usefull information in global list jx_f
         bop(dux_dp[[iexp]], c(2, 0, nb_ff), "=", mff)
         bop(dux_dp[[iexp]], c(2, nb_ff, nb_sc_tot), "=", dur_dsc[[iexp]])
         bop(dux_dp[[iexp]], c(2, nb_ff+nb_sc_tot, nb_f$nb_poolf), "=", mpf)
         jx_f$param=param
         jx_f$x_f[[iexp]]=x_f
         jx_f$dux_dp[[iexp]]=dux_dp[[iexp]]
      }
   }
   jx_f$x=x
   return(list(x=jx_f$x, lf=jx_f$lf))
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
   # or Metab:N+m, where m is emu weight M+m .
   
   if (length(x)==0) {
      return(NULL)
   }
   # prepare x as matrix
   nb_exp=NA
   if (is.vector(x) && !is.list(x)) {
      nm_x=names(x)
      x=as.matrix(x)
   } else if (is.matrix(x)) {
      nm_x=rownames(x)
   } else if (is.list(x)) {
      nm_x=rownames(x[[1]])
      nb_exp=length(x)
      nm_exp=names(x)
   } else {
      stop("x is an unknown structure. It must be vector, matrix or a list of matrices")
   }
   
   # is it emu or cumomer vector
   emu=TRUE
   if (is.na(strsplit(nm_x[1], emusep, fixed=TRUE)[[1]][2])) {
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
      if (is.na(nb_exp)) {
         res=x[unlist(sapply(nm_l, function(nm) which(spl[,1]==nm&spl[,2]==longest[nm]))),,drop=FALSE]
      } else {
         i=unlist(sapply(nm_l, function(nm) which(spl[,1]==nm&spl[,2]==longest[nm])))
         res=lapply(x, function(xx) xx[i,,drop=FALSE])
      }
      return(res)
   }

   # separate cumos by name and order by weight
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
   if (is.na(nb_exp)) {
      x=x[i,,drop=FALSE]
   } else {
      x=lapply(x, "[", i,,drop=FALSE)
   }
   spl=spl[,i,drop=FALSE]
   n=if (is.na(nb_exp)) nrow(x) else nrow(x[[1]])
   i=seq_len(n)
   icumo=as.integer(spl[2,])
   metabs=spl[1,]
   umetabs=union(metabs, NULL)
#cat("metabs:\n")
#print(metabs)
#cat("tbl:\n")
#print(tbl)
   # extract, order and convert each metab vector
   if (is.na(nb_exp)) {
      res=matrix(0., nrow=0, ncol=ncol(x))
      for (metab in umetabs) {
#cat(paste(metab,":\n",sep=""))
         im=metabs==metab
#print(d)
         o=order(icumo[im])
         # ordered cumomer vector with #0==1 component
         vcumo=rbind(1,x[im,,drop=FALSE][o,,drop=FALSE])
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
   } else {
      rres=lapply(seq_len(nb_exp), function(iexp) matrix(0., nrow=0, ncol=ncol(x[[iexp]])))
      names(rres)=names(x)
      for (iexp in seq_len(nb_exp)) {
         res=rres[[iexp]]
         for (metab in umetabs) {
            im=metabs==metab
            o=order(icumo[im])
            # ordered cumomer vector with #0==1 component
            vcumo=rbind(1,x[[iexp]][im,,drop=FALSE][o,,drop=FALSE])
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
         rres[[iexp]]=res
      }
      res=rres
   }
   return(res)
}

cumo2lab=function(x) {
   # converts cumomer vector to fraction of labeled isotopomer 1-i#0
   # separate cumos by name and order by weight
   x=as.matrix(x)
   n=nrow(x)
   nm_x=rownames(x)
   if (is.null(nm_x)) {
      return(NULL)
   }
   metabs=c(); # unique metab names
   spl=unlist(strsplit(nm_x,":",fix=TRUE))
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
      o=order(icumo[im])
      # ordered cumomer matrix with #0==1 component
      vcumo=rbind(1,x[im,,drop=FALSE][o,,drop=FALSE])
      clen=log2(nrow(vcumo))
      # labeled fraction
      lab=1-Vcumo2iso0(clen)%*%vcumo
      rownames(lab)=metab
      res=rbind(res, lab)
   }
   return(res)
}
cumo2iso=function(x) {
   # converts cumomer vector to isotopomer vector
   # separate cumos by name and order by weight
   x=as.matrix(x)
   n=nrow(x)
   nm_x=rownames(x)
   if (is.null(nm_x)) {
      return(NULL)
   }
   metabs=c(); # unique metab names
   spl=unlist(strsplit(nm_x,":",fix=TRUE))
   i=1:n
   icumo=as.integer(spl[2*i])
   metabs=spl[2*i-1]
   umetabs=union(metabs, NULL)
   # extract, order and convert each metab vector
   res=c()
   for (metab in umetabs) {
      im=metabs==metab
      o=order(icumo[im])
      # ordered cumomer matrix with #0==1 component
      vcumo=rbind(1,x[im,,drop=FALSE][o,,drop=FALSE])
      clen=round(log2(nrow(vcumo)))
      # labeled fraction
      lab=Tcumo2iso(clen)%*%vcumo
      rownames(lab)=paste(metab, seq_len(2**clen)-1, sep=":")
      res=rbind(res, lab)
   }
   return(res)
}
cumo_gradj=function(param, labargs) {
   # calculate gradient of cost function for cumomer minimization probleme
   if (any(ina <- is.na(jx_f$res)))
      grad=2*as.numeric(crossprod(jx_f$res[!ina], jx_f$jacobian[!ina,,drop=FALSE]))
   else
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
   nb_exp=labargs$nb_exp
   nblr=sapply(jx_f$dux_dp, nrow) # labeled residual lengths
   nblr_tot=sum(nblr)
   tmp=double(nblr_tot*length(param))
   redim(tmp, c(nblr_tot, length(param)))
   base=0
   for (iexp in seq_len(nb_exp)) {
      tmp[base+seq_len(nblr[iexp]),]=jx_f$dux_dp[[iexp]]
      base=base+nblr[iexp]
   }
   if (is.null(jx_f$udr_dp)) {
      jx_f$udr_dp=rbind(tmp, labargs$dufm_dp, labargs$dupm_dp)
   } else {
      jx_f$udr_dp[seq_len(nblr_tot),]=tmp
   }
   dimnames(jx_f$udr_dp)=list(labargs$nm$resid, labargs$nm$par)
   return(NULL)
}

fwrv2Abr=function(fwrv, spAbr, incu, nm_rcumo, getA=TRUE, getb=TRUE, emu=FALSE) {
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
#cat("fwrv2Abr 1\n")
   nb_c=spAbr$nb_c
   if (nb_c == 0) {
      A=spAbr$a
      b=spAbr$b
      return(list(A=if (getA) A else NULL, b=if (getb) b else NULL))
   }
   if (getA) {
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
#print(spAbr)
      spAbr$a$set_mat_data(x)
   }
   
   # construct a complete b vector
   if (getb) {
      ind_b=if (emu) spAbr[["ind_b_emu"]] else spAbr[["ind_b"]]
#      if (!is.matrix(ind_b))
#         ind_b=t(ind_b)
      nprodx=ncol(ind_b)-2-emu
      prodx=incu[c(ind_b[,2+emu+seq_len(nprodx)])]
      dim(prodx)=c(nrow(ind_b), nprodx)
      spAbr$bmat$v=fwrv[ind_b[,"indf"]]*arrApply(prodx, 2, "prod")
      spAbr$b$v=-col_sums(spAbr$bmat)
   }
   return(list(A=if (getA) spAbr$a else NULL, b=if (getb) spAbr$b else NULL))
}

fx2jr=function(fwrv, spAbr, nb, incu) {
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
   # first a_fx creation for further updates (a_fx can be of many different sizes
   # depending on the time scheme for ODE and possible parallel experiments)
   nm_a_fx=as.character(nco)
   if (is.null(spAbr$a_fxx)) {
      # create auxiliary data for a_fx
      spAbr$a_fxx=list()
   }   
   if (is.null(spAbr$a_fxx[[nm_a_fx]])) {
#cat("build", nm_a_fx, "\n", file=fclog)
      l=new.env()
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
      # cache ind and pos for (b_f-a_fx)
      l$bma_pos=if (l$b_f$nrow*l$b_f$ncol < 2251799813685248 && l$a_fx$nrow*l$a_fx$ncol < 2251799813685248) match(l$a_fx$i+l$a_fx$j*l$a_fx$nrow, l$b_f$i+l$b_f$j*l$b_f$nrow, nomatch=0L) else match_ij(l$a_fx$i, l$a_fx$j, l$b_f$i, l$b_f$j)
      l$bma_ind=which(l$bma_pos == 0L)
      
      # prepare b_x auxiliaries
      if (length(spAbr$ind_bx) > 0) {
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
      spAbr$a_fxx[[nm_a_fx]]=l
   }
   x=incu[(1L+nb$xi+nb$nbc_x[w])+seq_len(nb_xw),,drop=FALSE]
   aux=spAbr$a_fxx[[nm_a_fx]]
   if (emu) {
      redim(x, c(nb_c, w, nco))
      tmp=x[ind_a[,"ic0"]+1L,,,drop=FALSE]
      tmp[aux$ineg,,]=-tmp[aux$ineg,,]
   } else {
      tmp=x[ind_a[,"ic0"]+1L,,drop=FALSE]
      tmp[aux$ineg,]=-tmp[aux$ineg,]
   }
#browser()
   aux$xmat$v[]=tmp[aux$oa_x]
   aux$a_fx$v[]=slam::col_sums(aux$xmat)
   a_fx=aux$a_fx
   
   # prepare b_f
   # NB emu: b is shorter than xw by the last M+N vector which is added as (1-sum(lighter weights))
   ind_b=if(emu) spAbr$ind_b_emu else spAbr$ind_b
   nprodx=ncol(ind_b)-2-emu
   prodx=incu[c(ind_b[,2+emu+seq_len(nprodx)]),]
   dim(prodx)=c(nrow(ind_b), nprodx, nco)
   aux$b_fmat$v[]=-arrApply(prodx, 2, "prod")[aux$ob_f]
   aux$b_f$v[]=slam::col_sums(aux$b_fmat)
   b_f=aux$b_f
   
   # prepare b_x
   if (length(spAbr$ind_bx) > 0) {
      ind_bx=if (emu) spAbr$ind_bx_emu else spAbr$ind_bx
      nprodx=ncol(ind_bx)-3-emu
      if (nprodx >= 1) {
         prodx=incu[c(ind_bx[,3+emu+seq_len(nprodx)]),]
         dim(prodx)=c(nrow(ind_bx), nprodx, nco)
         aux$b_xmat$v[]=(-fwrv[ind_bx[,"indf"]]*arrApply(prodx, 2, "prod"))[aux$ob_x]
      } else {
         aux$b_xmat$v[]=(-fwrv[ind_bx[,"indf"]])[aux$ob_x]
      }
      aux$b_x$v[]=slam::col_sums(aux$b_xmat)
      b_x=aux$b_x
   } else {
      b_x=simple_triplet_zero_matrix(nrow=nb_xw*nco, ncol=nb$nbc_x[w])
   }
#browser()
   j_rhsw=stm_pm(b_f, a_fx, "-", aux$bma_pos, aux$bma_ind)
   if (is.null(aux$o_j)) {
      o=order(j_rhsw$i+j_rhsw$nrow*j_rhsw$j)
      aux$o_j=o
   }
   o=aux$o_j
   j_rhsw$i[]=j_rhsw$i[o]
   j_rhsw$j[]=j_rhsw$j[o]
   j_rhsw$v[]=j_rhsw$v[o]
   return(list(j_rhsw=j_rhsw, b_x=b_x))
}

put_inside=function(param, ui, ci, tol_in=1.e-10, tol_out=1.e-7, tol_desc=1.e-3) {
   # put param inside of feasible domain delimited by u%*%param >= ci
   nm_par=names(param)
   mes=""
   ineq=as.numeric(ui%*%param)-ci
   if (all(ineq>tol_in)) {
      # nothing to do, already inside and well inside
      return(param)
   }
   dp=ldp(as.matrix(ui), -ineq)
   if (!is.null(dp)) {
      # get new active inequalities
      ineqd=as.numeric(ui%*%(param+dp))-ci
      # check that we are not too far outside
      if (any(ineqd < -tol_out)) {
         param=NA
         attr(param, "mes")="Inequality system is ill-conditionned. Failed to solve."
         attr(param, "err")=1
         return(param)
      }
      iact=ineqd<=tol_in
#print(ineqd[iact])
      # solve an ldp pb to find non decreasing direction
      # for active inequalities
      ma=ui[iact,,drop=FALSE]
      na=sum(iact)
      # find closest vector to c(1,1,...) making the direction increasing
      tma=tcrossprod(ma)
      bet=ldp(tma, tol_desc - apply(tma, 1, sum))
      if (is.null(bet)) {
         param=param+dp
         attr(param, "mes")="Infeasible constraints for inside direction."
         attr(param, "err")=0
         return(param)
      }
      vn=crossprod(ma, 1.+bet)
      vn=vn/norm(vn)
      decr=(ui%*%vn)<0.
      alpha=((-ineqd)/(ui%*%vn))[decr]
      alpha=alpha[alpha>0]
      alpha=0.5*min(tol_desc, alpha)
      dpn=dp+alpha*vn
      # check that new dp is still inside
      if (any(ui%*%(param+dpn)-ci < 0.)) {
         attr(param, "err")=0 # just a warning
         attr(param, "mes")="Failed to put free parameters strictly inside of the feasible domain. They are left on the border."
         dpn=dp
      }
      names(param)=nm_par
      if (!is.null(mget("nb_ff", ifnotfound=list(NULL))[[1L]]) && nb_ff > 0) {
         i=abs(dpn[seq_len(nb_ff)])>=tol_in
         if (any(i)) {
            tmp=cbind(param[1:nb_ff], param[1:nb_ff]+dpn[1:nb_ff], dpn[1:nb_ff])
            dimnames(tmp)=list(nm_par[1:nb_ff], c("outside", "inside", "delta"))
            obj2kvh(tmp[i,,drop=FALSE], "Free fluxes put inside of feasible domain")
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
   
   res=df_dfl%stm%nb_f$dfl_dffg+df_dffd
   if (nb_fgr > 0) {
      i=res$j > nb_ff
      res$v[i]=res$v[i]*nb_f$mu
   }
   dimnames(res)=list(nm_list$fwrv, names(param)[c(i_ffn, i_ffx, i_fgn, i_fgx)])
   o=order(res$i+res$nrow*res$j)
   res$i[]=res$i[o]
   res$j[]=res$j[o]
   res$v[]=res$v[o]
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
      x=x[o,,drop=FALSE]
      if (!is.null(m) && nrow(x)==nrow(m)) {
         m=m[o,,drop=FALSE]
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
   return(list(ti=ti, usm=d[o,,drop=FALSE]))
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
   # 2016-09-25 sokol: adapted for arbitrary reaction length
   nw=length(spr)
   spemu=spr
   if (nw < 1) {
      return(spemu)
   }
   #x2tb_f=spr[[1]]$x2tb_f
   nme2iemu=seq_along(nm_inemu)
   names(nme2iemu)=nm_inemu
   ba_e=length(nm_inemu)
   
   # cumo weights are allways start from 1 to weigth max
   # if there is only one weight (max), the lower weights must be present
   # but corresponding matrices and vectors are of dim=0

   # for iwc==1 emu+M1 are identical to cumomers
   # and the system A*x=b for emu+M0 does not change as
   # all fluxes in A and b sum to 0.
   for (iwc in seq_len(nw)) {
      # iwc is the resulting fragment length
      sp=spr[[iwc]]
      spemu[[iwc]]=sp
      nb_c=sp$nb_c
      if (nb_c == 0) {
         next
      }
      ind_b=sp$ind_b
      nb_ind=nrow(ind_b)
      nprodx=ncol(ind_b)-2
      nm_c=nm_incu[c(ind_b[,2+seq_len(nprodx)])]
      iprodx=seq_len(nprodx)
      if (iwc > 1 && nprodx > 1) {
         ba_e=1+nb$xiemu
         # prepare names
         dim(nm_c)=c(nb_ind, nprodx)
         # get fragment length for each ind_x which is product of several terms
         flen=vapply(strsplit(nm_c, ":"), function(v) if (length(v) > 1) sumbit(as.numeric(v[2])) else 0, 1)
         dim(flen)=dim(nm_c)
         wid=apply(flen, 1, paste0, collapse=",")
         wdisp=list() # weight dispather helper
         for (ir in seq_len(nb_ind)) {
            v=flen[ir,]
            if (!is.null(wdisp[[wid[ir]]]))
               next
            m=Reduce(`%m+%`, v);
            wdisp[[wid[ir]]]=lapply(seq_len(iwc-1), function(i) which(m==i, arr.ind=TRUE)-1)
         }
      }
     
      for (iwe in seq_len(iwc)) {
         # iwe runs from m+0 to m+(iwc-1) to form all masses but last for
         # the current fragment length.
         if (iwe == 1 || nprodx == 1) {
            # For m+0 (iwe=1) vector b is the same in cumo and emu
            ind_b_emu=cbind(iwe=1, ind_b)
            ind_b_emu[,3+iprodx]=nme2iemu[paste(nm_c, iwe-1, sep="+")]
         } else {
            for (ir in seq_len(nb_ind)) {
#if (ir==38 && iwc==2 && iwe==2)
#   browser()
               addw=wdisp[[wid[ir]]][[iwe-1]]
               ie=nme2iemu[paste(rep(nm_c[ir,], each=nrow(addw)), addw, sep="+")] # row emu names
               dim(ie)=c(nrow(addw), nprodx)
               ind_b_emu=rbind(ind_b_emu, cbind(iwe, ind_b[ir,1], ind_b[ir,2], ie))
            }
         }
      }
      ind_b_emu[is.na(ind_b_emu)]=1 # ones stay ones
      spemu[[iwc]]$ind_b_emu=ind_b_emu
      # prepare b_x
      ind_bx=c()
      if (length(sp[["ind_bx"]]) > 0) {
         for (ix in iprodx) {
            i=ind_b_emu[,3+ix] > ba_e # exclude from differentiation plain input entries
            tmp=ind_b_emu[i,,drop=FALSE]
            ind_bx=rbind(ind_bx, tmp[,c(1:3,ix+3,3+iprodx[-ix])]) # move diff var to ic1 place
         }
         if (length(ind_bx)) {
            colnames(ind_bx)=c("iwe", "indf", "irow", "ic1", sprintf("indx%d", seq_len(nprodx-1)))
            ind_bx[,"ic1"]=ind_bx[,"ic1"]-ba_e
            ind_bx[,"irow"]=ind_bx[,"irow"]+(ind_bx[,"iwe"]-1)*nb_c
         }
      }
      spemu[[iwc]]$ind_bx_emu=ind_bx
      spemu[[iwc]]$nb_emul=sum(head(nb$emus, iwc-1))
   }
   return(spemu)
}

opt_wrapper=function(param, method, measurements, jx_f, labargs, trace=1) {
   oldmeas=labargs$measurements
   labargs$measurements=measurements
   labargs$jx_f=jx_f
   if (method == "BFGS") {
#browser()
      control=list(maxit=500, trace=trace)
      control[names(control_ftbl$BFGS)]=control_ftbl$BFGS
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-5, control,
         method="BFGS", outer.iterations=100, outer.eps=1e-08,
         labargs)
   } else if (method == "Nelder-Mead") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl[["Nelder-Mead"]])]=control_ftbl[["Nelder-Mead"]]
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="Nelder-Mead", outer.iterations=100, outer.eps=1e-08,
         labargs)
   } else if (method == "SANN") {
      control=list(maxit=1000, trace=trace)
      control[names(control_ftbl$sann)]=control_ftbl$sann
      res=constrOptim(param, cumo_cost, grad=cumo_gradj,
         ui, ci, mu = 1e-4, control,
         method="SANN", outer.iterations=100, outer.eps=1e-08,
         labargs)
   } else if (method == "nlsic") {
      control=list(trace=trace, btfrac=0.25, btdesc=0.1, maxit=50, errx=1.e-5,
         ci=list(report=F), history=FALSE, adaptbt=TRUE, sln=sln,
         maxstep=max(10.*sqrt(norm2(param)), 1.)
      )
      control[names(control_ftbl$default)]=control_ftbl$default
      control[names(control_ftbl$nlsic)]=control_ftbl$nlsic
      res=nlsic(param, lab_resid, 
         ui, ci, control, e=ep, eco=cp, flsi=lsi_fun,
         labargs)
   } else if (method == "ipopt") {
      control=list(max_iter=500, print_level=trace*5)
      control[names(control_ftbl$ipopt)]=control_ftbl$ipopt
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
   } else if (method == "pso") {
      control=list(
         type="SPSO2011",
         trace=trace,
         maxit=100,
         reltol=1.e-2
      )
      control[names(control_ftbl$pso)]=control_ftbl$pso
#print(control_ftbl)
#print(control)
      res=try(psoptim_ic(param, cumo_cost, labargs, mean=0.5, control=control, ui=ui, ci=ci), silent=TRUE)
#print(res)
      if (inherits(res, "try-error")) {
         res=list(err=1L, par=NULL, mes=attr(res, "condition")$message)
      } else {
         tmp=list(err=res$msgcode, par=res$par, mes=res$msg)
         res[c("msg", "par", "msg")]=NULL
         res=c(tmp, res) # preserve the rest of the fields: stats etc.
      }
   } else {
      stop_mes("Unknown minimization method '", method, "'", file=fcerr)
   }
   if (is.null(res$err)) {
      res$err=0L
   }
   labargs$measurements=oldmeas
   return(res)
}

# wrapper for Monte-Carlo simulations
mc_sim=function(imc) {
#cat("mc_sim imc=", imc, "\n")
#print(labargs$spa[[1]])
   #set.seed(seeds[imc])
   #cat(sort(ls(pos=1)), sep="\n", file=sprintf("tmp_%d.log", imc))
   labargs=get("labargs", envir=.GlobalEnv)
   for (item in c("nb_f", "measurements", "case_i", "dirres", "baseshort", "nb_exp")) {
      assign(item, labargs[[item]])
   }
   # random measurement generation
#cat("simlab=\n")
#print(head(as.double(refsim$simlab)))
   
   if (case_i) {
      meas_mc=lapply(seq_len(nb_exp), function(iexp) if (nb_f$nb_meas[iexp]) refsim$usm[[iexp]]+rnorm(n=length(refsim$usm[[iexp]]))*measurements$dev$labeled[[iexp]] else NULL)
   } else {
      meas_mc=lapply(seq_len(nb_exp), function(iexp) if (nb_f$nb_meas[iexp]) refsim$simlab[[iexp]]+rnorm(n=length(refsim$simlab[[iexp]]))*measurements$dev$labeled[[iexp]] else NULL)
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
#cat("imc=", imc, "\\n", sep="")
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
#rres=lab_resid(param, TRUE, labargs)
#d=diag(qr(rres$jacobian, LAPACK=TRUE)$qr)
#cat("param=", param, "\n", sep=" ")
#cat("d=", d, "\n", sep=" ")
   res=opt_wrapper(param, tail(methods, 1L), measurements_mc, loc_jx_f, labargs, trace=0)
   #save(res, file=sprintf("mc_%d.RData", imc))
   if (res$err && !is.null(res$mes) && nchar(res$mes) > 0) {
   #if (TRUE) {
      fclog=file(file.path(dirres, "tmp", sprintf("%s.%smc%d.log", baseshort, runsuf, imc)), "wb")
      cat((if (res$err) "***Error" else "***Warning"), " in Monte-Carlo i=", imc, ": ", res$mes, "\n", file=fclog, sep="")
      close(fclog)
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
cl_worker=function(funth=NULL, argth=NULL) {
   if ("labargs" %in% names(argth))
      labargs=argth$labargs
   tryCatch({
#cat("cl_worker imc=", imc, "\tmem=", sum(memuse()), "\n", sep="")
#print(ls(labargs$spa[[1]]))
      #if (is.null(labargs$spa[[1]]$a)) {
      #   labargs$spa=sparse2spa(labargs$spa)
#print(labargs$spa[[1]])
      #}
#print(labargs$spa[[1]]$a)
#browser()
      do.call(funth, argth)
#cat("i=", i, "\tmem=", sum(memuse()), "\n", sep="")
   },
   error=function(e) {
      traceback()
      print(e)
      stop(e)
   })
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
stm_pm=function(e1, e2, pm=c("+", "-"), pos=if (e1$nrow*e1$ncol < 2251799813685248 && e2$nrow*e2$ncol < 2251799813685248) match(e2$i+e2$j*e2$nrow, e1$i+e1$j*e1$nrow, nomatch=0L) else match_ij(e2$i, e2$j, e1$i, e1$j), ind=which(pos == 0L)) {
   # 2**51 is the max matrix size that fits the length of double mantissa (53 bits))
   if (pm == "+") {
      e1$v[pos] = e1$v[pos] + e2$v[pos > 0L]
      e1$v = c(e1$v, e2$v[ind])
   } else {
      e1$v[pos] = e1$v[pos] - e2$v[pos > 0L]
      e1$v = c(e1$v, -e2$v[ind])
   }
   e1$i = c(e1$i, e2$i[ind])
   e1$j = c(e1$j, e2$j[ind])
   return(e1)
}
# reorder indexes to accelerate sparse matrix construction
sparse2spa=function(spa) {
   for (ispa in seq_along(spa)) {
      l=spa[[ispa]];
      if (l$nb_c == 0) {
         l$a=Rmumps$new(integer(0), integer(0), double(0), l$nb_c)
         l$b=simple_triplet_matrix(i=integer(0), j=integer(0), v=double(0), nrow=l$nb_c, ncol=1)
         next
      }
      l$l=l;
      with(l, {
#browser();
      if (!is.matrix(ind_a))
         ind_a=t(ind_a)
      if (!is.matrix(ind_b))
         ind_b=t(ind_b)
      # prepare sparse xmat where col_sums(xmat) will give a$v
      iv0=ind_a[,"ir0"]+ind_a[,"ic0"]*nb_c
      o=order(iv0)
      ind_a=ind_a[o,]
      iv0=iv0[o]
      lrep=lrepx=rle(iv0)
      lrepx$values=seq_along(lrep$values)
      xmat=simple_triplet_matrix(i=unlist(lapply(lrep$lengths, seq_len)),
         j=inverse.rle(lrepx), v=rep(1, length(iv0)))
      iu0=lrep$values
      i=as.integer(iu0%%nb_c)
      j=as.integer(iu0%/%nb_c)
      l$a=Rmumps$new(i, j, rep(pi, length(iu0)), nb_c)
      if (!is.null(control_ftbl$mumps)) {
#browser()
         lapply(grep("^icntl_", names(control_ftbl$mumps), v=TRUE), function(nm) {
            i=suppressWarnings(as.integer(strsplit(nm, "_")[[1L]][2]))
            v=suppressWarnings(as.integer(control_ftbl$mumps[[nm]]))
            if (!is.na(i) && !is.na(v))
               l$a$set_icntl(v, i)
         })
      }
      l$iadiag=which(i==j)
      # prepare sparse bmat where col_sums(bmat) will give b$v
      if (emu) {
         if (!is.matrix(ind_b_emu))
            ind_b_emu=t(ind_b_b_emu)
         iv0=ind_b_emu[,"irow"]+(ind_b_emu[,"iwe"]-1)*nb_c-1
         nb_bcol=w
      } else {
         iv0=ind_b[,"irow"]-1
         nb_bcol=1
      }
      o=order(iv0)
      if (emu) {
         ind_b_emu=ind_b_emu[o,,drop=FALSE]
      } else {
         ind_b=ind_b[o,,drop=FALSE]
      }
#browser()
      iv0=iv0[o]
      lrep=lrepx=rle(iv0)
      lrepx$values=seq_along(lrep$values)
      bmat=simple_triplet_matrix(i=unlist(lapply(lrep$lengths, seq_len)),
         j=inverse.rle(lrepx), v=rep(1, length(iv0)))
      iu0=lrep$values
      i=as.integer(iu0%%nb_c)
      j=as.integer(iu0%/%nb_c)
      l$b=simple_triplet_matrix(i=i+1, j=j+1, v=rep(pi, length(iu0)), nrow=nb_c, ncol=nb_bcol)
   })}
   return(spa)
}
`%m+%`=function(m, n) outer(if (is.array(m)) m else seq(0, m), seq(0, n), "+")
