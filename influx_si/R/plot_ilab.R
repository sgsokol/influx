plot_mti=function(ti, x, m=NULL, ...) {
   # plot MS time curse curves x[icurve, itime] and points from m
   # x and m are supposed to have the same dimension (2) and organization
   nm=if (is.null(m)) rownames(x) else rownames(m)
   nb_curve=nrow(x)
   if (is.null(nm)) {
      nb_curve=nrow(x)
      nm=seq_len(nb_curve)
   } else {
      nb_curve=if (is.null(m)) nrow(x) else nrow(m)
      # strip the ftbl row number
      nm=sub(":[-0-9]+$", "", nm)
      # make names look like M_0, M_1 etc.
      if (length((nms <- strsplit(nm, ":", fixed=TRUE))[[1L]]) > 2L) {
         nm_leg=as.expression(sapply(nms, function(v) substitute(M[i], list(i=v[4]))))
      } else {
         nm_leg=as.expression(sapply(seq_len(nb_curve)-1, function(v) substitute(M[i], list(i=v))))
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
   plot(range(ti, tim), range(c(x,m[inna])), t="n", ylab="Label", xlab="Time", ...)
   if (!is.null(x))
      matplot(ti, t(x), t="l", lty=seq_len(nb_curve), col=seq_len(nb_curve), lwd=2, add=TRUE, ...)
   if (!is.null(m)) {
      # plot measured points
      for (i in seq_len(nrow(m))) {
         inna=which(!is.na(m[i,]))
         if (length(inna) == 0)
            next
         points(tim[inna], m[i,inna], col=i, cex=0.5, t="b", lwd=0.5, ...)
         if (nrow(m) == nb_curve && !is.null(x)) {
            # draw filled polygons between simulation and data
            tmax=max(tim[inna])
            iti=which(ti <= tmax)
            polygon(c(ti[iti],rev(tim[inna])), c(x[i,iti], rev(m[i,inna])), col=rgb(red=t(col2rgb(i)), alpha=31, max=255), bord=FALSE, ...)
         }
      }
   }
   legend("topright", legend=nm_leg, lty=1:nb_curve, col=1:nb_curve, lwd=2, cex=0.75, bg=rgb(0.97,0.97,0.97, 0.75))
}
plot_lti=function(ti, x, m=NULL, ...) {
   # plot label (other than MS) time curse curves x[icurve, itime] and points from m
   # x and m are supposed to have the same dimension (2) and organization
   nm=if (is.null(m)) rownames(x) else rownames(m)
   nb_curve=nrow(x)
   if (is.null(nm)) {
      nb_curve=nrow(x)
      nm=seq_len(nb_curve)
   } else {
      nb_curve=if (is.null(m)) nrow(x) else nrow(m)
      # strip the ftbl row number
      nm=sub(":[-0-9]+$", "", nm)
      # make names look like l:x01, l:x10 etc.
      if (length(nms <- strsplit(nm, ":", fixed=TRUE)[[1L]]) > 2L) {
         nm_leg=paste0(nms[c(1L,3L)], collapse=":")
      } else {
         nm_leg=as.characetr(seq_len(nb_curve))
      }
   }
   # x and m may have different time moments
   if (is.null(m)) {
      tim=ti
   } else {
      tim=as.numeric(colnames(m))
   }
   plot(range(ti, tim), range(x, m, na.rm=TRUE), t="n", ylab="Label", xlab="Time", ...)
   if (!is.null(x))
      matplot(ti, t(x), t="l", lty=seq_len(nb_curve), col=seq_len(nb_curve), lwd=2, add=TRUE, ...)
   if (!is.null(m)) {
      # plot measured points
      for (i in seq_len(nrow(m))) {
         inna=which(!is.na(m[i,]))
         if (length(inna) == 0)
            next
         points(tim[inna], m[i,inna], col=i, cex=0.5, t="b", lwd=0.5, ...)
         if (nrow(m) == nb_curve && !is.null(x)) {
            # draw filled polygons between simulation and data
            tmax=max(tim[inna])
            iti=which(ti <= tmax)
            polygon(c(ti[iti],rev(tim[inna])), c(x[i,iti], rev(m[i,inna])), col=rgb(red=t(col2rgb(i)), alpha=31, max=255), bord=FALSE, ...)
         }
      }
   }
   legend("topright", legend=nm_leg, lty=1:nb_curve, col=1:nb_curve, lwd=2, cex=0.75, bg=rgb(0.97,0.97,0.97, 0.75))
}
#browser()
if (write_res) {
   for (iexp in seq_len(nb_exp)) {
      usmf=jx_f$usmf[[iexp]] # unscaled simulated measurements
      me=measurements$vec$kin[[iexp]] # measured dynamic labeling data
      if (!is.null(measurements$outlier[[iexp]])) {
         iout=measurements$outlier[[iexp]]
         iout=iout[iout <= prod(dim(me))] # only outliers of label kinetics
         me[iout]=NA # measured dynamic labeling data
      }

      # get unique MS fragment names
      nm_selm=natsort(grep("^m:", if (is.null(rownames(me))) rownames(usmf) else rownames(me), v=TRUE))
      pdf(sprintf("%s/%s.pdf", dirres, nm_exp[iexp]))
      #browser()
      if (length(nm_selm) > 0) {
         plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
         text(0.5, 0.05, lab="MS measurements", cex=2)
         nmf=unique(apply(structure(sapply(strsplit(nm_selm, ":", fixed=TRUE), "[", 1:4)[2:3,], dim=c(2L, length(nm_selm))), 2, paste0, sep="", collapse=":"))
         for (metf in nmf) {
            i=grep(sprintf("m:%s:", metf), nm_selm, fix=TRUE, v=TRUE)
            #isim=pmatch(sapply(strsplit(i, ":", fixed=TRUE), function(v) paste0(v[-length(v)], collapse=":")), rownames(usmf))
            isim=grep(sprintf("m:%s:", metf), rownames(usmf), fix=TRUE, v=TRUE)
            mf=strsplit(metf, ":", fixed=TRUE)[[1L]]
            met=mf[1L]
            fr=mf[2L]
            metlen=clen[strsplit(met, "+", fixed=TRUE)[[1]][1]]
            if (is.na(fr) || fr == paste(seq_len(metlen), collapse=",") || fr == sprintf("1~%d", metlen)) {
               mainlab=met
            } else {
               mainlab=sprintf("%s:%s", met, fr)
            }
            plot_mti(tifull[[iexp]][-1L], usmf[isim,,drop=FALSE], me[i,,drop=FALSE], main=mainlab, ylim=0:1)
         }
      }
      # plot non measured MS from mid
      if (exists("mid")) {
         nm_simm=rownames(mid[[iexp]])
         nmm=if (length(nm_selm) > 0L) unique(sapply(strsplit(nm_selm, ":", fixed=TRUE), "[", 1:4)[2,]) else character(0L)
         nmmid=if (length(nm_simm) > 0L) unique(sapply(strsplit(nm_simm, "+", fixed=TRUE), "[", 1L)) else character(0L)
         if (emu) {
            nmmid=sapply(strsplit(nmmid, ":", fixed=TRUE), "[", 1L)
         }
         nmp=natsort(setdiff(nmmid, nmm))
         if (length(nmp)) {
            plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
            text(0.5, 0.05, lab="MS simulations", cex=2)
         }
         for (met in nmp) {
            if (emu) {
               i=grep(sprintf("^%s:", met), nm_simm, v=TRUE)
               # take fragments
               fr=unique(sapply(strsplit(i, "[+:]"), "[", 2L))
               for (f in fr) {
                  i=grep(sprintf("^%s:%s\\+", met, f), nm_simm, v=TRUE)
                  fi=as.integer(f)
                  mainlab=if (fi == 2**clen[met]-1) met else sprintf("%s:%s", met, fr)
                  plot_mti(tifull[[iexp]][-1L], mid[[iexp]][i,,drop=FALSE], NULL, main=mainlab, ylim=0:1)
               }
            } else {
               i=grep(sprintf("^%s\\+", met), nm_simm, v=TRUE)
               plot_mti(tifull[[iexp]][-1L], mid[[iexp]][i,,drop=FALSE], NULL, main=met, ylim=0:1)
            }
         }
         # get unique label (!=MS) fragment names
         nm_sell=natsort(grep("^[^m]:", if (is.null(rownames(me))) rownames(usmf) else rownames(me), v=TRUE))
         if (length(nm_sell) > 0) {
            plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
            text(0.5, 0.05, lab="NMR measurements", cex=2)
            nmf=unique(apply(structure(sapply(strsplit(nm_sell, ":", fixed=TRUE), "[", 1:4)[1:3,], dim=c(3L, length(nm_sell))), 2L, paste0, sep="", collapse=":"))
            for (metf in nmf) {
               i=grep(sprintf("%s:", metf), nm_sell, fix=TRUE, v=TRUE)
               #isim=pmatch(sapply(strsplit(i, ":", fixed=TRUE), function(v) paste0(v[-length(v)], collapse=":")), rownames(usmf))
               isim=grep(sprintf("%s:", metf), rownames(usmf), fix=TRUE, v=TRUE)
               mf=strsplit(metf, ":", fixed=TRUE)[[1L]]
               met=mf[2L]
               fr=mf[3L]
               metlen=clen[strsplit(met, "+", fixed=TRUE)[[1]][1]]
               if (is.na(fr) || fr == paste(seq_len(metlen), collapse=",") || fr == sprintf("1~%d", metlen)) {
                  mainlab=met
               } else {
                  mainlab=sprintf("%s:%s", met, fr)
               }
               plot_lti(tifull[[iexp]][-1L], usmf[isim,,drop=FALSE], me[i,,drop=FALSE], main=mainlab, ylim=0:1)
            }
         }
      }
      dev.off()
   }
}
