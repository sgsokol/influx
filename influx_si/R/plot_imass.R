plot_ti=function(ti, x, m=NULL, ...) {
   # plot time curse curves x[icurve, itime] and points from m
   # x and m are supposed to have the same dimension and organization
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
      if (length((nms <- strsplit(nm, ":", fixed=TRUE))[[1]]) > 2) {
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
#browser()
cat("***Warning: 'plot_imass.R' is deprecated. Use 'plot_ilab.R' instead.\n", file=fclog)
for (iexp in seq_len(nb_exp)) {
   usmf=jx_f$usmf[[iexp]] # unscaled simulated measurements
   me=measurements$vec$kin[[iexp]] # measured dynamic labeling data
   if (!is.null(measurements$outlier[[iexp]])) {
      iout=measurements$outlier[[iexp]]
      iout=iout[iout <= prod(dim(me))] # only outliers of label kinetics
      me[iout]=NA # measured dynamic labeling data
   }

   # get unique fragment names
   nm_sel=grep("^m:", if (is.null(rownames(me))) rownames(usmf) else rownames(me), v=TRUE)
   nm_sel=sort(nm_sel)
   pdf(sprintf("%s/%s.pdf", dirres, nm_exp[iexp]))
   if (length(nm_sel) > 0) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="MS measurements", cex=2)
      nmf=unique(apply(sapply(strsplit(nm_sel, ":", fixed=TRUE), "[", 1:4)[2:3,], 2, paste0, sep="", collapse=":"))
      for (metf in nmf) {
         i=grep(sprintf("m:%s:", metf), nm_sel, fix=TRUE, v=TRUE)
         #isim=pmatch(sapply(strsplit(i, ":", fixed=TRUE), function(v) paste0(v[-length(v)], collapse=":")), rownames(usmf))
         isim=grep(sprintf("m:%s:", metf), rownames(usmf), fix=TRUE, v=TRUE)
         plot_ti(tifull[[iexp]][-1L], usmf[isim,,drop=FALSE], me[i,,drop=FALSE], main=strsplit(metf, ":")[[1]][1], ylim=0:1)
      }
   }
   # plot non measured metabs from mid
   nm_sim=rownames(mid[[iexp]])
   nmm=if (length(nm_sel)) unique(sapply(strsplit(nm_sel, ":", fixed=TRUE), "[", 1:4)[2,]) else character(0L)
   nmmid=unique(sapply(strsplit(nm_sim, "+", fixed=TRUE), "[", 1L))
   if (emu) {
      nmmid=sapply(strsplit(nmmid, ":", fixed=TRUE), "[", 1L)
   }
   nmp=sort(setdiff(nmmid, nmm))
   if (length(nmp)) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="MS simulations", cex=2)
   }
   for (met in nmp) {
      if (emu) {
         i=grep(sprintf("^%s:", met), nm_sim, v=TRUE)
         # take fragments
         fr=unique(sapply(strsplit(i, "[+:]"), "[", 2L))
         for (f in fr) {
            i=grep(sprintf("^%s:%s\\+", met, f), nm_sim, v=TRUE)
            fi=as.integer(f)
            mainlab=if (fi == 2**clen[met]-1) met else sprintf("%s #%s", met, int2bit(fi, clen[met]))
            plot_ti(tifull[[iexp]][-1L], mid[[iexp]][i,,drop=FALSE], NULL, main=mainlab, ylim=0:1)
         }
      } else {
         i=grep(sprintf("^%s\\+", met), nm_sim, v=TRUE)
         plot_ti(tifull[[iexp]][-1L], mid[[iexp]][i,,drop=FALSE], NULL, main=met, ylim=0:1)
      }
   }
   dev.off()
}
