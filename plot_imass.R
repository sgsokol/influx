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
      # make names look like Metab:frag+0 for m0 Metab:frag+1 for m1 etc.
      nm=matrix(unlist(strsplit(nm, ":")), byrow=T, ncol=4)[,c(2:4), drop=F]
      nm=paste(paste(nm[,1L], nm[,2L], sep=":"), nm[,3L], sep="+")
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
      matplot(ti, t(x), t="l", lty=seq_len(nb_curve), col=seq_len(nb_curve), lwd=2, add=T, ...)
   if (!is.null(m)) {
      # plot measured points
      for (i in seq_len(nrow(m))) {
         inna=which(!is.na(m[i,]))
         if (length(inna) == 0)
            next
         points(tim[inna], m[i,inna], col=i, cex=0.5, t="b", lwd=0.5, ...)
         if (nrow(m) == nb_curve && !is.null(x)) {
            # draw filled polygons between simulation and data
            polygon(c(ti,rev(tim[inna])), c(x[i,], rev(m[i,inna])), col=rgb(red=t(col2rgb(i)), alpha=31, max=255), bord=F, ...)
         }
      }
   }
   legend("topright", legend=nm, lty=1:nb_curve, col=1:nb_curve, lwd=2, cex=0.75, bg=rgb(0.97,0.97,0.97, 0.75))
}

for (iexp in seq_len(nb_exp)) {
   usmf=jx_f$usmf[[iexp]] # unscaled simulated measurements
   me=measurements$vec$kin[[iexp]] # measured dynamic labeling data
   if (!is.null(measurements$outlier[[iexp]])) {
      iout=measurements$outlier[[iexp]]
      iout=iout[iout <= prod(dim(me))] # only outliers of label kinetics
      me[iout]=NA # measured dynamic labeling data
   }

   # get unique metab names
   nm_sel=grep("^m:", if (is.null(rownames(me))) rownames(usmf) else rownames(me), v=T)
   nm_sel=sort(nm_sel)
   if (length(nm_sel) > 0) {
      nm=unique(sapply(strsplit(nm_sel, ":"), "[", 1:4)[2,])
      pdf(sprintf("%s/%s.pdf", dirw, nm_exp[iexp]))
      for (met in nm) {
         i=grep(sprintf("m:%s:", met), nm_sel, fix=T, v=T)
         isim=pmatch(sapply(strsplit(i, ":"), function(v) paste0(v[-length(v)], collapse=":")), rownames(usmf))
         plot_ti(tifull[[iexp]][-1L], usmf[isim,,drop=FALSE], me[i,,drop=F], main=met, ylim=0:1)
      }
      dev.off()
   }
}
