plot_ms=function(x, m=NULL, dev=NULL, ...) {
   # plot stacked ms bars from x and m
   # x and m are supposed to have the same dimension and organization
   n=length(x)
   nm=names(x)
   
   if (length(m)) {
      m=matrix(m, nrow=n)
      nc=ncol(m)
      colnames(m)=rep("meas", nc)
   } else {
      nc=1
   }
   if (length((nms <- strsplit(nm, ":", fixed=TRUE))[[1]]) > 2) {
      nm_leg=as.expression(sapply(nms, function(v) substitute(M[i], list(i=v[4]))))
   } else {
      nm_leg=as.expression(sapply(seq_len(n)-1, function(v) substitute(M[i], list(i=v))))
   }
   xlim=range(0, if (length(m) && !is.null(dev)) min(colSums(m, na.rm=TRUE)-dev[1]) else 0, sum(x), if (length(m) && !is.null(dev)) max(col_sums(m, na.rm=TRUE))+dev else 0)*1.1
   bc=barplot(cbind(sim=x, m), width=0.1*diff(xlim), horiz=TRUE, xlab="MS fraction", ylim=c(0, 0.25*nc), xlim=xlim, asp=TRUE, col=co,
      legend.text=nm_leg, args.legend=list(cex=0.75, horiz=TRUE), ...)
   if (!is.null(m) && !is.null(dev)) {
      for (j in seq(ncol(m))) {
         cume=cumsum(m[,j])
         h=seq(-1, 1, len=n)*0.03
         segments(cume-dev, bc[j+1]+h, cume+dev, bc[j+1]+h)
         arrows(cume-dev, bc[j+1]+h, cume+dev, bc[j+1]+h, length=0.05, angle=80, code=3)
      }
   }
}
plot_lab=function(x, m=NULL, dev=NULL, ...) {
   # plot stacked bars of LABEL_MEASUREMENTS
   # x and m are supposed to have the same dimension and organization
   n=length(x)
   nm=names(x)
   if (length((nms <- strsplit(nm, ":", fixed=TRUE))[[1]]) > 1) {
      nm_leg=sapply(nms, "[", 3)
   } else {
      nm_leg=paste("#", seq_len(n), sep="")
   }
   xlim=range(0, sum(x), sum(m)+dev)*1.1
   bc=barplot(cbind(sim=x, meas=m), width=0.08*diff(xlim), horiz=TRUE, xlab="Label fraction", ylim=xlim*0.25, xlim=xlim, asp=TRUE, col=co,
      legend.text=nm_leg, args.legend=list(cex=0.75, horiz=TRUE), ...)
   if (!is.null(m) && !is.null(dev)) {
      cume=cumsum(m)
      h=seq(-1, 1, len=n)*0.03
      segments(cume-dev, bc[2]+h, cume+dev, bc[2]+h)
      arrows(cume-dev, bc[2]+h, cume+dev, bc[2]+h, length=0.05, angle=80, code=3)
   }
}
plot_peak=function(x, m=NULL, dev=NULL, ...) {
   # plot stacked bars of PEAK_MEASUREMENTS
   # x and m are supposed to have the same dimension and organization
   n=length(x)
   nm=names(x)
   if (length((nms <- strsplit(nm, ":", fixed=TRUE))[[1]]) > 1) {
      nm_leg=sapply(nms, "[", 4)
   } else {
      nm_leg=seq_len(n)
   }
   bc=barplot(cbind(sim=x, meas=m), width=0.1, horiz=TRUE, xlab="Peak fraction", ylim=c(0, 0.25), xlim=range(0, sum(x), sum(m)+dev, if (length(m)) m[1]-dev else 0)*1.1, asp=TRUE, col=co,
      legend.text=nm_leg, args.legend=list(cex=0.75, horiz=TRUE), ...)
   if (!is.null(m) && !is.null(dev)) {
      cume=cumsum(m)
      h=seq(-1, 1, len=n)*0.03
      segments(cume-dev, bc[2]+h, cume+dev, bc[2]+h)
      arrows(cume-dev, bc[2]+h, cume+dev, bc[2]+h, length=0.05, angle=80, code=3)
   }
}
plot_flux=function(x, m=NULL, dev=NULL, ...) {
   # plot FLUX_MEASUREMENTS
   # x and m are supposed to have the same dimension and organization
   n=length(x)
   xlim=range(0, x, m-dev, m+dev)*1.1
   bc=barplot(cbind(sim=x, meas=m), width=0.05*diff(xlim), horiz=TRUE, xlab="Flux value", ylim=c(-0.1, 0.2)*diff(xlim), xlim=xlim, col=co, ...)
   if (!is.null(m) && !is.null(dev)) {
      cume=cumsum(m)
      segments(cume-dev, bc[2], cume+dev, bc[2])
      arrows(cume-dev, bc[2], cume+dev, bc[2], length=0.1, angle=80, code=3)
   }
}
plot_pool=function(x, m=NULL, dev=NULL, ...) {
   # plot METAB_MEASUREMENTS
   # x and m are supposed to have the same dimension and organization
   n=length(x)
   xlim=range(0, x, m-dev, m+dev)*1.1
   w=0.1*diff(xlim)
   bc=barplot(cbind(sim=x, meas=m), width=w, horiz=TRUE, xlab="Metabolite concentration", ylim=c(-0.08, 0.17)*20*w, xlim=xlim, col=co, ...)
   if (!is.null(m) && !is.null(dev)) {
      cume=cumsum(m)
      segments(cume-dev, bc[2], cume+dev, bc[2])
      arrows(cume-dev, bc[2], cume+dev, bc[2], length=0.1, angle=80, code=3)
   }
}
# colors from http://tools.medialab.sciences-po.fr/iwanthue/ (param: 20, soft (k-means)))
co=c(
"#45aecf",
"#d0502b",
"#7f63d7",
"#66b344",
"#c34db6",
"#b9b53d",
"#4e71c0",
"#da9234",
"#8e9bdf",
"#8c8c42",
"#cb86d4",
"#437e43",
"#d03f7c",
"#59bf8f",
"#d0454e",
"#8a5191",
"#965e2a",
"#de82a2",
"#df926c",
"#a84f5b"
)
if (is.null(.GlobalEnv$jx_f) || is.null(jx_f$simlab)) {
   stop_mes("plot_smeas.R: simulated data are not available. Plotting skipped", file=fcerr)
}
for (iexp in seq_len(nb_exp)) {
   sim=jx_f$simlab[[iexp]]
   me=measurements$vec$labeled[[iexp]] # measured stationary ms data

   # get unique met-fragment names
   nm_sel=grep("^m:", if (is.null(names(me))) names(sim) else names(me), v=TRUE)
   pdf(sprintf("%s/%s.pdf", dirw, nm_exp[iexp]), width=8, height=6)
   if (length(nm_sel) > 0) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="MS measurements\n(error bars=±2*dev)", cex=2)
      nmf=sort(unique(apply(sapply(strsplit(nm_sel, ":", fixed=TRUE), "[", 1:4)[2:3,], 2, paste0, sep="", collapse=":")))
      for (metf in nmf) {
         i=grep(sprintf("m:%s:", metf), nm_sel, fixed=TRUE, v=TRUE)
         # count repeated fragments
         nbf=length(grep(sprintf("m:%s:0", metf), i, fixed=TRUE, v=TRUE))
         if (nbf > 1) {
            isim=i[seq(length(i)/nbf)]
         } else {
            isim=i
         }
         mf=strsplit(metf, ":")[[1]]
         met=mf[1]
         fr=mf[2]
         metlen=clen[mets_in_res[[iexp]][i[1]]]
         if (fr == paste(seq_len(metlen), collapse=",") || fr == sprintf("1~%d", metlen)) {
            mainlab=met
         } else {
            mask=rep("0", metlen)
            mask[eval(parse(text=paste0("c(", sub("~", ":", fr, fixed=TRUE), ")")))]="1"
            mainlab=sprintf("%s #%s", met, paste0(mask, collapse=""))
         }
         plot_ms(sim[isim], me[i], 2*measurements$dev$labeled[[iexp]][i], main=mainlab)
      }
   }
   # plot MS of non measured metabs
#browser()
   if (!is.null(.GlobalEnv$mid)) {
      nm_sim=names(mid[,iexp])
      # get "pyr" from "m:pyr:1~3:0:171" (on ms)
      nmm=unique(sapply(strsplit(nm_sel, ":", fixed=TRUE), "[", 1:4)[2,])
      # get "pyr" from "pyr:7+0" (on simulated)
      nmmid=unique(sapply(strsplit(nm_sim, "+", fixed=TRUE), "[", 1))
      if (emu)
         nmmid=unique(sapply(strsplit(nmmid, ":", fixed=TRUE), "[", 1))
      nmp=sort(setdiff(nmmid, nmm))
      if (length(nmp)) {
         plot(0:1, c(0,0.1), type="n", axes=FALSE, xlab="", ylab="")
         text(0.5, 0.05, lab="MS simulations", cex=2)
      }
      for (met in nmp) {
         if (emu) {
            i=grep(sprintf("^%s:", met), nm_sim, v=TRUE)
            # take fragments
            fr=unique(sapply(strsplit(i, "[+:]"), "[", 2))
            for (f in fr) {
               i=grep(sprintf("^%s:%s\\+", met, f), nm_sim, v=TRUE)
               fi=as.integer(f)
               mainlab=if (fi == 2**clen[met]-1) met else sprintf("%s #%s", met, int2bit(fi, clen[met]))
               plot_ms(mid[i, iexp], NULL, NULL, main=mainlab)
            }
         } else {
            i=grep(sprintf("^%s\\+[0-9]+$", met), nm_sim, v=TRUE)
            plot_ms(mid[i, iexp], NULL, NULL, main=met)
         }
      }
   }
   # LABEL_MEASURMENTS
   nm_sel=grep("^l:", if (is.null(names(me))) names(sim) else names(me), v=TRUE)
   if (length(nm_sel) > 0) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="Label measurements\n(error bars=±2*dev)", cex=2)
      nm=sort(unique(sapply(strsplit(nm_sel, ":", fixed=TRUE), "[", 2)))
      for (met in nm) {
         i=grep(sprintf("l:%s:", met), nm_sel, fixed=TRUE, v=TRUE)
         plot_lab(sim[i], me[i], 2*measurements$dev$labeled[[iexp]][i], main=met)
      }
   }
   # PEAK_MEASURMENTS
   nm_sel=grep("^p:", if (is.null(names(me))) names(sim) else names(me), v=TRUE)
   if (length(nm_sel) > 0) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="Peak measurements\n(error bars=±2*dev)", cex=2)
      nmp=sort(unique(apply(sapply(strsplit(nm_sel, ":", fixed=TRUE), "[", 2:3), 2, paste, collapse=":")))
      for (metp in nmp) {
         i=grep(sprintf("p:%s:", metp), nm_sel, fix=TRUE, v=TRUE)
         plot_peak(sim[i], me[i], 2*measurements$dev$labeled[[iexp]][i], main=metp)
      }
   }
   # FLUX_MEASURMENTS
   nm_sel=names(fmn)
   if (length(nm_sel) > 0) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="Flux measurements\n(error bars=±2*dev)", cex=2)
      nmf=sort(substring(nm_sel, 5))
      for (f in nmf) {
         i=grep(sprintf("^.\\.n\\.%s$", f), nm_sel, v=TRUE)
         plot_flux(jx_f$simfmn[i], fmn[i], 2*measurements$dev$flux[i], main=f)
      }
   }
   # POOL_MEASURMENTS
   nm_sel=names(measurements$vec$pool)
   if (length(nm_sel) > 0) {
      plot(0:1, c(0, 0.1), type="n", axes=FALSE, xlab="", ylab="")
      text(0.5, 0.05, lab="Metabolite pool measurements\n(error bars=±2*dev)", cex=2)
      nmp=sort(substring(nm_sel, 4))
      for (p in nmp) {
         i=grep(sprintf("pm:%s", p), nm_sel, fixed=TRUE, v=TRUE)
         plot_pool(jx_f$simpool[i], measurements$vec$pool[i], 2*measurements$dev$pool[i], main=p)
      }
   }
  dev.off()
}
