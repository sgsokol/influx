source("aide_plot.R")
# read simulated and real measurements
f="chlamy_i_comp_10_res.kvh"; td=get_usm(f); nmm=rownames(td$usm); ti=td$ti; me=get_labcin("chlamy_pulse_comp_final.txt", nmm)
# plot similated vs measured data for a given metabolite
i=grep(":Glc", nmm, fix=T); plot_ti(ti, td$usm[i,,drop=F], me[i,,drop=F])

# read simulated cumomers
xs=kvh_get_matrix(f, "simulated cumomers"); tis=as.numeric(colnames(xs)); nmx=rownames(xs); ms=cumo2mass(xs); nmma=rownames(ms)
# plot cumomers of weight 1 (their index is an integer power of 2)
ix=grep("^Sed.*_cy", nmma, fix=F); plot_ti(tis, ms[ix,])

# write one metabolite per file pdf
mets=unique(sapply(nmm, function(n) strsplit(n, "[:+]")[[1]][2]))
sapply(mets, function(m) {
   cat(m, "\n");
   i=grep(join("", c(":", m)), nmm, fix=TRUE);
   if (length(i) > 0) {
      pdf(join("", c(m, ".pdf")));
      plot_ti(td$ti, td$usm[i,,drop=FALSE], me[i,,drop=FALSE], main=m);
      dev.off()
   } else {
      cat("not found: ", m, "\n")
   }
})
