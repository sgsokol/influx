# CRAN packages required for influx_si
options(error=function(e) q("no"))
preq=c("nnls", "rmumps", "arrApply", "slam", "limSolve", "multbxxc")
# check what is already installed
iinst=sapply(preq, requireNamespace, quietly = TRUE)
repos=getOption("repos")
repos=if (!length(repos) || repos == "@CRAN@") "https://cloud.r-project.org/" else repos
# install lacking packages
plack=preq[!iinst]
if (length(plack)) {
   #print(interactive())
   #cat("installing: ", plack, "\n")
   install.packages(plack, repos=repos)
}
# update already present
pold=preq[iinst]
if (length(pold)) {
   #cat("updating: ", pold, "\n")
   update.packages(oldPkgs=pold, repos=repos, ask=FALSE)
}
q("no")
