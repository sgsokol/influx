# CRAN packages required for influx_si
preq=c("nnls", "rmumps", "arrApply", "slam", "limSolve")
# already installed packages
pinst=rownames(installed.packages())

# install lacking packages
plack=setdiff(preq, pinst)
if (length(plack))
   install.packages(plack, repos="https://cloud.r-project.org/")
# update already present
pold=setdiff(preq, plack)
if (length(pold))
   update.packages(oldPkgs=pold, repos="https://cloud.r-project.org/")

# install multbxxc from github
if (.Platform$OS.type == "windows") {
   install.packages("https://github.com/sgsokol/multbxxc/raw/master/multbxxc_1.0.zip", repos=NULL)
} else if (Sys.info()[["sysname"]] == "Darwin") {
   install.packages("https://github.com/sgsokol/multbxxc/raw/master/multbxxc_1.0.tgz", repos=NULL)
} else {
   install.packages("https://github.com/sgsokol/multbxxc/raw/master/multbxxc_1.0.tar.gz", repos=NULL)
}
