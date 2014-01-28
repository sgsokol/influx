# 2014-01-08 sokol
# translate expa format to FTBL/EQUALITY/NET section + FLUX/{NET,XCH}
# May be, the result file needs to be manually edited, e.g. some
# reactions must be removed or some reactions must be replaced by numerical values

# NB. In metabolite and reaction names, "-" are replaced by "_"

# usage: R --vanilla --slave --args file.expa < expa2ftbl.R > file.ftbl_eq

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
   cat("usage: R --vanilla --slave --args file.expa < expa2ftbl.R > file.ftbl_eq_net\n", file=stderr())
   stop("No expa file name provided.")
}
if (length(args) > 1) {
   cat("usage: R --vanilla --slave --args file.expa < expa2ftbl.R > file.ftbl_eq_net\n", file=stderr())
   stop("Only one expa file is expected.")
}

#cat(args, "\n", file=stderr())
stoech=list() # indexed by reac name, item=vector of named stoechiom coeffs
all_metabs=c()
start_proc=F
cnx=file(args, "rb")
while(length(line<-readLines(cnx, 1L)) > 0L) {
   #cat("got: ", line, "\n")
   if (line=="(Internal Fluxes)") {
      start_proc=T
      next
   }
   if (line=="(Exchange Fluxes)") {
      break
   }
   if (!start_proc) {
      next
   }
   # process (Internal Fluxes) here
   items=strsplit(line, "\t")[[1]]
   #print(items)
   nit=length(items)
   if (nit==0) {
      next
   }
   reac=gsub("-", "_", items[1])
   coefs=as.numeric(items[seq(3, by=2, len=(nit-2)/2)])
   mets=gsub("-", "_", items[seq(4, by=2, len=(nit-2)/2)])
   names(coefs)=mets
   all_metabs=c(all_metabs, mets)
   stoech[[reac]]=coefs
}
all_metabs=sort(unique(all_metabs))
all_reacs=sort(names(stoech))
nmet=length(all_metabs)
nreac=length(stoech)

# make stoechiometric matrix by making balance of each metabolite
stmat=matrix(0., nrow=nmet, ncol=nreac)
dimnames(stmat)=list(all_metabs, all_reacs)
for (ireac in seq(along=all_reacs)) {
   reac=all_reacs[[ireac]]
   coefs=stoech[[reac]]
   stmat[names(coefs),ireac]=coefs
}

# get carbonated metabolites from file.ftbl (.expa is removed and replaced by .ftbl)
fbase=substring(args, 1L, nchar(args)-5L)
system(sprintf("ftbl2netan.py %s > %s.netan", fbase, fbase))
metcarb=system(sprintf("grep -e ^metabs %s.netan", fbase), intern=T)
reaccarb=system(sprintf("grep -e ^reac %s.netan", fbase), intern=T)
# strip unsueful chars at the beginnig and the end
if (length(metcarb) > 0L && nchar(metcarb) > 16L) {
   metcarb=strsplit(substring(metcarb, 14, nchar(metcarb)-3), "', '")[[1]]
} else {
   metcarb=""
}
# strip unsueful chars at the beginnig and the end
if (length(reaccarb) > 0L && nchar(reaccarb) > 14L) {
   reaccarb=strsplit(substring(reaccarb, 12, nchar(reaccarb)-3), "', '")[[1]]
} else {
   reaccarb=""
}

# strip out some metabolites
# - already present in the NETWORK section
# - input
# - output
# - leading to linearly dependent equations

# remove NETWORK metabs
meteq=setdiff(all_metabs, metcarb)
stmat=stmat[meteq,,drop=F]

# remove input and output metabs
inp_metabs=which(apply(stmat, 1L, function(row) sum(row>0)==0L))
inp_flux=unlist(apply(stmat[inp_metabs,,drop=F], 1L, function(row) names(which(row<0))))
inp_flux=setdiff(inp_flux, reaccarb)
out_metabs=which(apply(stmat, 1L, function(row) sum(row<0)==0L))
out_flux=unlist(apply(stmat[out_metabs,,drop=F], 1L, function(row) names(which(row>0))))
out_flux=setdiff(out_flux, reaccarb)

noncarb_dep=colnames(stmat)
noncarb_dep=setdiff(noncarb_dep, reaccarb)
noncarb_dep=setdiff(noncarb_dep, inp_flux)
noncarb_dep=setdiff(noncarb_dep, out_flux)

if (length(inp_metabs) > 0 || length(out_metabs) > 0) {
   stmat=stmat[-c(inp_metabs, out_metabs),,drop=F]
}

qs=qr(t(stmat), LAPACK=T)
d=abs(diag(qs$qr))
rank=sum(d > d[1]*1.e-10)
if (rank == 0L) {
   stop(sprintf("No linearly independent equalities in the '%s' file.", args))
}
stmat=stmat[qs$pivot[1:rank],]

all_metabs=rownames(stmat)

# write the EQUALITY section on stdout
# header
cat("EQUALITIES\n\tNET\n\t\tVALUE\tFORMULA\n")

for (met in all_metabs) {
   row=stmat[met,]
   iposi=which(row > 0)
   rposi=names(iposi)
   cposi=sapply(names(row[iposi]), function(m) {val=row[m]; if (val != 1) return(paste(val, "*", names(val), sep="")) else return(m)})
   # metabolite name in comment
   cat("\t\t// ", met, "\n", sep="")
   cat("\t\t0\t")
   if (length(cposi)) {
      cat(paste(cposi, collapse="+"))
   }
   inega=which(row < 0)
   rnega=names(inega)
   cnega=sapply(names(row[inega]), function(m) {val=row[m]; if (val != -1) return(paste(-val, "*", names(val), sep="")) else return(m)})
   if (length(cnega)) {
      cat("-", paste(cnega, collapse="-"), sep="")
   }
   cat("\n")
}

# write the FLUX section on stdout
# header
cat("FLUX\n\tNET\n\t\tNAME\tFCD\tVALUE(F/C)\tED_WEIGHT\tLOW(F)\tINC(F)\tUP(F)\n")
# NET subsection
for (fl in unique(c(inp_flux, out_flux, noncarb_dep)) {
   cat("\t\t", fl, "\tD\n", sep="")
}
# XCH subsection
cat("\tXCH\n\t\tNAME\tFCD\tVALUE(F/C)\tED_WEIGHT\tLOW(F)\tINC(F)\tUP(F)\n")
for (fl in unique(c(inp_flux, out_flux, noncarb_dep))) {
   cat("\t\t", fl, "\tC\t0\n", sep="")
}
