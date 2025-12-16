#' store sparse matrices in <pref>_01.mtx, <pref>_02.mtx, etc
#'
#' Files <pref>_[0-9]+.mtx are silently cleaned up before writing if exist.
#' Files are written in MTX format, i.e. each row is a triplicate (i, j, x)
#' The first row is a triplicate (nrow, ncol, nb_nonzero)
#' The number in file name is the cumomer weight.
mat2file=function(param, labargs, pref=paste0("cumo_", baseshort, "_"), dres=file.path(dirres, "tmp")) {
    pbase=file.path(dres, pref)
    fli=list.files(dres, pattern=paste0(pref, "[0-9]+\\.mtx"))
    suppressWarnings(file.remove(fli))
    ndig=floor(log10(length(spa)))+1L
    dfmt=sprintf("%%0%dd", ndig)
    lf=param2fl(param, labargs)
    for (w in seq_along(labargs$spa)) {
        fname=file.path(paste0(pbase, sprintf(dfmt, w), ".mtx"))
        atri=fwrv2Abr(lf$fwrv, labargs$spa[[w]], NULL, getA=TRUE, getb=FALSE, labargs$emu)$A$triplet()
        mat2mtx(atri, fname)
    }
}
#browser()
mat2file(param, labargs)
