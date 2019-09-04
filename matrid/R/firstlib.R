.First.lib <- function(lib,pkg)
{
   library.dynam("matrid",pkg,lib)
   cat("matrid 1.0 loaded\n")
}
