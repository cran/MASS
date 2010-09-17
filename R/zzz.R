.noGenerics <- TRUE

.onAttach <- function(...) require("stats")

.onUnload <- function(libpath)
    library.dynam.unload("MASS", libpath)

if(getRversion() < "2.6.0") nzchar <- function(x) nchar(x, "b") > 0
