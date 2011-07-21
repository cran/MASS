.noGenerics <- TRUE

.onUnload <- function(libpath)
    library.dynam.unload("MASS", libpath)

if(getRversion() < "2.14") {
    ## dummies
    dev.hold <- dev.flush <- function(level) {}
}
