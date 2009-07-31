if(!nzchar(Sys.getenv("MASS_TESTING"))) q("no")
unlink("scripts", recursive = TRUE)
dir.create("scripts")
Sys.unsetenv("R_TESTS")
setwd("scripts")

runone <- function(f)
{
    message("  Running ", sQuote(basename(f)))
    outfile <- paste(basename(f), "out", sep = "")
    failfile <- paste(outfile, "fail", sep=".")
    unlink(c(outfile, failfile))
    cmd <- paste(shQuote(file.path(R.home(), "bin", "R")),
                 "CMD BATCH --no-save",
                 shQuote(f), shQuote(outfile))
    res <- system(cmd)
    if (res) {
        cat(tail(readLines(outfile), 20), sep="\n")
        file.rename(outfile, failfile)
        return(1L)
    }
    0L
}


library(MASS)
dd <- system.file("scripts", package="MASS")
files <- list.files(dd, pattern="\\.R$", full.names=TRUE)
res <- 0L
for(f in files) res <- res + runone(f)

proc.time()

if(res) stop(gettextf("%d scripts failed", res))
