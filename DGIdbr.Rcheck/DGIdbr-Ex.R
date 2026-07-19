pkgname <- "DGIdbr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('DGIdbr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("DGIdbr")
### * DGIdbr

flush(stderr()); flush(stdout())

### Name: DGIdbr
### Title: Main entry point for DGIdbr package
### Aliases: DGIdbr

### ** Examples

## Not run: 
##D DGIdbr(mode = "group", base_tables = ".", group_filename = "group.csv", base_out = ".")
## End(Not run)



cleanEx()
nameEx("drug_card")
### * drug_card

flush(stderr()); flush(stdout())

### Name: drug_card
### Title: Look up a drug: targets, mechanism, indications
### Aliases: drug_card

### ** Examples

## Not run: 
##D drug_card("ASPIRIN")
##D drug_card("CISPLATIN")
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
