## .................................................................................
## Purpose: Wrapper script for running fomes from command line
##
## Author: Nick Brazeau
##
## Date: 30 January, 2023
##
## Notes:
## .................................................................................

#............................................................
# dependencies
#...........................................................
deps <- c("fomes", "goodegg", "optparse", "readr", "magrittr", "stringr")
deps <- !sapply(deps, function(x){x %in% installed.packages()[,1]} ) # note this a named vector

# catch fomes remote
if(deps["fomes"]) {
  if (!"remotes" %in% installed.packages()[,1]){
    install.packages("remotes")
  }
  remotes::install_github("TriCEM/fomes")
  deps <- deps[names(deps) != "fomes"]
}

# catch goodegg remote
if(deps["goodegg"]) {
  if (!"remotes" %in% installed.packages()[,1]){
    install.packages("remotes")
  }
  remotes::install_github("nickbrazeau/goodegg")
  deps <- deps[names(deps) != "goodegg"]
}

# rest of deps
if (any(deps)) {
  install.packages(names(deps)[deps])
}

#......................
# call dependencies
#......................
library(fomes)
library(goodegg)
library(optparse)
library(readr)
library(magrittr)

#............................................................
# parse CL inputs
#...........................................................
option_list=list(
  make_option(c("-I", "--input"),
              type = "character", default = NULL,
              help = paste("Input path for YAML Map file controlling Gillespie SIR-NE Model Simulating",
                           "Assumes format of ' Iseed: 5 ' (i.e. param: input) for each model parameter in the `sim_Gillespie_SIR` function"),
              metavar = "character"),
  make_option(c("-O", "--output"),
              type = "character", default = NULL,
              help = paste("Output filename to write result of Gillespie SIR-NE Model Simulation",
                           "NB, the model output will be appended with .modelout.rds at the end of your filename path",
                           "and the final epidemic size will be appended with .finalsize.txt"
              ),
              metavar = "character")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if(is.null(opt$input)){
  print_help(opt_parser)
  stop("Missing input argument", call. = FALSE)
}

if(is.null(opt$output)){
  print_help(opt_parser)
  stop("Missing output argument", call. = FALSE)
}



#............................................................
# read YAML, extract, and make assertions
#TODO metaprogramming way to make var names from rows?
#...........................................................
modelYAML <- readr::read_delim(file = opt$input,
                               delim = ": ",
                               col_names = F) %>%
  magrittr::set_colnames(c("param", "val"))
#......................
# extract
#......................
Iseed <- as.numeric(modelYAML$val[modelYAML$param == "Iseed"])
N <- as.numeric(modelYAML$val[modelYAML$param == "N"])
beta <- as.numeric(unlist(stringr::str_split(modelYAML$val[modelYAML$param == "beta"], ",", simplify = F)))
dur_I <- as.numeric(modelYAML$val[modelYAML$param == "dur_I"])
init_contact_mat <- modelYAML$val[modelYAML$param == "init_contact_mat"]
rho <- as.numeric(modelYAML$val[modelYAML$param == "rho"])
initNC <- as.numeric(modelYAML$val[modelYAML$param == "initNC"])
term_time <- as.numeric(modelYAML$val[modelYAML$param == "term_time"])
return_contact_matrices <- as.logical(modelYAML$val[modelYAML$param == "return_contact_matrices"])
#......................
# catch
#......................
if (init_contact_mat != "NULL") {
  init_contact_mat <- readr::read_table(init_contact_mat, col_names = F)
  init_contact_mat <- as.matrix(init_contact_mat)
} else {
  init_contact_mat <- NULL # liftover
}


#............................................................
# run function and send out
#...........................................................
ret <- sim_Gillespie_SIR(Iseed = Iseed,
                         N = N,
                         beta = beta,
                         dur_I = dur_I,
                         init_contact_mat = init_contact_mat,
                         rho = rho,
                         initNC = initNC,
                         term_time = term_time,
                         return_contact_matrices = return_contact_matrices)

# save model
saveRDS(ret, paste0(opt$output, ".modelout.rds"))
# save final epidemic size
saveRDS(summary(ret)$FinalEpidemicSize, paste0(opt$output, "finalsize.txt"))
