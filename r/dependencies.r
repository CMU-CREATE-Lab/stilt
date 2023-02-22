# STILT Dependency Loader
# For documentation, see https://github.com/uataq/stilt
# Ben Fasoli

if (!'stilt_wd' %in% ls())
  stilt_wd <- getwd()

perf0 <- 0
perf1 <- 0
perf2 <- 0
perf3 <- 0
perf4 <- 0

perf0 <- perf0 - as.numeric(Sys.time())
# Source r/src R scripts
rsc <- dir(file.path(stilt_wd, 'r', 'src'), pattern = '.*\\.r$', full.names = T)
invisible(lapply(rsc, source))
perf0 <- perf0 + as.numeric(Sys.time())

perf1 <- perf1 - as.numeric(Sys.time())
# Load external libraries
if (!'lib.loc' %in% ls()) lib.loc <- NULL
libs <- load_libs('dplyr',
                  'ncdf4',
                  'parallel',
                  'raster',
                  'rslurm',
                  lib.loc = lib.loc)
perf1 <- perf1 + as.numeric(Sys.time())

perf2 <- perf2 - as.numeric(Sys.time())
# Load permute fortran dll for footprint matrix permutation
permute_exe <- file.path(stilt_wd, 'r/src/permute.so')
if (!file.exists(permute_exe))
  stop('calc_footprint(): failed to find permute.so in r/src/')
dyn.load(permute_exe)
perf2 <- perf2 + as.numeric(Sys.time())

perf3 <- perf3 - as.numeric(Sys.time())
# Validate arguments and load dependencies if necessary
if ((!class(projection) == 'function') && ('projection' %in% ls()))
  validate_projection(projection)
if ('varsiwant' %in% ls())
  validate_varsiwant(varsiwant)
if (all(c('xmn', 'xmx', 'ymn', 'ymx') %in% ls()))
  validate_footprint_extent(xmn, xmx, ymn, ymx)
perf3 <- perf3 + as.numeric(Sys.time())

# Disable grouping message from dplyr >=1.0.0
options(dplyr.summarise.inform = F)

message("<<<    dependencies.r, perf0: ", perf0)
message("<<<    dependencies.r, perf1: ", perf1)
message("<<<    dependencies.r, perf2: ", perf2)
message("<<<    dependencies.r, perf3: ", perf3)