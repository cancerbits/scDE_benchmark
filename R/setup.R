CONFIG_FILE <- 'config.yaml'

# load project-specific parameters
config <- yaml::read_yaml(CONFIG_FILE)

# options regarding parallelization
#Sys.setenv(OMP_NUM_THREADS = 1)
#Sys.setenv(OMP_THREAD_LIMIT = 1)
#Sys.setenv(OPENBLAS_NUM_THREADS = 1)
#Sys.setenv(MKL_NUM_THREADS = 1)

RhpcBLASctl::blas_set_num_threads(config$OPENBLAS_NUM_THREADS)
options(future.globals.maxSize = 8 * 1024 ^ 3)
options(parallelly.fork.enable = TRUE)
future::plan(strategy = 'future::multicore', workers = config$R_FUTURE_WORKERS)

# knitr and package options
knitr::opts_chunk$set(comment = NA, fig.width = 7, fig.height = 4, out.width = '70%',
                      warning = TRUE, error = TRUE, echo = TRUE, message = TRUE,
                      dpi = 300, rows.print=25)
options(ggrepel.max.overlaps = Inf)

# set seed, store time, source other scripts with function definitions
set.seed(8569205)

SETUP_TIME <- proc.time()
