CONFIG_FILE <- 'config.yaml'

# load project-specific parameters
config <- yaml::read_yaml(CONFIG_FILE)

# configure multi-threading
RhpcBLASctl::blas_set_num_threads(config$OPENBLAS_NUM_THREADS)
RhpcBLASctl::omp_set_num_threads(config$OMP_NUM_THREADS)
options(future.globals.maxSize = 8 * 1024 ^ 3)
options(parallelly.fork.enable = TRUE)
future::plan(strategy = 'future::multicore', workers = config$R_FUTURE_WORKERS)

# knitr and package options
knitr::opts_chunk$set(comment = NA, fig.width = 7, fig.height = 4, out.width = '70%',
                      warning = TRUE, error = TRUE, echo = TRUE, message = TRUE,
                      dpi = 300, rows.print=25)
options(ggrepel.max.overlaps = Inf)

# set seed, store time
set.seed(8569205)

SETUP_TIME <- proc.time()
