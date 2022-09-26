# the idea is to source this setup script at the beginning of 
# all rmarkdown documents 
# it will read the project config file and set default options

# load project-specific parameters
config <- yaml::read_yaml(file = 'config.yaml')

# set knitr options
knitr::opts_chunk$set(comment = NA, fig.width = 7, fig.height = 4, out.width = '70%',
                      warning = TRUE, error = TRUE, echo = TRUE, message = TRUE,
                      dpi = 100)

# set some other package-specific options
options(ggrepel.max.overlaps = Inf)

options(future.plan = 'sequential', 
        future.globals.maxSize = 8 * 1024 ^ 3)

ggplot2::theme_set(ggplot2::theme_bw(base_size = 12))

# set a random seed
set.seed(993751)

# store the current CPU and real time
SETUP_TIME <- proc.time()

# source any scripts with commonly used functions
source('R/utilities.R')
