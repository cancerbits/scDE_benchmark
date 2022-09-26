# this top-level script may be used to run (knit) the individual
# Rmd files

# set up output path for the reports
config <- yaml::read_yaml("config.yaml")
report_dir <- file.path(config$out_root)

# example 1: knit one specific file
rmarkdown::render(input = 'Rmd/example.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# example 2: knit one specific file and pass parameters
rmarkdown::render(input = 'Rmd/example_with_parameters.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env(),
                  params = list(mean = 100, sd = 100),
                  output_file = 'example_with_parameters_100_100')
