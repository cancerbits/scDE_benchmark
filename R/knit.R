# this top-level script may be used to run (knit) the individual
# Rmd files

# set up output path for the reports
config <- yaml::read_yaml("config.yaml")
report_dir <- file.path(config$out_root, 'reports')
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# download and process the single-cell immune data
rmarkdown::render(input = 'Rmd/download_and_setup_immune_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# download and process the bulk immune data
rmarkdown::render(input = 'Rmd/download_and_setup_blueprint_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# simulate data using muscat
rmarkdown::render(input = 'Rmd/simulate_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# run the pipelines on the simulated data
rmarkdown::render(input = 'Rmd/run_pipelines_on_simulated_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# summarize results on the simulated data
rmarkdown::render(input = 'Rmd/summarize_results_on_simulated_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# rank pipelines on simulated data
rmarkdown::render(input = 'Rmd/pipeline_ranking_on_simulated_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# run the pipelines on the immune data
rmarkdown::render(input = 'Rmd/run_pipelines_on_immune_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# summarize results on the immune data
rmarkdown::render(input = 'Rmd/summarize_results_on_immune_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())


# rank pipelines on immune data
rmarkdown::render(input = 'Rmd/pipeline_ranking_on_immune_data.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# combined ranking of pipelines
rmarkdown::render(input = 'Rmd/pipeline_ranking_combined.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())

# create supplemental figures
rmarkdown::render(input = 'Rmd/supp_figures.Rmd',
                  output_dir = report_dir,
                  knit_root_dir = config$project_root,
                  envir = new.env())
