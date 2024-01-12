# this top-level script may be used to run (knit) the individual
# Rmd files

# set up output path for the reports
config <- yaml::read_yaml("config.yaml")
report_dir <- file.path(config$out_root, 'reports')
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# set the names of the squair et al. data comparison files
fnames <- c('Hagai2018_mouse-lps.rds', 'Hagai2018_rat-lps.rds', 'Hagai2018_pig-lps.rds',
            'CanoGamez2020_Naive-iTreg.rds', 'CanoGamez2020_Naive-Memory_Resting.rds', 
            'CanoGamez2020_Th0-iTreg_Memory.rds', 'CanoGamez2020_Th0-Th17_Memory.rds', 
            'CanoGamez2020_Th0-Th2_Memory.rds', 'Angelidis2019_24m.rds',
            'Angelidis2019_3m.rds', 'Angelidis2019_alvmac.rds', 'Angelidis2019_pneumo.rds')
fnames <- c('Angelidis2019_24m.rds', 'Angelidis2019_3m.rds', 
            'CanoGamez2020_Naive-Memory_Resting.rds', 'CanoGamez2020_Naive-iTreg.rds', 
            'CanoGamez2020_Th0-Th17_Memory.rds', 'CanoGamez2020_Th0-iTreg_Memory.rds',
            'Hagai2018_mouse-lps.rds', 'Hagai2018_pig-lps.rds', 'Hagai2018_rat-lps.rds')


# download and process the single-cell immune data
# rmarkdown::render(input = 'Rmd/download_and_setup_immune_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# download and process the bulk immune data
# rmarkdown::render(input = 'Rmd/download_and_setup_blueprint_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# setup Squair et al. data
# rmarkdown::render(input = 'Rmd/setup_squair_et_al_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())
# rmarkdown::render(input = 'Rmd/squair_et_al_data_report.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# simulate data using muscat
# rmarkdown::render(input = 'Rmd/simulate_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

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
# rmarkdown::render(input = 'Rmd/run_pipelines_on_immune_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# run the pipelines on the squair et al data
# for (fname in fnames) {
#   ds_name <- gsub('\\.rds|\\.Rds', '', fname)
#   rmarkdown::render(input = 'Rmd/run_pipelines_on_squair_et_al_data.Rmd',
#                     output_dir = report_dir,
#                     knit_root_dir = config$project_root,
#                     envir = new.env(),
#                     params = list(fname = fname),
#                     output_file = sprintf('run_pipelines_on_squair_et_al_data_%s.html', ds_name))
# }

# summarize results on the immune data
# rmarkdown::render(input = 'Rmd/summarize_results_on_immune_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# summarize results on the squair et al data
# for (fname in fnames) {
#   ds_name <- gsub('\\.rds|\\.Rds', '', fname)
#   rmarkdown::render(input = 'Rmd/summarize_results_on_squair_et_al_data.Rmd',
#                     output_dir = report_dir,
#                     knit_root_dir = config$project_root,
#                     envir = new.env(),
#                     params = list(fname = fname),
#                     output_file = sprintf('summarize_results_on_squair_et_al_data_%s.html', ds_name))
# }

# rank pipelines on immune data
# rmarkdown::render(input = 'Rmd/pipeline_ranking_on_immune_data.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# rank pipelines on squair et al. data
# for (fname in fnames) {
#   ds_name <- gsub('\\.rds|\\.Rds', '', fname)
#   rmarkdown::render(input = 'Rmd/pipeline_ranking_on_squair_et_al_data.Rmd',
#                     output_dir = report_dir,
#                     knit_root_dir = config$project_root,
#                     envir = new.env(),
#                     params = list(fname = fname),
#                     output_file = sprintf('pipeline_ranking_on_squair_et_al_data_%s.html', ds_name))
# }

# combined ranking of pipelines
# rmarkdown::render(input = 'Rmd/pipeline_ranking_combined.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())

# create supplemental figures
# rmarkdown::render(input = 'Rmd/supp_figures.Rmd',
#                   output_dir = report_dir,
#                   knit_root_dir = config$project_root,
#                   envir = new.env())
