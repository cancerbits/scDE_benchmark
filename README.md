# Code repository for: Paper

Christoph Hafemeister and Florian Halbritter

1: St. Anna Children's Cancer Research Institute (CCRI), Vienna, Austria

**Abstract**

Abstract

## Repository structure

* `project.Dockerfile` defines the environment used to carry out all experiments
* `config.yaml` is used to set paths and multi-threading parameters
* `R/` holds R function definitions and misc utility scripts
* `Rmd/` holds R markdown documents for the individual steps of the project
* `bash/` holds shell scripts to build and run the docker image, and to parse the config file
* `pipelines/` holds two csv files defining the names and parameterizations of transformations and DE methods

## Reproducing the results

The file `R/knit.R` calls all `Rmd/*.Rmd` files in order to reproduce the analysis.

## Links

**Paper:** [https://doi.org](https://doi.org)

**Data files:** Preprocessed data can be found on zenodo at [https://doi.org/10.5281/zenodo.7751830](https://doi.org/10.5281/zenodo.7751830).
