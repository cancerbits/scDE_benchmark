# Code repository for manuscript: Single-cell RNA-seq differential expression tests within a sample should use pseudo-bulk data of pseudo-replicates

Christoph Hafemeister and Florian Halbritter

St. Anna Children's Cancer Research Institute (CCRI), Vienna, Austria

**Abstract**

Single-cell RNA sequencing (scRNA-seq) has become a standard approach to investigate molecular differences between cell states. Comparisons of bioinformatics methods for the count matrix transformation (normalization) and differential expression (DE) analysis of these data have already highlighted recommendations for effective between-sample comparisons and visualization. Here, we examine two remaining open questions: (i) What are the best combinations of data transformations and statistical test methods, and (ii) how do pseudo-bulk approaches perform in single-sample designs? We evaluated the performance of 343 DE pipelines (combinations of eight types of count matrix transformations and ten statistical tests) on simulated and real-world data, in terms of precision, sensitivity, and false discovery rate. We confirm superior performance of pseudo-bulk approaches without prior transformation. For within-sample comparisons, we advise the use of three pseudo-replicates, and provide a simple R package [DElegate](https://github.com/cancerbits/DElegate) to facilitate application of this approach.

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

**Paper:** [bioRxiv preprint](https://doi.org/10.1101/2023.03.28.534443)

**R package:** [DElegate](https://github.com/cancerbits/DElegate)

**Data files:** Preprocessed data can be found on zenodo at [https://doi.org/10.5281/zenodo.7751830](https://doi.org/10.5281/zenodo.7751830).
