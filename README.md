# project_template
Template repository for genomics data analysis projects

# Usage

## Option 1
Use the "Use this template" button on github.com

## Option 2
Clone this repo, remove the .git folder, go from there

# Customizing the template

## config.yaml

Use this file to save project-specific parameters, e.g. set paths. This way all parameters may be stored in one place and can be loaded by scripts easily.

## Docker

This is assuming that you are using a project-specific docker image (you don't have to). If you want to use the `run_docker*` scripts, you will have to edit the `project_docker` line in `config.yaml` to point to the image for the project.

- Edit `project.Dockerfile` 
- Build image: `docker build -t cancerbits/dockr:project_name - < project.Dockerfile` (replace `project_name`)

`run_docker_rstudio.sh` expects two parameters: the port number and a password. E.g. `./run_docker_rstudio.sh 8881 abcd` will create a docker container where RStudio Server is mapped to port 8881 and the login password is abcd.
