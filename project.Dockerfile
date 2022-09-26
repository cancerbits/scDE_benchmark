# pull base image
FROM rocker/tidyverse:4.1.3

# who maintains this image
LABEL maintainer Humble Hornet "humble.hornet@ccri.at"
LABEL version 4.1.3-v1

# change some permissions
RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

# The default CRAN mirror is set to a snapshot from 2022-04-21
# Bioconductor is at version 3.14

# install a couple of R packages that we commonly use
#RUN apt-get -y update && apt-get -y install \
#  libglpk40 \
#  libxt6 \
#  && apt-get clean \
#  && rm -rf /tmp/* /var/tmp/*
#RUN R -e "BiocManager::install(c('markdown', 'Seurat', 'SeuratObject', 'sparseMatrixStats', 'edgeR', 'apeglm', 'DESeq2', 'fgsea', 'hypeR', 'patchwork'))"
#RUN R -e "remotes::install_github(repo = 'cancerbits/canceRbits', ref = '1f3dd76')"

# more package installation below this line
