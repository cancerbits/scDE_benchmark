# pull base image (rocker/tidyverse 4.3.1 from June 27, 2023)
FROM rocker/tidyverse@sha256:35d47a2d3bacf0ada4fe5498e8b6870e8c0d883b7978c1dcad317f91cbfb7d0d

# who maintains this image
LABEL maintainer Christoph Hafemeister "christoph.hafemeister@ccri.at"
LABEL version 4.3.1-v1

# change some permissions
RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

# Change default CRAN mirror to a particular snapshot in time
ENV CRAN=https://packagemanager.posit.co/cran/__linux__/jammy/2023-08-14
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-08-14'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site
# To install everything from source, use the repo below
#ENV CRAN=https://packagemanager.rstudio.com/cran/2023-08-14
#RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2023-08-14'), download.file.method = 'libcurl')" >> ${R_HOME}/etc/Rprofile.site

# update packages
RUN R -e "update.packages(ask = FALSE)"
RUN R -e "BiocManager::install(ask = FALSE, update = TRUE, force = FALSE)"

# Add some libraries that we need
RUN apt-get -y update && apt-get -y install \
  libglpk40 \
  libxt6 \
  libxt-dev \
  libgeos-dev \
  libcairo2-dev \
  libhdf5-dev \
  cmake \
  patch \
  libgsl-dev \
  && apt-get clean \
  && rm -rf /tmp/* /var/tmp/*

# add nlopt (for muscat R package)
RUN git clone -b v2.7.1 --single-branch https://github.com/stevengj/nlopt.git \
  && cd nlopt \
  && mkdir build \
  && cd build \
  && cmake .. \
  && make \
  && make install \
  && cd ../.. \
  && rm -rf nlopt

# R package installation
RUN R -e "BiocManager::install(c('markdown', 'hdf5r', 'Seurat', 'SeuratObject', 'sparseMatrixStats', 'edgeR', 'apeglm', 'DESeq2', 'hypeR', 'patchwork', 'R.utils'))"
RUN R -e "remotes::install_github(repo = 'cancerbits/canceRbits', ref = '8f154d6')"
RUN R -e "devtools::install_github('ChristophH/riffle', ref = '8091c53')"
RUN R -e "devtools::install_github('const-ae/transformGamPoi', ref = 'f39f5e1')"
RUN R -e "BiocManager::install(c('caret', 'PRROC', 'scran', 'muscat', 'Rfast', 'GO.db'))"
RUN R -e "devtools::install_github('immunogenomics/presto', ref = '052085d')"

# miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh \
  && /bin/sh Miniconda3-py38_4.12.0-Linux-x86_64.sh -b -p /opt/miniconda3 \
  && rm -f Miniconda3-py38_4.12.0-Linux-x86_64.sh
RUN chmod -R a+rw /opt/miniconda3
ENV PATH=/opt/miniconda3/bin:$PATH

# install celltypist python package (this also installs scanpy, anndata and other python dependencies)
ENV CELLTYPIST_FOLDER="/opt/celltypist"
ENV PATH="$CELLTYPIST_FOLDER/bin:$PATH"
RUN python -m venv $CELLTYPIST_FOLDER && \
    pip install wheel --no-cache-dir &&\
    pip install celltypist --no-cache-dir && \
    celltypist --update-models

# more R packages
RUN R -e "BiocManager::install(c('MAST'))"
RUN R -e "devtools::install_github(repo = 'stephenturner/annotables', ref = '631423c')"
RUN R -e "BiocManager::install(c('scDD', 'fastcluster'))"
