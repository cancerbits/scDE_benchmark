#!/bin/bash

# Note, that this should be executed on the host, not within the project container

# parse config parameters:
source bash/parse_yaml.sh
eval $(parse_yaml config.yaml CONF_)

outdir=${CONF_out_root_host}/squair_et_al_data
mkdir -p ${outdir}
cd ${outdir}

wget -O squair_et_al_sc_rnaseq.tar.gz  https://zenodo.org/record/7659806/files/sc_rnaseq.tar.gz?download=1
wget -O squair_et_al_bulk_rnaseq.tar.gz  https://zenodo.org/record/7659806/files/bulk_rnaseq.tar.gz?download=1

tar -xvzf ${outdir}/squair_et_al_sc_rnaseq.tar.gz
tar -xvzf ${outdir}/squair_et_al_bulk_rnaseq.tar.gz

rm -rf squair_et_al_sc_rnaseq.tar.gz
rm -rf squair_et_al_bulk_rnaseq.tar.gz
