#!/bin/bash

if [ $# -eq 0 ]
  then
    echo "Argument 1 should be a port number, argument 2 a string to use as password."
    exit 1
fi

# parse config parameters:
source bash/parse_yaml.sh
eval $(parse_yaml config.yaml CONF_)

# ids defined in image for the rstudio user
uid=1000
gid=1000
# subid ranges on host
subuidSize=$(( $(podman info --format "{{ range .Host.IDMappings.UIDMap }}+{{.Size }}{{end }}" ) - 1 ))
subgidSize=$(( $(podman info --format "{{ range .Host.IDMappings.GIDMap }}+{{.Size }}{{end }}" ) - 1 ))


podman run -d --rm \
  --name ${CONF_project_name}_${USER} \
  --uidmap $uid:0:1 --uidmap 0:1:$uid --uidmap $(($uid+1)):$(($uid+1)):$(($subuidSize-$uid)) \
  --gidmap $gid:0:1 --gidmap 0:1:$gid --gidmap $(($gid+1)):$(($gid+1)):$(($subgidSize-$gid)) \
  --group-add=keep-groups \
  -p ${1}:8787 \
  -e RUNROOTLESS=false \
  -e PASSWORD=${2} \
  -e UMASK=002 \
  -e TZ=Europe/Vienna \
  --volume=${CONF_project_root_host}:${CONF_project_root} \
  --volume=${CONF_out_root_host}:${CONF_out_root} \
  ${CONF_project_docker}
