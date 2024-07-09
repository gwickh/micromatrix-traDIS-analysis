#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20240709
#version    :1.0.0
#desc       :Script to mount HPC volume
#usage		:bash mount_volume.sh <directory name>
#===========================================================================================================

if -d ~/$1
    echo "$1 exists, mounting //smb.qib-hpc-data.ciscloud/research-groups/Mark-Webber/ to $1"
else
    echo "creating $1"
    echo "$1 exists, mounting //smb.qib-hpc-data.ciscloud/research-groups/Mark-Webber/ to $1"

sudo \
    mount.cifs \
    //smb.qib-hpc-data.ciscloud/research-groups/Mark-Webber/ \
    ~/$1 \
    -o domain=nr4,user=wickhamg,uid=$(id -u),gid=$(id -g),mfsymlinks