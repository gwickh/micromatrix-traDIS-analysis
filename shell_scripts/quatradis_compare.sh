#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20240704
#version    :1.0.0
#desc       :Script to run QuaTraDIS
#usage		:bash quatradis_compare.sh <directory containing fastq files> <reference file>
#===========================================================================================================
source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment

envpath="$(sudo find ~ -maxdepth 3 -name envs)"
if [ -d $envpath/quatradis ] 
then
	echo "QuaTraDIS conda env present" 
else
	echo "creating conda env: quatradis" 
	conda create -y -n quatradis -c bioconda -c conda-forge
    conda activate quatradis
    mamba install -c conda-forge -c bioconda bioconductor-edger -y
    mamba install -c conda-forge -c bioconda quatradis -y
    conda deactivate
fi

conda activate quatradis
tradis -V
tradis pipeline compare \
	--output_dir $1/analysis_output \
	--annotations $2 \
	--condition_files $3/* \
	--control_files $4/*
conda deactivate