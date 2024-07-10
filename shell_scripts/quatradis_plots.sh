#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20240704
#version    :1.0.0
#desc       :Script to run QuaTraDIS
#usage		:bash quatradis_plots.sh <directory containing fastq files> <reference file>
#===========================================================================================================
source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment

envpath="$(sudo find ~ -maxdepth 3 -name envs)"
if [ -d $envpath/quatradis ] 
then
	echo "QuaTraDIS conda env present" 
else
	echo "creating conda env: quatradis" 
	conda create -y -n quatradis python=3.11 -c bioconda -c conda-forge
    conda activate quatradis
    mamba install -c conda-forge -c bioconda bioconductor-edger -y
    mamba install -c conda-forge -c bioconda quatradis -y
    conda deactivate
fi

ls $1/*fastq.gz | sed -r 's/^.+\///' > $1/my_fastq_list.txt
mkdir -p \
	$1/plot_files/condition_files \
	$1/plot_files/control_files \
	$1/analysis_output

conda activate quatradis
tradis pipeline multiple \
	--output_dir $1/plot_files \
	$1/my_fastq_list.txt \
	$2
