#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20240704
#version    :1.0.0
#desc       :Script to run QuaTraDIS
#usage		:bash quatradis.sh <directory containing fastq files>
#===========================================================================================================
source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment

alert_banner() {
	echo ""
	echo "####################################################################################################"
	echo ""
	echo "$alert"	
	echo ""
	echo "####################################################################################################"
	echo "" 
}


envpath="$(sudo find ~ -maxdepth 3 -name envs)"
if [ -f $envpath/bin/quatradis ] 
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

conda activate quatradis
tradis --help