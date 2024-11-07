#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20241106
#version    :1.0.0
#desc       :Script to run albatradis
#usage		:bash albatradis.sh <EMBL file> <dir containing plotfiles>
#===========================================================================================================
echo "finding conda base environment"
source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment
echo "conda base environment located"

echo "finding albatradis conda environment"
envpath="$(sudo find ~ -maxdepth 3 -name envs)" #check if albatradis conda env exists
if [ -d $envpath/albatradis ] 
then
	echo "albatradis conda env present" 
else
	echo "creating conda env: quatradis" 
	conda create -y -n albatradis -c bioconda -c conda-forge albatradis
fi

conda activate albatradis
for k in $2/*
do
    #check if text file with plotfile filenames exists
    if [ -f $k/$(basename $k).txt ]
    then
        #if so, then delete
        rm $k/$(basename $k).txt
    fi
    #create text file with plotfile filenames
    ls -d $k/*plot.gz > $k/.$(basename $k).txt && mv $k/.$(basename $k).txt $k/$(basename $k).txt
    #pass filenames to albatradis
    echo "running albatradis on $k" 
    xargs \
        -a $k/$(basename $k).txt \
        albatradis -a $1
    mv ./output $k/ #mv output from current dir to working dir
    #prepend dirname to output filenames
    for j in $k/output/*
    do 
        mv $j $(dirname $j)/$(basename $k)_$(basename $j)
    done  
done
conda deactivate