#!/bin/bash

# help page
if [ "$#" != "3" ]; then
	echo "Usage: `basename $0` [FASTQ_R1] [FASTQ_R2] [TYPE]"
  	echo "Incorrect number of arguments"
	echo "- Required: FASTQ_R1 = FASTQ PAIRED-END READ 1"
	echo "- Required: FASTQ_R2 = FASTQ PAIRED-END READ 2"
	echo "- Required: TYPE = Target type (SNV | Fusion | DUX4 | FocDel | IGH)"
	exit 0
fi

# Setup all directories needed to run code and for results
BASE=`pwd`
TARGET_BASE=${BASE}/ALL_targets
VIRTUAL_ENV=${BASE}/.virtualenvs/km

if [ "$3" == "SNV" ]; then
	TARGET=${TARGET_BASE}/SNV
elif [ "$3" == "Fusion" ]; then
	TARGET=${TARGET_BASE}/Fusion
elif [ "$3" == "FocDel" ]; then
	TARGET=${TARGET_BASE}/focal_deletions
elif [ "$3" == "DUX4" ]; then
	TARGET=${TARGET_BASE}/DUX4
elif [ "$3" == "IGH" ]; then
	TARGET=${TARGET_BASE}/IGH_fusion 
else 
	echo "Not a valid argument"
	echo "[TYPE] needs to be one of (SNV | Fusion | DUX4 | FocDel | IGH)"
	exit 0
fi

# Directory setup
OUTPUTS=${BASE}/output
SAMPLE=`basename $1 _1.fastq.gz`

# km and jellyfish installed in virtual environment.
# must source the environment
source ${VIRTUAL_ENV}/bin/activate

# If statment that identifies the presence of `basename $1 .fastq.gz`.jf in the output dir 
if [ -f ${OUTPUTS}/${SAMPLE}/*.jf ]; then
	echo "Jellyfish Sample file exists"
else 
	# Make result directory if it doesnt exist
	mkdir -p ${OUTPUTS}/${SAMPLE}
	
	# Run jellyfish
	echo -e "generating count table for ${SAMPLE}"
	jellyfish count -m 31 -o ${OUTPUTS}/${SAMPLE}/countTable31.jf -s 100M -t 12 -C -L 2 <(zcat $1) <(zcat $2)
fi

# Run km find mutation
wc -l ${OUTPUTS}/${SAMPLE}/countTable31.jf
km find_mutation ${TARGET}/ ${OUTPUTS}/${SAMPLE}/countTable31.jf > ${OUTPUTS}/${SAMPLE}/${SAMPLE}_"$3".txt
