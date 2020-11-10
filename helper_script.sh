#!/bin/bash

### Run km on multiple samples against all target sets ###
  ### Also runs Rscript to collate & filter output ###

# help page
if [ "$#" != "1" ]; then
	echo "Usage: `basename $0` [SAMPLE_SHEET]"
  	echo "Incorrect number of arguments"
	echo "- Required: SAMPLE_SHEET = Tab-delimited Text file containing FASTQ R1 and R2 locations"
	echo ""
	echo "<READ1 PATH> <--TAB--> <READ2 PATH>"
	echo ""
	exit 0
fi

# Set base directory
BASE=`pwd`

### For loop; extracts sample data from txt file
while read line; do

  R1=$(echo $line | awk '{ print $1 }')
  R2=$(echo $line | awk '{ print $2 }')
  
  # Run_km
  bash run_km.sh ${R1} ${R2} Fusion
  bash run_km.sh ${R1} ${R2} DUX4
  bash run_km.sh ${R1} ${R2} SNV
  bash run_km.sh ${R1} ${R2} FocDel
  bash run_km.sh ${R1} ${R2} IGH
  
  # Directory setup for running Rscript
  OUTPUTS=${BASE}/output
  SAMPLE=`basename ${R1} _1.fastq.gz`
  
  # Run filter_km_output.R
  Rscript ${BASE}/bin/filter_km_output.R ${OUTPUTS}/${SAMPLE}
  
done < $1
  
  
