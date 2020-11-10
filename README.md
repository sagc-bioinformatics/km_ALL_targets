
# `km` analysis with SAHMRI ALL targets

Jimmy Breen : jimmy.breen@sahrmi.com
Jacqueline Rehn : jacqueline.rehn@sahmri.com

Utilises k-mer variant detection software 'km' developed by a canadian bioinformatic group (https://github.com/iric-soft/km) and input targets designed for all gene fusions/variants

## Install jellyfish and km

1. To install and run the km script we'll need to create a conda environment to install km and jellyfish in. Jellyfish needs additional python bindings to work with km, so we have created an install procedure to follow.

```
conda create -n km python=3.7 virtualenv -y
conda activate km
```

2. Next we will ....
```
git clone []
```

3. Now we can install the km/jellyfish components using this `install.sh` script. 

```
bash install.sh
```


Install R libraries:
- tidyverse
- Biostrings
- here

## Quick Start



To run

```
(km)$ bash run_km.sh
Usage: run_km.sh [FASTQ_R1] [FASTQ_R2] [TYPE]
Incorrect number of arguments
- Required: FASTQ_R1 = FASTQ PAIRED-END READ 1
- Required: FASTQ_R2 = FASTQ PAIRED-END READ 2
- Required: TYPE = Target type (SNV | Fusion | DUX4 | FocDel | IGH)
```



```
(km)$ bash run_km.sh AYA-0682-201240_1.fastq.gz AYA-0682-201240_2.fastq.gz SNV
generating count table for AYA-0682-201240
936314 /homes/jimmy.breen/bioinfoCore/Cancer/ALL/km_test_jacqui/output/AYA-0682-201240/countTable31.jf

```

## Notes

- There will be `PathQuant.py:125: RuntimeWarning: invalid value encountered in greater_equal` warnings on km but these can be ignored
- Summary Rscript requires `tidyverse`, `Biostrings` and `here` r packages and are installed at the start of the script. If there are any issues, it might be worth running on your local copy of R (rather than conda env).
