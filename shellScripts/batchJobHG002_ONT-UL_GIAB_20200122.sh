#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem-per-cpu 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12

##PATH=$R_HOME/bin:$PATH

snakemake --cores $(nproc) data/qualityControl/dataInput/HG002_ONT-UL_GIAB_20200122_fastqc.html data/qualityControl/readFilteringAndTrimming/HG002_ONT-UL_GIAB_20200122_fastqc.html data/qualityControl/assemblyAndErrorCorrection/HG002_ONT-UL_GIAB_20200122.txt data/qualityControl/alignment/HG002_ONT-UL_GIAB_20200122.txt data/variantAnnotation/structuralVariants/HG002_ONT-UL_GIAB_20200122.tsv data/variantAnnotation/singleNucleotidePolymorphisms/HG002_ONT-UL_GIAB_20200122.vcf data/variantCalling/transposableElements/HG002_ONT-UL_GIAB_20200122.vcf --rerun-incomplete

sbatch shellScripts/batchJobHG002_GRCh38_ONT-UL_GIAB_20200204.sh