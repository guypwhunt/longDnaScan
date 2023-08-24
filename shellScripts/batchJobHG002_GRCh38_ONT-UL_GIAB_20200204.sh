#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem-per-cpu 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 12

##PATH=$R_HOME/bin:$PATH

snakemake --cores $(nproc) data/qualityControl/dataInput/HG002_GRCh38_ONT-UL_GIAB_20200204_fastqc.html data/qualityControl/readFilteringAndTrimming/HG002_GRCh38_ONT-UL_GIAB_20200204_fastqc.html data/qualityControl/assemblyAndErrorCorrection/HG002_GRCh38_ONT-UL_GIAB_20200204.txt data/qualityControl/alignment/HG002_GRCh38_ONT-UL_GIAB_20200204.txt data/variantAnnotation/structuralVariants/HG002_GRCh38_ONT-UL_GIAB_20200204.tsv data/variantAnnotation/singleNucleotidePolymorphisms/HG002_GRCh38_ONT-UL_GIAB_20200204.vcf data/variantCalling/transposableElements/HG002_GRCh38_ONT-UL_GIAB_20200204.vcf --rerun-incomplete

sbatch shellScripts/batchJobHG002_ONT-UL_GIAB_20200122.sh