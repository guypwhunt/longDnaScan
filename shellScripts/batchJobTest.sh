#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH -p cpu
#SBATCH --mem-per-cpu 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24

##PATH=$R_HOME/bin:$PATH

snakemake --cores $(nproc) data/qualityControl/dataInput/test_fastqc.html data/qualityControl/readFilteringAndTrimming/test_fastqc.html data/qualityControl/assemblyAndErrorCorrection/test.txt data/qualityControl/alignment/test.txt data/variantAnnotation/structuralVariants/test.tsv data/variantAnnotation/singleNucleotidePolymorphisms/test.vcf data/variantCalling/transposableElements/test.vcf

sbatch shellScripts/batchJobTest.sh