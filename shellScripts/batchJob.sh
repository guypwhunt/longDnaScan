#!/bin/sh
#SBATCH --time=47:00:00
#SBATCH -p cpu
#SBATCH --mem-per-cpu 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24

PATH=$R_HOME/bin:$PATH

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/longDnaScan/

conda activate longDnaScan 

snakemake -np data/qualityControl/dataInput/HG002_ONT-UL_GIAB_20200122_fastqc.html data/qualityControl/readFilteringAndTrimming/HG002_ONT-UL_GIAB_20200122_fastqc.html data/qualityControl/assemblyAndErrorCorrection/HG002_ONT-UL_GIAB_20200122.txt data/qualityControl/alignment/HG002_ONT-UL_GIAB_20200122.txt data/variantAnnotation/structuralVariants/HG002_ONT-UL_GIAB_20200122.tsv data/variantAnnotation/singleNucleotidePolymorphisms/HG002_ONT-UL_GIAB_20200122.vcf

snakemake --cores $(nproc) data/qualityControl/dataInput/HG002_ONT-UL_GIAB_20200122_fastqc.html data/qualityControl/readFilteringAndTrimming/HG002_ONT-UL_GIAB_20200122_fastqc.html data/qualityControl/assemblyAndErrorCorrection/HG002_ONT-UL_GIAB_20200122.txt data/qualityControl/alignment/HG002_ONT-UL_GIAB_20200122.txt data/variantAnnotation/structuralVariants/HG002_ONT-UL_GIAB_20200122.tsv data/variantAnnotation/singleNucleotidePolymorphisms/HG002_ONT-UL_GIAB_20200122.vcf