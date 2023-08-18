#!/bin/sh
#SBATCH --time=47:00:00
#SBATCH -p cpu
#SBATCH --mem-per-cpu 100G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20

PATH=$R_HOME/bin:$PATH

echo "started"

minimap2 --MD -ax map-ont -K 800000M -t 20 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/development/data/dataInput/referenceGenome/hg38.fa data/alignment/HG002_GRCh38_ONT-UL_GIAB_20200204.fastq | \
samtools view --threads 20 -bS | samtools sort --threads 20 -o data/alignment/HG002_GRCh38_ONT-UL_GIAB_20200204.bam

echo "alignment finished"

samtools index -@ 20 data/alignment/HG002_GRCh38_ONT-UL_GIAB_20200204.bam

echo "index finished"

samtools depth data/alignment/HG002_GRCh38_ONT-UL_GIAB_20200204.bam > data/qualityControl/alignment/HG002_GRCh38_ONT-UL_GIAB_20200204.txt

echo "depth finished"

python python/quality_control__aligment.py HG002_GRCh38_ONT-UL_GIAB_20200204

echo "completed"