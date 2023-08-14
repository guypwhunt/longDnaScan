#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem 200G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

source activate miniStructuralVariantAnalysis 

## Concatenate all fastq files into a single file
cat /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/*.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq

## Run fastQC quality check
fastqc /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq --outdir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/FastQC -t 2
##fastqc /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/*.fastq --outdir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/FastQC

## Run pycoQC quality check
pycoQC -f /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/sequencing_summary.txt -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/PycoQC/pycoQC.html

## Run MinION_QC quality check
 MinIONQC.R -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/sequencing_summary.txt -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/QualityControl/MinION_QC

## PoreChop Adapter Removal
porechop -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/AdapterRemoval/PoreChop/PoreChop.fastq --discard_middle

## nanofilt Trim Reads and Filtering
NanoFilt -l 500 --headcrop 10 </scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/AdapterRemoval/PoreChop/PoreChop.fastq> /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq

## Minimap & Miniasm Genome Assembly
minimap2 -x ava-ont /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq | gzip -1 > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz

miniasm -f /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa

awk '/^S/{print ">"$2"\n"$3}' /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

## Flye Genome Assembly
flye --nano-raw /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq --genome-size 1m --out-dir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye

## Shasta Genome Assembly (This Fails due to memory limit)
cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/
sed -n '1~4s/^@/>/p;2~4p' /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/all_guppy.fasta

/scratch/prj/sgdp_nanopore/software/shasta-Linux-0.11.1 --input  /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/all_guppy.fasta --config Nanopore-UL-May2022 --threads 10 

## Assembly-Stats Genome Assembly Quality Control
assembly-stats /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/AssemblyStats/MinimapAndMiniasm.txt
assembly-stats /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye/assembly.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/AssemblyStats/Flye.txt
assembly-stats /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/AssemblyStats/Shasta.txt

## Dnadiff (Reference Genome Alignment) Genome Assembly Quality Control
cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/Dnadiff/MinimapMiniasm
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/Dnadiff/Flye
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye/assembly.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssemblyQualityControl/Dnadiff/Shasta
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta

## Racon Error Correction
minimap2 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/minimap.racon.paf
racon /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/minimap.racon.paf /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

minimap2 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye/assembly.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/flye.racon.paf
racon /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/flye.racon.paf /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Flye/assembly.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/flye.racon.consensus.fasta

# This Errors
minimap2 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/shasta.racon.paf
racon /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/shasta.racon.paf /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta -f > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/shasta.racon.consensus.fasta 

## Minipolish Error Correction
minipolish -t 10 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/minipolish/minimapMiniasm.gfa
awk '/^S/{print ">"$2"\n"$3}' /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/minipolish/minimapMiniasm.gfa > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/minipolish/minimapMiniasm.fasta

## Medaka Error Correction (couldn't install Medaka)
#### Nanopolish Error Correction (Need to Include) ####

## Pilon Error Correction 
bwa index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta
bwa index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/flye.racon.consensus.fasta
#bwa index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/shasta.racon.consensus.fasta 

bwa mem -t 10 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta -x ont2d /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.sam
bwa mem -t 10 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/flye.racon.consensus.fasta -x ont2d /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.sam
#bwa mem -t 10 /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/shasta.racon.consensus.fasta -x ont2d /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/Guppy/all_guppy.fastq > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.sam

samtools view -Sb /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.bam
samtools view -Sb /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.bam
#samtools view -Sb /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.bam

samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.bam
samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.bam
#samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.bam

samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam
samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.sorted.bam
#samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.sorted.bam

java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta --bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam --threads 10 --outdir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm
java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/flye.racon.consensus.fasta --bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye_racon.sorted.bam --threads 10 --outdir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye
#java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/shasta.racon.consensus.fasta --bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta_racon.sorted.bam --threads 10 --outdir /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta

## Error Correction Align to Reference Genome
cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrectionQualityControl/Dnadiff/Racon
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrectionQualityControl/Dnadiff/Minipolish
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/minipolish/minimapMiniasm.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/miniasm
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm/pilon.fasta

cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/flye
dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye/pilon.fasta

#cd /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/shasta
#dnadiff -p dnadiff /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/shasta/pilon.fasta


## Sniffles Variant Calling
minimap2 --MD -a  /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.sam
minimap2 --MD -a  /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/minipolish/minimapMiniasm.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.sam
minimap2 --MD -a  /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/miniasm/pilon.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.sam
minimap2 --MD -a  /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/ErrorCorrection/Pilon/flye/pilon.fasta > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.sam

# Convert to bam file
samtools view -bS /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.bam
samtools view -bS /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.bam
samtools view -bS /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.bam
samtools view -bS /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.sam > /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.bam

# Sort the bam file
samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.bam
samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.bam
samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.bam
samtools sort -o /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.sorted.bam /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.bam

# create an index file
samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.sorted.bam
samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.sorted.bam
samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam
samtools index /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.sorted.bam

sniffles -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.sorted.bam -v /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/racon.vcf --allow-overwrite
sniffles -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.sorted.bam -v /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/minipolish.vcf --allow-overwrite
sniffles -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam -v /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/miniasm_pilon.vcf --allow-overwrite
sniffles -i /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.sorted.bam -v /scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Output/VariantCalling/Sniffles/flye_pilon.vcf --allow-overwrite

source deactivate miniStructuralVariantAnalysis 


source activate tldr

