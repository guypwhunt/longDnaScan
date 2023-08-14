#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem 200G
##SBATCH --nodes 11

PATH=$R_HOME/bin:$PATH

source activate miniStructuralVariantAnalysis 

inputDirectory='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/Guppy/'
outputDirectory='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/'
genomeReferencePath='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/WorkflowTest/Input/ReferenceGenome/chr17.fasta'
genomeReferencePath='/scratch/prj/herv_project/UCSC_hg38/hg38.fa'
threads=100

cd $outputDirectory

mkdir -p Output/AdapterRemoval/PoreChop
mkdir -p Output/ErrorCorrection/Medaka
mkdir -p Output/ErrorCorrection/minipolish
mkdir -p Output/ErrorCorrection/Pilon/flye
mkdir -p Output/ErrorCorrection/Pilon/miniasm
mkdir -p Output/ErrorCorrection/Pilon/shasta
mkdir -p Output/ErrorCorrection/Racon
mkdir -p Output/ErrorCorrectionQualityControl/Dnadiff/Minipolish
mkdir -p Output/ErrorCorrectionQualityControl/Dnadiff/Pilon
mkdir -p Output/ErrorCorrectionQualityControl/Dnadiff/Racon
mkdir -p Output/GenomeAssembly/Flye
mkdir -p Output/GenomeAssembly/MinimapMiniasm
mkdir -p Output/GenomeAssembly/Shasta
mkdir -p Output/GenomeAssemblyQualityControl/Dnadiff/Minipolish
mkdir -p Output/GenomeAssemblyQualityControl/Dnadiff/Pilon
mkdir -p Output/GenomeAssemblyQualityControl/Dnadiff/Racon
mkdir -p Output/GenomeAssemblyQualityControl/Dnadiff/Flye
mkdir -p Output/GenomeAssemblyQualityControl/AssemblyStats

### Update to basecalling
mkdir -p Output/Guppy

mkdir -p Output/QualityControl/FastQC
mkdir -p Output/QualityControl/MinION_QC
mkdir -p Output/QualityControl/PycoQC
mkdir -p Output/ReadTrimmingAndFiltering/NanoFilt
mkdir -p Output/VariantCalling/Sniffles
mkdir -p Output/VariantCalling/TLDR/Racon
mkdir -p Output/VariantCalling/TLDR/Pilon


## Concatenate all fastq files into a single file
cat "${inputDirectory}"*.fastq > "${outputDirectory}"Output/Guppy/all_guppy.fastq

### Quality Control (QC) ###
## Run fastQC quality check
fastqc "${outputDirectory}"Output/Guppy/all_guppy.fastq --outdir "${outputDirectory}"Output/QualityControl/FastQC -t $threads

## Run pycoQC quality check
pycoQC -f "${outputDirectory}"Input/Guppy/sequencing_summary.txt -o "${outputDirectory}"Output/QualityControl/PycoQC/pycoQC.html

## Run MinION_QC quality check
 MinIONQC.R -i "${outputDirectory}"Input/Guppy/sequencing_summary.txt -o "${outputDirectory}"Output/QualityControl/MinION_QC

### Read Filtering & Trimming ###
## PoreChop Adapter Removal
porechop -i "${outputDirectory}"Output/Guppy/all_guppy.fastq -o "${outputDirectory}"Output/AdapterRemoval/PoreChop/PoreChop.fastq --discard_middle

## nanofilt Trim Reads and Filtering
NanoFilt -l 500 --headcrop 10 <"${outputDirectory}"Output/AdapterRemoval/PoreChop/PoreChop.fastq> "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq

### Assembly ###
## Minimap & Miniasm Genome Assembly
minimap2 -x ava-ont "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq | gzip -1 > "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz

miniasm -f "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz > "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa

awk '/^S/{print ">"$2"\n"$3}' "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa > "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

## Flye Genome Assembly 
## This is for raw reads
## flye --nano-raw "${outputDirectory}"Output/Guppy/all_guppy.fastq --genome-size 1m --out-dir "${outputDirectory}"Output/GenomeAssembly/Flye
## This is for corrected reads
##flye --nano-corr "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq --genome-size 1m --out-dir "${outputDirectory}"Output/GenomeAssembly/Flye

## Shasta Genome Assembly (This Fails due to memory limit)
##cd "${outputDirectory}"Output/GenomeAssembly/Shasta/
##sed -n '1~4s/^@/>/p;2~4p' "${outputDirectory}"Output/Guppy/all_guppy.fastq > "${outputDirectory}"Output/GenomeAssembly/Shasta/all_guppy.fasta
##sed -n '1~4s/^@/>/p;2~4p' "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > "${outputDirectory}"Output/GenomeAssembly/Shasta/all_guppy.fasta

##/scratch/prj/sgdp_nanopore/software/shasta-Linux-0.11.1 --input  "${outputDirectory}"Output/GenomeAssembly/Shasta/all_guppy.fasta --config Nanopore-UL-May2022 --threads $threads --memoryBacking 2M

### Assembly QC ###
## Assembly-Stats Genome Assembly Quality Control
assembly-stats "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > "${outputDirectory}"Output/GenomeAssemblyQualityControl/AssemblyStats/MinimapAndMiniasm.txt
##assembly-stats "${outputDirectory}"Output/GenomeAssembly/Flye/assembly.fasta > "${outputDirectory}"Output/GenomeAssemblyQualityControl/AssemblyStats/Flye.txt
##assembly-stats "${outputDirectory}"Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta > "${outputDirectory}"Output/GenomeAssemblyQualityControl/AssemblyStats/Shasta.txt

## Dnadiff (Reference Genome Alignment) Genome Assembly Quality Control
cd "${outputDirectory}"Output/GenomeAssemblyQualityControl/Dnadiff/MinimapMiniasm
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

### Continue HERE!!!!
##cd "${outputDirectory}"Output/GenomeAssemblyQualityControl/Dnadiff/Flye
##dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/GenomeAssembly/Flye/assembly.fasta

##cd "${outputDirectory}"Output/GenomeAssemblyQualityControl/Dnadiff/Shasta
##dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta

## Racon Error Correction
minimap2 "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > "${outputDirectory}"Output/ErrorCorrection/Racon/minimap.racon.paf
racon "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/ErrorCorrection/Racon/minimap.racon.paf "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

##minimap2 "${outputDirectory}"Output/GenomeAssembly/Flye/assembly.fasta "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > "${outputDirectory}"Output/ErrorCorrection/Racon/flye.racon.paf
##racon "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/ErrorCorrection/Racon/flye.racon.paf "${outputDirectory}"Output/GenomeAssembly/Flye/assembly.fasta > "${outputDirectory}"Output/ErrorCorrection/Racon/flye.racon.consensus.fasta

## This Errors
##minimap2 "${outputDirectory}"Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > "${outputDirectory}"Output/ErrorCorrection/Racon/shasta.racon.paf
##racon "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/ErrorCorrection/Racon/shasta.racon.paf "${outputDirectory}"Output/GenomeAssembly/Shasta/ShastaRun/Assembly.fasta -f > "${outputDirectory}"Output/shasta.racon.consensus.fasta 

## Minipolish Error Correction
##minipolish -t $threads "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa > "${outputDirectory}"Output/ErrorCorrection/minipolish/minimapMiniasm.gfa
##awk '/^S/{print ">"$2"\n"$3}' "${outputDirectory}"Output/ErrorCorrection/minipolish/minimapMiniasm.gfa > "${outputDirectory}"Output/ErrorCorrection/minipolish/minimapMiniasm.fasta

## Medaka Error Correction (couldn't install Medaka)
#### Nanopolish Error Correction (Need to Include) ####

### Assembly Error Correction ###
## Pilon Error Correction 
bwa index "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta
##bwa index "${outputDirectory}"Output/ErrorCorrection/Racon/flye.racon.consensus.fasta
#bwa index "${outputDirectory}"Output/ErrorCorrection/Racon/shasta.racon.consensus.fasta 

bwa mem -t $threads "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta -x ont2d "${outputDirectory}"Output/Guppy/all_guppy.fastq > "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sam
##bwa mem -t $threads "${outputDirectory}"Output/ErrorCorrection/Racon/flye.racon.consensus.fasta -x ont2d "${outputDirectory}"Output/Guppy/all_guppy.fastq > "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.sam
#bwa mem -t $threads "${outputDirectory}"Output/ErrorCorrection/Racon/shasta.racon.consensus.fasta -x ont2d "${outputDirectory}"Output/Guppy/all_guppy.fastq > "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.sam

samtools view -Sb "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sam > "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.bam
##samtools view -Sb "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.sam > "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.bam
#samtools view -Sb "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.sam > "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.bam

samtools sort -o "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.bam
##samtools sort -o "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.sorted.bam "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.bam
#samtools sort -o "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.sorted.bam "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.bam

samtools index "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam
##samtools index "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.sorted.bam
#samtools index "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.sorted.bam

java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta --bam "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam --threads 10 --outdir "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm
##java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome "${outputDirectory}"Output/ErrorCorrection/Racon/flye.racon.consensus.fasta --bam "${outputDirectory}"Output/ErrorCorrection/Pilon/flye_racon.sorted.bam --threads 10 --outdir "${outputDirectory}"Output/ErrorCorrection/Pilon/flye
#java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome "${outputDirectory}"Output/ErrorCorrection/Racon/shasta.racon.consensus.fasta --bam "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta_racon.sorted.bam --threads 10 --outdir "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta

### Assembly Error Correction QC ###
## Error Correction Align to Reference Genome
cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Racon
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

##cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Minipolish
##dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/minipolish/minimapMiniasm.fasta

cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/miniasm
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm/pilon.fasta

##cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/flye
##dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/flye/pilon.fasta

#cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/shasta
#dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/shasta/pilon.fasta

### Alignment ###
## Sniffles Variant Calling
minimap2 --MD -a  "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta > "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sam
##minimap2 --MD -a  "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/minipolish/minimapMiniasm.fasta > "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.sam
minimap2 --MD -a  "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm/pilon.fasta > "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sam
##minimap2 --MD -a  "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/flye/pilon.fasta > "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.sam

# Convert to bam file
samtools view -bS "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sam > "${outputDirectory}"Output/VariantCalling/Sniffles/racon.bam
##samtools view -bS "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.sam > "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.bam
samtools view -bS "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sam > "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.bam
##samtools view -bS "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.sam > "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.bam

# Sort the bam file
samtools sort -o "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sorted.bam "${outputDirectory}"Output/VariantCalling/Sniffles/racon.bam
##samtools sort -o "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.sorted.bam "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.bam
samtools sort -o "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.bam
##samtools sort -o "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.sorted.bam "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.bam

# create an index file
samtools index "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sorted.bam
##samtools index "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.sorted.bam
samtools index "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam
##samtools index "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.sorted.bam

sniffles -i "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sorted.bam -v "${outputDirectory}"Output/VariantCalling/Sniffles/racon.vcf --allow-overwrite
##sniffles -i "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.sorted.bam -v "${outputDirectory}"Output/VariantCalling/Sniffles/minipolish.vcf --allow-overwrite
sniffles -i "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam -v "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.vcf --allow-overwrite
##sniffles -i "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.sorted.bam -v "${outputDirectory}"Output/VariantCalling/Sniffles/flye_pilon.vcf --allow-overwrite

source deactivate miniStructuralVariantAnalysis 

source activate tldr
cd "${outputDirectory}"Output/VariantCalling/TLDR/Racon
tldr -b "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sorted.bam -e /users/k20064105/tldr/ref/teref.ont.human.fa -r "${genomeReferencePath}" --color_consensus

cd "${outputDirectory}"Output/VariantCalling/TLDR/Pilon
tldr -b "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sorted.bam -e /users/k20064105/tldr/ref/teref.ont.human.fa -r "${genomeReferencePath}" --color_consensus

