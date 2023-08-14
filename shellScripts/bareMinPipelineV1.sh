#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH -p cpu
#SBATCH --mem 800G
#SBATCH --nodes 4

PATH=$R_HOME/bin:$PATH

source activate miniStructuralVariantAnalysis 

inputDirectory='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/MNDTest/Input/Guppy/'
outputDirectory='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/MNDTest/'
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

cp -a /scratch/prj/sgdp_nanopore/Projects/BR22_00043_MND_sequencing/A0690_02/20230413_1551_1A_PAM24731_67c4936d/fastq_pass/ $inputDirectory

gunzip "${inputDirectory}"fastq_pass/*.fastq.gz

## Concatenate all fastq files into a single file
cat "${inputDirectory}"fastq_pass/*.fastq > "${outputDirectory}"Output/Guppy/all_guppy.fastq

rm -r "${inputDirectory}"fastq_pass/*.fastq

### Quality Control (QC) ###
## Run fastQC quality check
fastqc "${outputDirectory}"Output/Guppy/all_guppy.fastq --outdir "${outputDirectory}"Output/QualityControl/FastQC -t $threads

echo ""
echo "fastqc complete"
echo ""

### Read Filtering & Trimming ###
## PoreChop Adapter Removal
porechop -i "${outputDirectory}"Output/Guppy/all_guppy.fastq -o "${outputDirectory}"Output/AdapterRemoval/PoreChop/PoreChop.fastq --discard_middle

echo ""
echo "porechop complete"
echo ""

rm "${outputDirectory}"Output/Guppy/all_guppy.fastq

## nanofilt Trim Reads and Filtering
NanoFilt -l 500 --headcrop 10 <"${outputDirectory}"Output/AdapterRemoval/PoreChop/PoreChop.fastq> "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq

echo ""
echo "NanoFilt complete"
echo ""

rm "${outputDirectory}"Output/AdapterRemoval/PoreChop/PoreChop.fastq

### Assembly ###
## Minimap & Miniasm Genome Assembly
minimap2 -x ava-ont "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq | gzip -1 > "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz

miniasm -f "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/minimap.paf.gz > "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa

awk '/^S/{print ">"$2"\n"$3}' "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.gfa > "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta


echo ""
echo "Assembly complete"
echo ""

### Assembly QC ###
## Assembly-Stats Genome Assembly Quality Control
assembly-stats "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > "${outputDirectory}"Output/GenomeAssemblyQualityControl/AssemblyStats/MinimapAndMiniasm.txt

## Dnadiff (Reference Genome Alignment) Genome Assembly Quality Control
cd "${outputDirectory}"Output/GenomeAssemblyQualityControl/Dnadiff/MinimapMiniasm
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta

## Racon Error Correction
minimap2 "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq > "${outputDirectory}"Output/ErrorCorrection/Racon/minimap.racon.paf
racon "${outputDirectory}"Output/ReadTrimmingAndFiltering/NanoFilt/NanoFilt.fastq "${outputDirectory}"Output/ErrorCorrection/Racon/minimap.racon.paf "${outputDirectory}"Output/GenomeAssembly/MinimapMiniasm/miniasm.fasta > "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

### Assembly Error Correction ###
## Pilon Error Correction 
bwa index "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

bwa mem -t $threads "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta -x ont2d "${outputDirectory}"Output/Guppy/all_guppy.fastq > "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sam

samtools view -Sb "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sam > "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.bam

samtools sort -o "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.bam

samtools index "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam

java -Xmx200G -jar /users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar --genome "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta --bam "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm_racon.sorted.bam --threads 10 --outdir "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm

### Assembly Error Correction QC ###
## Error Correction Align to Reference Genome
cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Racon
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta

##cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Minipolish
##dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/minipolish/minimapMiniasm.fasta

cd "${outputDirectory}"Output/ErrorCorrectionQualityControl/Dnadiff/Pilon/miniasm
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm/pilon.fasta

### Alignment ###
minimap2 --MD -a  "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Racon/miniasm.racon.consensus.fasta > "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sam
minimap2 --MD -a  "${genomeReferencePath}" "${outputDirectory}"Output/ErrorCorrection/Pilon/miniasm/pilon.fasta > "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sam

# Convert to bam file
samtools view -bS "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sam > "${outputDirectory}"Output/VariantCalling/Sniffles/racon.bam
samtools view -bS "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.sam > "${outputDirectory}"Output/VariantCalling/Sniffles/miniasm_pilon.bam

# Sort the bam file
samtools sort -o "${outputDirectory}"Output/VariantCalling/Sniffles/racon.sorted.bam "${outputDirectory}"Output/VariantCall