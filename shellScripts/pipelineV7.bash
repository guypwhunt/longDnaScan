#!/bin/sh
#SBATCH --time=48:00:00
#SBATCH -p cpu
#SBATCH --mem-per-cpu=100G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=5


# Notes to self
# Add ability to keep or delete intermediate files
# Add ability to select step (based on functionality)

PATH=$R_HOME/bin:$PATH

source activate miniStructuralVariantAnalysis 

inputDirectory='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/MNDTest/Input/Guppy/'
#inputDirectory='/scratch/prj/sgdp_nanopore/Projects/BR22_00043_MND_sequencing/A0690_02/20230413_1551_1A_PAM24731_67c4936d/fastq_pass/'
outputDirectory='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/MNDTest/Output/'
genomeReferencePath='/scratch/prj/sgdp_nanopore/Projects/09_structural_variant_benchmarking_dataset/MNDTest/Input/ReferenceGenome/chr17.fasta'
#genomeReferencePath='/scratch/prj/herv_project/UCSC_hg38/hg38.fa'
pilonPath='/users/k20064105/Miniconda3/pkgs/pilon-1.24-hdfd78af_0/share/pilon-1.24-0/pilon.jar'
threads=20
memoryApplocation="800G"

cd $outputDirectory

mkdir -p qualityControl/preReadFilteringAndTrimming/
mkdir -p qualityControl/postReadFilteringAndTrimming/
mkdir -p qualityControl/preAssemblyErrorCorrection/
mkdir -p qualityControl/postAssemblyErrorCorrection/
mkdir -p qualityControl/alignment/
mkdir -p readFilteringAndTrimming/
mkdir -p assembly/
mkdir -p input/
mkdir -p alignment/
mkdir -p variantCalling/structuralVariants/
mkdir -p variantCalling/transposons/
mkdir -p variantCalling/structuralVariantsAnnotation/
mkdir -p variantCalling/transposonsAnnotation/

### Data Input ###
cp -a "${inputDirectory}". "${outputDirectory}"input/
## Concatenate all fastq files into a single file
gunzip -q "${outputDirectory}"input/*.fastq
cat "${outputDirectory}"input/*.fastq > "${outputDirectory}"readFilteringAndTrimming/all.fastq
gzip -q "${outputDirectory}"readFilteringAndTrimming/all.fastq
rm -r "${outputDirectory}"input
echo ""
echo "Data Input Complete"
echo ""

### Quality Control (QC) ###
## Run fastQC quality check
fastqc "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz --outdir "${outputDirectory}"qualityControl/preReadFilteringAndTrimming/ -t $threads
echo ""
echo "Pre-Read Filtering And Trimming Quality Control Complete"
echo ""

### Read Filtering & Trimming ###
## PoreChop Adapter Removal
porechop -i "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz -o "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz --discard_middle -t $threads
## nanofilt Trim Reads and Filtering
gunzip -q "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz
NanoFilt -l 500 --headcrop 10 --maxlength 100000 <"${outputDirectory}"readFilteringAndTrimming/all.fastq> "${outputDirectory}"readFilteringAndTrimming/a.fastq
rm "${outputDirectory}"readFilteringAndTrimming/all.fastq
mv "${outputDirectory}"readFilteringAndTrimming/a.fastq "${outputDirectory}"readFilteringAndTrimming/all.fastq
gzip -q "${outputDirectory}"readFilteringAndTrimming/all.fastq
echo ""
echo "Read Filtering And Trimming Complete"
echo ""

## Re-Run fastQC quality check
fastqc "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz --outdir "${outputDirectory}"qualityControl/postReadFilteringAndTrimming/ -t $threads
echo ""
echo "Post-Read Filtering And Trimming Quality Control Complete"
echo ""

### Assembly ###
#### FLYE TEST #####
flye --nano-raw "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz --out-dir "${outputDirectory}"assembly -t $threads --resume
mv "${outputDirectory}"assembly/assembly.fasta "${outputDirectory}"assembly/all.fasta
gzip -q "${outputDirectory}"assembly/all.fasta
rm "${outputDirectory}"assembly/all.fasta

## Minimap & Miniasm Genome Assembly
#minimap2 -x ava-ont -K $memoryApplocation -t $threads "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz | gzip -1 > "${outputDirectory}"assembly/all.paf.gz
# minimap2 -x ava-pb -K $memoryApplocation -t $threads "${outputDirectory}"readFilteringAndTrimming/all.fastq "${outputDirectory}"readFilteringAndTrimming/all.fastq | gzip -1 > "${outputDirectory}"assembly/all.paf.gz
#miniasm -f "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz "${outputDirectory}"assembly/all.paf.gz > "${outputDirectory}"assembly/all.gfa
#awk '/^S/{print ">"$2"\n"$3}' "${outputDirectory}"assembly/all.gfa > "${outputDirectory}"assembly/all.fasta
#gzip -q "${outputDirectory}"assembly/all.fasta
#rm "${outputDirectory}"assembly/all.gfa
#echo ""
#echo "Assembly Complete"
#echo ""

### Assembly QC ###
## Assembly-Stats Genome Assembly Quality Control
gunzip "${outputDirectory}"assembly/all.fasta
assembly-stats "${outputDirectory}"assembly/all.fasta > "${outputDirectory}"qualityControl/preAssemblyErrorCorrection/assemblyStats.txt
## Dnadiff (Reference Genome Alignment) Genome Assembly Quality Control
cd "${outputDirectory}"qualityControl/preAssemblyErrorCorrection
dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"assembly/all.fasta
gzip -q "${outputDirectory}"assembly/all.fasta
echo ""
echo "Pre-Assembly Error Correction Quality Control Complete"
echo ""

### Assembly Error Correction ###
## Racon Error Correction
#minimap2 -K 800000M -t $threads "${outputDirectory}"assembly/all.fasta.gz "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz > "${outputDirectory}"assembly/all.paf
#rm "${outputDirectory}"assembly/all.paf.gz
#gzip -q "${outputDirectory}"assembly/all.paf
#racon -t $threads "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz "${outputDirectory}"assembly/all.paf.gz "${outputDirectory}"assembly/all.fasta.gz > "${outputDirectory}"assembly/a.fasta
#mv "${outputDirectory}"assembly/a.fasta "${outputDirectory}"assembly/all.fasta
#rm "${outputDirectory}"assembly/all.fasta.gz
#gzip -q "${outputDirectory}"assembly/all.fasta
### Pilon Error Correction 
#bwa index "${outputDirectory}"assembly/all.fasta.gz
#bwa mem -t $threads "${outputDirectory}"assembly/all.fasta.gz -x ont2d "${outputDirectory}"readFilteringAndTrimming/all.fastq.gz > "${outputDirectory}"assembly/all.sam
#gzip -q "${outputDirectory}"assembly/all.sam
#samtools view --threads $threads -Sb "${outputDirectory}"assembly/all.sam.gz > "${outputDirectory}"assembly/all.bam
#samtools sort --threads $threads -o "${outputDirectory}"assembly/all.bam "${outputDirectory}"assembly/all.bam
#samtools index -@ $threads "${outputDirectory}"assembly/all.bam

#java -Xmx"${memoryApplocation}" -jar "${pilonPath}" --genome "${outputDirectory}"assembly/all.fasta.gz --nanopore "${outputDirectory}"assembly/all.bam --threads $threads --outdir "${outputDirectory}"assembly
#java -Xmx"${memoryApplocation}" -jar "${pilonPath}" --genome "${outputDirectory}"assembly/all.fasta --pacbio "${outputDirectory}"assembly/all.bam --threads $threads --outdir "${outputDirectory}"assembly
#mv "${outputDirectory}"assembly/pilon.fasta "${outputDirectory}"assembly/all.fasta
#rm "${outputDirectory}"assembly/all.fasta.gz
#gzip "${outputDirectory}"assembly/all.fasta
#rm -r "${outputDirectory}"readFilteringAndTrimming
#echo ""
#echo "Assembly Error Correction Complete"
#echo ""

### Assembly Error Correction QC ###
## Error Correction Align to Reference Genome
#gunzip "${outputDirectory}"assembly/all.fasta
#assembly-stats "${outputDirectory}"assembly/all.fasta > "${outputDirectory}"qualityControl/postAssemblyErrorCorrection/assemblyStats.txt
## Dnadiff (Reference Genome Alignment) Genome Assembly Quality Control
#cd "${outputDirectory}"qualityControl/postAssemblyErrorCorrection
#dnadiff -p dnadiff "${genomeReferencePath}" "${outputDirectory}"assembly/all.fasta
#gzip -q "${outputDirectory}"assembly/all.fasta
#echo ""
#echo "Post-Assembly Error Correction Quality Control Complete"
#echo ""

### Alignment ###
minimap2 --MD -ax map-ont -K 800000M -t $threads "${genomeReferencePath}" "${outputDirectory}"assembly/all.fasta.gz > "${outputDirectory}"alignment/all.sam
gzip -q "${outputDirectory}"alignment/all.sam
# Convert to bam file
samtools view --threads $threads -bS "${outputDirectory}"alignment/all.sam.gz > "${outputDirectory}"alignment/all.bam
# Sort the bam file
samtools sort --threads $threads -o "${outputDirectory}"alignment/all.bam "${outputDirectory}"alignment/all.bam
# create an index file
samtools index -@ $threads "${outputDirectory}"alignment/all.bam
rm -r "${outputDirectory}"assembly
echo ""
echo "Alignment Complete"
echo ""

### Variant Analysis ###
## Structural Variant Calling
sniffles -t $threads -i "${outputDirectory}"alignment/all.bam -v "${outputDirectory}"variantCalling/structuralVariants/all.vcf --allow-overwrite
source deactivate miniStructuralVariantAnalysis 
## Transposon Calling
source activate tldr
cd "${outputDirectory}"variantCalling/transposons
tldr -b "${outputDirectory}"alignment/all.bam -e /users/k20064105/tldr/ref/teref.ont.human.fa -r "${genomeReferencePath}" --color_consensus
echo ""
echo "Variant Analysis Complete"
echo ""