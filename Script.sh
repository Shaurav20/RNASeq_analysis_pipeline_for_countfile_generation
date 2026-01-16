#!/bin/bash

# ==============================================================================
# RNA-Seq Analysis Pipeline for HCT116 Colorectal Carcinoma Samples
# Project: SRP603504 - Comparison of untreated vs PMA/Ionomycin treated cells
# Samples: 
#   SRR34712398 - Untreated HCT116 (Cancer sample)
#   SRR34712393 - PMA/Ionomycin treated HCT116
# Analysis steps: SRA download → QC → Trimming → Alignment → Quantification
# Reference genome: GRCh38 (chromosome 1 only for testing)
# ==============================================================================

# ==============================================================================
# SECTION 0: SETUP AND INSTALLATION
# ==============================================================================

# Note: Windows PowerShell commands for initial setup (run in PowerShell as Administrator):

# Install Windows Subsystem for Linux (WSL) and Ubuntu
# wsl --install -d Ubuntu

# After Ubuntu installation, set username and password when prompted

# Once Ubuntu is installed, you can run this bash script from WSL/Ubuntu terminal

# Update and upgrade Ubuntu system
echo "Updating Ubuntu system packages..."
sudo apt update && sudo apt upgrade -y

# Install essential bioinformatics tools via apt
echo "Installing basic bioinformatics tools..."
sudo apt install -y wget curl build-essential

# Install Miniconda for Python package management
echo "Installing Miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
rm miniconda.sh

# Initialize conda for bash shell
echo "Initializing conda..."
source $HOME/miniconda/etc/profile.d/conda.sh

# Add conda to PATH permanently
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# Configure conda channels
echo "Configuring conda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Create dedicated RNA-seq environment
echo "Creating rnaseq_clean environment..."
conda create -n rnaseq_clean python=3.10 -y

# Activate the environment
conda activate rnaseq_clean

# Install bioinformatics tools in the environment
echo "Installing bioinformatics tools..."
conda install -y -c bioconda \
    hisat2 \
    subread \
    fastqc \
    sra-tools \
    trimmomatic

# Verify installations
echo "Verifying installations..."
hisat2 --version
featureCounts -v
fastqc --version

# ==============================================================================
# SECTION 1: DATA DOWNLOAD AND PREPARATION
# ==============================================================================

# Create project directory structure
echo "Creating project directory structure..."
mkdir -p /mnt/d/Projects/New_folder_3/{Biosample_1,Biosample_2,Results,Reference}

# Download reference genome (chromosome 1 only for testing)
echo "Downloading reference genome (chromosome 1)..."
cd /mnt/d/Projects/New_folder_3/Reference
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz

# Download GTF annotation file
echo "Downloading GTF annotation file..."
wget https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz
gunzip Homo_sapiens.GRCh38.110.gtf.gz

# Note: For full genome analysis, download the complete reference:
# wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

# ==============================================================================
# SECTION 2: BIOSAMPLE 1 PROCESSING (Untreated HCT116 - SRR34712398)
# ==============================================================================

# Note: The prefetch step should be run in Windows PowerShell first:
# prefetch SRR34712398

echo "Processing Biosample 1 (SRR34712398)..."
cd /mnt/d/Projects/New_folder_3/Biosample_1

# Extract FASTQ from prefetched SRA file
# Note: Assuming SRA file is already prefetched in Windows directory
echo "Extracting FASTQ files from prefetched SRA..."
fasterq-dump SRR34712398 --split-files --threads 4 --progress --outdir .

# Quality control with FastQC
echo "Running FastQC on raw reads..."
fastqc SRR34712398_1.fastq SRR34712398_2.fastq

# Trimming with Trimmomatic (using conda-installed version)
echo "Trimming reads with Trimmomatic..."
# Note: Create polyG.fa adapter file for overrepresented sequences
echo ">polyG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" > polyG.fa

trimmomatic PE -phred33 \
SRR34712398_1.fastq SRR34712398_2.fastq \
SRR34712398_1_paired.fastq SRR34712398_1_unpaired.fastq \
SRR34712398_2_paired.fastq SRR34712398_2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
ILLUMINACLIP:polyG.fa:2:30:10 \
HEADCROP:12 SLIDINGWINDOW:4:20 TRAILING:10 MINLEN:36

# Quality control on trimmed reads
echo "Running FastQC on trimmed reads..."
fastqc SRR34712398_1_paired.fastq SRR34712398_2_paired.fastq

# Build HISAT2 index for chromosome 1
echo "Building HISAT2 index for chromosome 1..."
hisat2-build /mnt/d/Projects/New_folder_3/Reference/Homo_sapiens.GRCh38.dna.chromosome.1.fa GRCh38_chr1

# Align trimmed reads to chromosome 1 index
echo "Aligning reads with HISAT2..."
hisat2 -q --phred33 -p 8 --dta -x GRCh38_chr1 \
-1 SRR34712398_1_paired.fastq \
-2 SRR34712398_2_paired.fastq \
-S test_aligned.sam \
--summary-file test_alignment_summary.txt

# Prepare chromosome 1 GTF file
echo "Preparing chromosome 1 GTF file..."
cd /mnt/d/Projects/New_folder_3/Reference
grep "^1\s" Homo_sapiens.GRCh38.110.gtf > chr1_body.gtf
grep "^#" Homo_sapiens.GRCh38.110.gtf > chr1_header.gtf
cat chr1_header.gtf chr1_body.gtf > Homo_sapiens.GRCh38.110.chr1.gtf

# Count reads at gene level
echo "Counting reads at gene level..."
cd /mnt/d/Projects/New_folder_3/Biosample_1
featureCounts -T 8 -p -t exon -g gene_id \
-a ../Reference/Homo_sapiens.GRCh38.110.chr1.gtf \
-o test_gene_counts.txt test_aligned.sam

# Prepare clean gene count file for R analysis
tail -n +3 test_gene_counts.txt | cut -f1,7 > cancer_gene_counts.txt

# ==============================================================================
# SECTION 3: BIOSAMPLE 2 PROCESSING (PMA/Ionomycin treated HCT116 - SRR34712393)
# ==============================================================================

# Note: The prefetch step should be run in Windows PowerShell first:
# prefetch SRR34712393

echo "Processing Biosample 2 (SRR34712393)..."
cd /mnt/d/Projects/New_folder_3/Biosample_2

# Extract FASTQ from prefetched SRA file
echo "Extracting FASTQ files from prefetched SRA..."
fasterq-dump SRR34712393 --split-files --threads 4 --progress --outdir .

# Quality control
fastqc SRR34712393_1.fastq SRR34712393_2.fastq

# Trimming
echo "Trimming reads with Trimmomatic..."
# Note: Create polyG1.fa adapter file for overrepresented sequences
echo ">polyG1
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" > polyG1.fa

trimmomatic PE -phred33 \
SRR34712393_1.fastq SRR34712393_2.fastq \
SRR34712393_1_paired.fastq SRR34712393_1_unpaired.fastq \
SRR34712393_2_paired.fastq SRR34712393_2_unpaired.fastq \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
ILLUMINACLIP:polyG1.fa:2:30:10 \
HEADCROP:12 SLIDINGWINDOW:4:20 TRAILING:10 MINLEN:36

# Quality control on trimmed reads
fastqc SRR34712393_1_paired.fastq SRR34712393_2_paired.fastq

# Align using the same index
hisat2 -q --phred33 -p 8 --dta -x ../Biosample_1/GRCh38_chr1 \
-1 SRR34712393_1_paired.fastq \
-2 SRR34712393_2_paired.fastq \
-S test_aligned.sam \
--summary-file test_alignment_summary.txt

# Count reads at gene level
featureCounts -T 8 -p -t exon -g gene_id \
-a ../Reference/Homo_sapiens.GRCh38.110.chr1.gtf \
-o test_gene_counts.txt test_aligned.sam

# Prepare clean gene count file for R analysis
tail -n +3 test_gene_counts.txt | cut -f1,7 > normal_gene_counts.txt

# ==============================================================================
# SECTION 4: ORGANIZE RESULTS
# ==============================================================================

echo "Organizing results in Results folder..."
cd /mnt/d/Projects/New_folder_3

# Create Results directory if not exists
mkdir -p Results

# Copy essential result files to Results folder
echo "Copying essential result files..."

# Biosample 1 results
cp Biosample_1/cancer_gene_counts.txt Results/
cp Biosample_1/test_alignment_summary.txt Results/cancer_test_alignment_summary.txt

# Biosample 2 results
cp Biosample_2/normal_gene_counts.txt Results/
cp Biosample_2/test_alignment_summary.txt Results/normal_test_alignment_summary.txt

# ==============================================================================
# SECTION 5: CLEANUP AND FINAL REPORT
# ==============================================================================

# Clean up large intermediate files (optional)
echo "Cleaning up intermediate files..."
rm -f Biosample_1/*.sam Biosample_2/*.sam
rm -f Biosample_1/*unpaired.fastq Biosample_2/*unpaired.fastq

# Print completion report
echo "==========================================="
echo "RNA-Seq Pipeline Completed Successfully!"
echo "==========================================="
echo "Project Directory: /mnt/d/Projects/New_folder_3/"
echo ""
echo "IMPORTANT: Prefetch steps should be completed in Windows PowerShell first:"
echo "1. Open PowerShell as Administrator"
echo "2. Run: prefetch SRR34712398"
echo "3. Run: prefetch SRR34712393"
echo ""
echo "Essential Results Files (in Results/ folder):"
echo "1. Count files for differential expression analysis:"
echo "   - cancer_gene_counts.txt"
echo "   - normal_gene_counts.txt"
echo ""
echo "2. Alignment statistics:"
echo "   - cancer_test_alignment_summary.txt"
echo "   - normal_test_alignment_summary.txt"
echo ""
echo "To use in R for differential expression analysis:"
echo "1. Load count files: cancer_gene_counts.txt, normal_gene_counts.txt"
echo "2. Use DESeq2 or edgeR for analysis"
echo "==========================================="

# Deactivate conda environment
conda deactivate

