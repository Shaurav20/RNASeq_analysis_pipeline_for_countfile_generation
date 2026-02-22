# RNASeq_analysis_pipeline_for_countfile_generation

## Project Overview

This project implements a complete RNA-Seq analysis pipeline for comparing untreated and PMA/Ionomycin-treated HCT116 colorectal carcinoma cell lines. The analysis processes two samples from the SRA study SRP603504 using a hybrid Windows/WSL workflow.

## Project Structure

```
/mnt/d/Projects/New_folder_3/
├── README.md                   # Project documentation
├── Script.docx                 # Main analysis pipeline script
├── Reference/                  # Genome reference files
│   ├── Homo_sapiens.GRCh38.dna.chromosome.1.fa
│   ├── Homo_sapiens.GRCh38.110.gtf
│   └── Homo_sapiens.GRCh38.110.chr1.gtf
├── Biosample_1/                # Untreated HCT116 (SRR34712398)
│   ├── SRR34712398_1.fastq
│   ├── SRR34712398_2.fastq
│   ├── SRR34712398_1_paired.fastq
│   ├── SRR34712398_2_paired.fastq
│   ├── SRR34712398_1_unpaired.fastq
│   ├── SRR34712398_2_unpaired.fastq
│   ├── FastQC_reports_raw/     # FastQC reports for raw reads
│   ├── FastQC_reports_trimmed/ # FastQC reports for trimmed reads
│   ├── GRCh38_chr1.*           # HISAT2 index files
│   ├── test_aligned.sam
│   ├── test_alignment_summary.txt
│   ├── test_gene_counts.txt
│   └── cancer_gene_counts.txt
├── Biosample_2/                # PMA/Ionomycin treated HCT116 (SRR34712393)
│   ├── SRR34712393_1.fastq
│   ├── SRR34712393_2.fastq
│   ├── SRR34712393_1_paired.fastq
│   ├── SRR34712393_2_paired.fastq
│   ├── SRR34712393_1_unpaired.fastq
│   ├── SRR34712393_2_unpaired.fastq
│   ├── FastQC_reports_raw/     # FastQC reports for raw reads
│   ├── FastQC_reports_trimmed/ # FastQC reports for trimmed reads
│   ├── test_aligned.sam
│   ├── test_alignment_summary.txt
│   ├── test_gene_counts.txt
│   └── normal_gene_counts.txt
├── Results/                    # Final analysis results
│   ├── cancer_gene_counts.txt
│   ├── normal_gene_counts.txt
│   ├── cancer_test_aligned.sam
│   ├── normal_test_aligned.sam
│   ├── cancer_test_alignment_summary.txt
│   ├── normal_test_alignment_summary.txt
│   └── metadata.csv
```

## Pipeline Steps

**Data Acquisition:** Download SRA files using prefetch and fasterq-dump

**Quality Control:** FastQC analysis of raw and trimmed reads

**Read Trimming:** Adapter and quality trimming with Trimmomatic

**Reference Preparation:** HISAT2 index building for GRCh38 chromosome 1

**Alignment:** Read alignment with HISAT2

**Quantification:** Gene counting with featureCounts

**Analysis Preparation:** Formatting count files for R/DESeq2

## Installation & Setup

### Prerequisites

Windows 10/11 with PowerShell administrator access

Minimum 50GB free disk space

8GB RAM minimum (16GB recommended)

### Setup Instructions

**1. Windows Setup (PowerShell as Administrator)**

Within windows powershell - Enable WSL and install Ubuntu

```
wsl --install -d Ubuntu 

# Set up SRA Toolkit (Windows installation)
# Download from: https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
```

**2. Data Download (PowerShell)**

Download SRA files (run in PowerShell)
```
prefetch SRR34712398
prefetch SRR34712393

# Verify files are downloaded
# Files will be in: C:\Users\[Username]\ncbi\public\sra\
```

**3. Run Analysis Pipeline (WSL/Ubuntu)**

Execute the complete pipeline script in WSL/Ubuntu terminal.

## Problems Encountered and Solutions

### Problem 1: Windows Bioinformatics Tool Compatibility

**Issue:** Initial attempts to run bioinformatics tools (bowtie2, samtools, HISAT2) directly in Windows PowerShell failed due to:

Missing dependencies

Incompatible builds for Windows architecture

Environment variable conflicts

Slow performance with large genome indexing

**Solution:** Implemented a hybrid Windows/WSL approach:

Use Windows PowerShell for data download (prefetch/fasterq-dump)

Use WSL/Ubuntu for computational analysis

Leverage Linux-native bioinformatics tools with better performance and compatibility

### Problem 2: Genome Index Building Failures

**Issue:** Multiple attempts to build genome indices failed:

bowtie2-build on full GRCh38 genome stalled indefinitely

Memory exhaustion during index building

Incorrect file permissions in Windows-mounted directories

Timeouts due to large genome size

**Solution:**

Chromosome 1 Only Strategy: Initially built index for chromosome 1 only to test pipeline

WSL Optimization: Moved reference files to WSL-native storage (/home/ instead of /mnt/d/)

Resource Management: Added --offrate and --bmax parameters for bowtie2

Alternative Aligner: Switched to HISAT2 which has better memory management

### Problem 3: Conda Environment Issues

**Issue:** Multiple conda environment problems:

Channel conflicts between bioconda and conda-forge

Package version incompatibilities

Broken dependencies after installation

SAMtools installation failures

**Solution:**

Clean Environment: Created fresh rnaseq_clean environment

Strict Channel Priority: Set conda config --set channel_priority strict

Sequential Installation: Installed packages in dependency order

Alternative Installation Methods: Used mamba for faster dependency resolution

### Problem 4: WSL/Windows File System Performance

**Issue:** Poor performance when accessing Windows files from WSL:

Very slow file I/O on /mnt/d/ mounted drives

Permission issues with Windows-created files

File locking conflicts between Windows and WSL

**Solution:**

Work in WSL Native Storage: Perform computation in /home/user/ directories

Copy Before Processing: Move files from /mnt/d/ to WSL-native storage before analysis

Use WSL-native Tools: Install bioinformatics tools in WSL, not Windows

Avoid Concurrent Access: Don't access same files simultaneously from Windows and WSL

### Problem 5: Quality Control Challenges

**Issue:** Initial FastQC reports showed several quality issues:

Overrepresented poly-G sequences (likely technical artifacts)

Adapter contamination

Poor quality bases at read ends

**Solution:**

Custom Trimming: Added custom poly-G sequence to Trimmomatic adapters file

Head Crop: Implemented HEADCROP:12 to remove poor quality bases

Multi-step Trimming: Used dual ILLUMINACLIP for both adapters and poly-G sequences

Post-trimming QC: Run FastQC after trimming to verify improvement

### Problem 6: Alignment Rate Issues

**Issue:** Low alignment rates (~14%) when using chromosome 1 only:

Most reads couldn't align to single chromosome

Expected for whole-transcriptome data aligned to partial genome

**Solution:**

Pipeline Validation: Confirmed pipeline works despite low alignment rate

Documentation: Clearly noted this is expected behavior for testing

Production Readiness: Pipeline validated and ready for full genome analysis

Alternative Explanation: Low alignment rate is correct for chromosome 1-only alignment

## Key Technical Decisions

**1. Hybrid Windows/WSL Architecture**

Decision: Use Windows for data download, WSL for computation
Rationale:

SRA Toolkit works well in Windows

Bioinformatics tools perform better in Linux

Best of both operating systems

**2. Chromosome 1 Testing Strategy**

Decision: Test pipeline on chromosome 1 before full genome
Rationale:

Faster iteration during development

Lower resource requirements

Validates pipeline before committing to full analysis

**3. HISAT2 over bowtie2**

Decision: Use HISAT2 as primary aligner
Rationale:

Better memory management

Built for RNA-Seq with splice-aware alignment

More reliable index building

**4. featureCounts for Quantification**

Decision: Use featureCounts from subread package
Rationale:

Fast and efficient

Works directly with SAM files (no BAM conversion needed)

Well-documented and widely used

## Performance Metrics

Read Processing: ~56 million read pairs per sample

Trimming Efficiency: ~95% read survival rate

Alignment Rate: ~14% (chromosome 1 only - expected)

## Next Steps for Full Analysis

**1. Download Complete GRCh38 Genome Reference**

Obtain the full genome FASTA file from Ensembl or UCSC

Use complete annotation file (GTF) for comprehensive gene counting

**2. Build Full HISAT2 Genome Index**

Create index for entire genome instead of chromosome 1 only

This will significantly improve alignment rates (>70% expected)

**3. Run Complete Alignment**

Align all reads to full genome reference

Perform quality assessment on complete alignment

**4. Perform Differential Expression Analysis with DESeq2**

Note: For statistically meaningful DESeq2 analysis:

-Minimum 3 biological replicates per condition are required

-Current analysis uses only 2 samples (1 per condition) for pipeline testing only

-In production, collect count files from multiple samples per condition

-DESeq2 requires replicates for variance estimation and reliable statistics

## Conclusion

This pipeline successfully overcomes multiple technical challenges through a hybrid Windows/WSL approach, careful resource management, and iterative testing. The validated pipeline is now ready for production-scale analysis on complete genomes with appropriate sample replication.

## License
This project is licensed under the MIT License. You are free to use, modify, and distribute this code with appropriate attribution.

## Author
Shaurav Bhattacharyya

- GitHub: @Shaurav20 (https://github.com/Shaurav20)
- LinkedIn: Shaurav Bhattacharyya (https://www.linkedin.com/in/shaurav-bhattacharyya-a347b5156/)
- Email: shaurav.feb@gmail.com

## Acknowledgements
SRA Toolkit – NCBI

FastQC – Babraham Bioinformatics

Trimmomatic – Bolger, Lohse, Usadel (2014)

HISAT2 – Kim, Langmead, Salzberg (2015)

featureCounts – Liao, Smyth, Shi (2014)

SAMtools – Danecek et al. (2021)

GRCh38 – Genome Reference Consortium

Ensembl – Release 110
