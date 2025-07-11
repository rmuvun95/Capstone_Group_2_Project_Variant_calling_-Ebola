#!/usr/bin/bash

# ==============================================================================
# Bioinformatics Pipeline: QC, Assembly
# Author: Group Two(Raissa, Renatha, John and Liya)
# Date: 2025-07-10
# Version: 1.0
# Description:
#   This script orchestrates a complete bioinformatics workflow:
#     1. Fastp trimming and filtering.
#     2. FastQC on trimmed reads.
#     3. Aligning using BWA and Samtools.
#     4. Variant calling using BCF Tools.
#   It automatically detects samples from filenames: *_R1.fastq.gz and *_R2.fastq.gz
# ==============================================================================

#Set pipeline error

set -euo pipefail

#Set colors

GREEN='\033[1;32m'
RED='\033[1;31m'
YELLOW='\033[1;33m'
BLUE='\033[1;36m'
BOLD='\033[1m'
NC='\033[0m'

## Print start of the process and date time
echo "==========================" 
START_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e  "${GREEN}$START_DATE_TIME : RUNNING Variant Calling Analysis for Viral Pathogens${NC}" 
 echo "==========================" 

#Setting input and output directories
INPUT_READS_DIR=""
OUTPUT_BASE_DIR="./results/"
REF_SEQ=""
THREADS=""

usage() {
  printf "${BLUE}Bioinformatics Pipeline: QC → Alignment → Variant Calling${NC}\n"
  printf "${BLUE}------------------------------------------------------${NC}\n"
  printf "${BLUE}This script performs the following for each sample:${NC}\n"
  printf "${BLUE}  1. FastQC on raw reads${NC}\n"
  printf "${BLUE}  2. Trimming and filtering with fastp${NC}\n"
  printf "${BLUE}  3. FastQC and MultiQC summary report generation on trimmed reads${NC}\n"
  printf "${BLUE}  4. Alignment with samtools${NC}\n"
  printf "${BLUE}  5. Variant Calling using Bcftools ${NC}\n"
  printf "${BLUE}  6. Compressing and Indexing VCF${NC}\n"
  printf "${BLUE}Samples are auto-detected from *_1.fastq.gz in the input folder.${NC}\n"
  printf "\n"
  printf "${BLUE}Usage: %s -i <input_reads_folder> -o <output_folder> -r <reference_folder> [-t <threads>]\n" "$(basename "$0")"
  printf "  -i   Folder containing raw FASTQ files (e.g., raw_reads/)\n"
  printf "  -o   Output folder for all results (e.g., results/)\n"
  printf "  -r   reference folder to align to (e.g., ref_folder/)\n"
  printf "  -t   Number of threads for parallel processes [default: 4]\n"
  printf "  -h   Show this help and exit${NC}\n"
  exit 1
}
# Parsing input parameters
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) INPUT_READS_DIR="$2"; shift 2 ;;
    -o|--outdir) OUTPUT_BASE_DIR="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    -r|--ref) REF_SEQ="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) printf "${RED}❌ ERROR: Unknown option: %s${NC}\n" "$1"; usage ;;
  esac
done

# Check required inputs
if [[ -z "$INPUT_READS_DIR"  || -z "$OUTPUT_BASE_DIR" || -z "$REF_SEQ" || -z "$THREADS" ]]; then
  printf "${RED}❌ ERROR: -i (input reads folder),  -o (output_folder), -r (reference_folder) , -t (number of threads) are required.${NC}\n"
  usage
fi

# Validate input directory
if [[ ! -d "$INPUT_READS_DIR" ]]; then
  printf "${RED}❌ ERROR: Input reads folder '%s' not found or not accessible inside the container.${NC}\n" "$INPUT_READS_DIR"
  printf "${YELLOW}   ➤ Ensure you used '--bind /local/path/to/data:/data' and your reads are in '/local/path/to/data/raw_reads_folder'.${NC}\n"
  exit 1
fi

# Validate reference directory
if [[ ! -d "$REF_SEQ" ]]; then
  printf "${RED}❌ ERROR: Reference sequence folder '%s' not found or not accessible inside the container.${NC}\n" "$REF_SEQ"
  printf "${YELLOW}   ➤ Ensure you used '--bind /local/path/to/data:/data' and your refrences are in '/local/path/to/data/reference_folder'.${NC}\n"
  exit 1
fi


# Check tools
for tool in fastqc fastp multiqc bwa samtools bcftools seqkit tabix; do
  if ! command -v "$tool" &> /dev/null; then
    printf "${RED}❌ ERROR: Required tool '%s' not found in PATH.${NC}\n" "$tool"
    printf "${YELLOW}   ➤ Install the required tool before proceeding.${NC}\n"
    exit 1
  fi
done

# Prepare output subfolders
RAW_FASTQC_DIR="$OUTPUT_BASE_DIR/qc"
TRIMMED_READS_DIR="$OUTPUT_BASE_DIR/trimmed"
TRIMMED_FASTQC_DIR="$OUTPUT_BASE_DIR/qc/fastqc_post_trim"
ALIGNMENT_DIR_SAM="$OUTPUT_BASE_DIR/sam"
ALIGNMENT_DIR="$OUTPUT_BASE_DIR/bam"
VARIANT_CALLING_DIR="$OUTPUT_BASE_DIR/vcf"
LOGS_DIR="./logs/"
PIPELINE_REPORTS="$LOGS_DIR/reports.txt"
mkdir -p "$RAW_FASTQC_DIR" "$TRIMMED_READS_DIR" "$TRIMMED_FASTQC_DIR" \
         "$ALIGNMENT_DIR_SAM" "$ALIGNMENT_DIR" "$VARIANT_CALLING_DIR" "$LOGS_DIR"
> "$PIPELINE_REPORTS"

# Prompting message for starting the analysis

printf "${BOLD}${BLUE}🔧 Bioinformatics Pipeline Started at [%s]${NC}\n" "$(date)"
printf "${BLUE}📂 Input Reads: %s${NC}\n" "$INPUT_READS_DIR"
printf "${BLUE}📦 Output Base Directory: %s${NC}\n" "$OUTPUT_BASE_DIR"
printf "${BLUE}📦 Reference Directory: %s${NC}\n" "$REF_SEQ"
printf "${BLUE}⚙️ Threads: %s${NC}\n" "$THREADS"
printf "\n${BLUE}--- Detecting and processing samples from filenames ---${NC}\n"

# Indexing the reference 

find "$REF_SEQ" -name '*.fa' | sort | while read -r FA_PATH; do
  FA_FILE=$(basename "$FA_PATH")
  echo "Indexing $FA_FILE"
  bwa index "$FA_PATH"
done


find "$INPUT_READS_DIR" -name '*_1.fastq.gz' | sort | while read -r R1_PATH; do
  R1_FILE=$(basename "$R1_PATH")
  SAMPLE_ID="${R1_FILE%%_1.fastq.gz}"
  echo "$SAMPLE_ID" "$R1_FILE"
  R2_FILE="${SAMPLE_ID}_2.fastq.gz"

  printf "\n${BOLD}----------------------------------------------------${NC}\n"
  printf "${BLUE}Processing Sample: ${BOLD}%s${NC}\n" "$SAMPLE_ID"
  printf "${BLUE}  R1: %s${NC}\n" "$R1_FILE"
  printf "${BLUE}  R2: %s${NC}\n" "$R2_FILE"

#Specifying output
INPUT_R1_PATH="${INPUT_READS_DIR}/${R1_FILE}"
  INPUT_R2_PATH="${INPUT_READS_DIR}/${R2_FILE}"
  TRIMMED_R1_PATH="${TRIMMED_READS_DIR}/${SAMPLE_ID}_1_trimmed.fastq.gz"
  TRIMMED_R2_PATH="${TRIMMED_READS_DIR}/${SAMPLE_ID}_2_trimmed.fastq.gz"
  FASTP_JSON_REPORT="${TRIMMED_READS_DIR}/${SAMPLE_ID}.fastp.json"
  FASTP_HTML_REPORT="${TRIMMED_READS_DIR}/${SAMPLE_ID}.fastp.html"
  ALIGNMENT_DIR_SAM="${ALIGNMENT_DIR_SAM}/${SAMPLE_ID}.sam"
  ALIGNMENT_DIR="${ALIGNMENT_DIR}/${SAMPLE_ID}.bam"
  VARIANT_CALLING_DIR="${VARIANT_CALLING_DIR}/${SAMPLE_ID}.vcf.gz"

#Step 1: Quality Control with fastqc

printf "${BLUE}▶ Step 1: Running FastQC on raw reads ${SAMPLE_ID}...${NC}\n"
  if [[ -f "$INPUT_R1_PATH" && -f "$INPUT_R2_PATH" ]]; then
    fastqc -t "$THREADS" -o "$RAW_FASTQC_DIR" "$INPUT_R1_PATH" "$INPUT_R2_PATH" &> "$LOGS_DIR/${SAMPLE_ID}_fastqc_raw.log"
    printf "${GREEN}  FastQC (raw) completed.${NC}\n"
  fi

 #Step 2: Trimming and filtering with fastp

# We are trimming quality < 30
# We are trimming overrepresented and poly G's with the flags -p nad -g respectfully 
# We are trimming reads less than length of 90 with the flag -l
# We are trimming 5 reads from the tail from read 1 and read 2 with the flags -t and -T respectfully
# We are calling the threads we specified earlier with the flag -w

  printf "${BLUE}▶ Step 2: Running fastp for trimming ${SAMPLE_ID}...${NC}\n"
  fastp -i "$INPUT_R1_PATH" -I "$INPUT_R2_PATH" \
        -o "$TRIMMED_R1_PATH" -O "$TRIMMED_R2_PATH" \
        --json "$FASTP_JSON_REPORT" --html "$FASTP_HTML_REPORT" \
        -q 30 -p -g -l 90 -t 5 -T 5 -w "$THREADS" \ 
        &> "$LOGS_DIR/${SAMPLE_ID}_fastp.log"
  printf "${GREEN}  fastp trimming completed.${NC}\n"

# Running fastqc on the trimmed reads

  printf "${BLUE}▶ Step 3: Running FastQC on trimmed reads ${SAMPLE_ID}...${NC}\n"
   if [[ -f "$TRIMMED_R1_PATH" && -f "$TRIMMED_R2_PATH" ]]; then
  fastqc -t "$THREADS" -o "$TRIMMED_FASTQC_DIR" "$TRIMMED_R1_PATH" "$TRIMMED_R2_PATH" &> "$LOGS_DIR/${SAMPLE_ID}_fastqc_trimmed.log"
  printf "${GREEN}  FastQC (trimmed) completed.${NC}\n"
  fi

# Running alignment


printf "${BLUE}▶ Step 4: Running alignment with BWA for ${SAMPLE_ID}...${NC}\n"

BAM_FILE="${ALIGNMENT_DIR}/${SAMPLE_ID}.bam"
echo "$BAM_FILE"
# Align and convert to BAM
bwa mem -t "$THREADS" "$REF_SEQ" "$TRIMMED_R1_PATH" "$TRIMMED_R2_PATH" | \
  samtools view -Sb > "$BAM_FILE"

# Index BAM
samtools index "$BAM_FILE"

printf "${GREEN}  Alignment and BAM indexing completed for ${SAMPLE_ID}.${NC}\n"




# Variant calling with bcftools

printf "${BLUE}▶ Step 5: Running variant calling with bcftools ${SAMPLE_ID}...${NC}\n"
bcftools mpileup -f "FA_PATH" "$ALIGNMENT_DIR" | bcftools call --ploidy 1 -mv -Oz -o "$VARIANT_CALLING_DIR"

# Indexing vcf file

tabix -p vcf "$VARIANT_CALLING_DIR"

done


echo "==========================" 
END_DATE_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo -e  "${GREEN}$END_DATE_TIME :  Pipeline completed successfully for all samples!${NC}" 
echo "==========================" 