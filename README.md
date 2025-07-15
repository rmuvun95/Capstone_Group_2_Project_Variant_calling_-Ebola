# ACDC Bioinformatics Training
## Capstone Project: Variant_calling_Ebola
### 🧬 Bioinformatics Pipeline steps: QC → Alignment → Variant Calling

### 🎯 Objective

This guide walks you through building and running a Singularity container for a variant calling pipeline on raw sequencing data.

You will:
- Build a custom container with essential tools (FastQC, fastp, BWA, bcftools, MultiQC)
- Use Miniforge and `mamba` to manage environments via `ebovar.yml`
- Execute the full pipeline using a wrapper script inside the container

---

### 1. 🔧 Tools Installed

The pipeline includes:

| Tool        | Purpose                          |
|-------------|----------------------------------|
| **FastQC**  | Raw read quality control         |
| **fastp**   | Trimming and filtering           |
| **BWA**     | Sequence alignment               |
| **Bcftools**| Variant calling                  |
| **MultiQC** | Aggregated reporting             |

All dependencies are managed through a Conda environment defined in `ebovar.yml`.

---

### 2. ⚙️ Container Definition Overview

Your Singularity definition file (`ebovar.def`) includes:

- Ubuntu 22.04 as the base OS
- System dependency installation via APT
- Miniforge + Mamba installation for Conda environment setup
- Environment activation for downstream use
- Pipeline script placement (`eboVar.sh`)

---

### 3. 📂 File Structure

Ensure the following files are present before building:

```
project/
├── ebovar.def          # Singularity definition file
├── ebovar.yml          # Conda environment definition
└── eboVar.sh           # Executable pipeline script
```
---

### 4. 🏗️ Building the Container

Run this from the directory containing `ebovar.def`:

```bash
sudo apptainer build ebovar.sif ebovar.def
```

---

### 5. ▶ Running the Pipeline

Run the container with your raw reads, output folder, and reference genome mounted:

```bash
apptainer run --bind $(pwd):/data ebovar.sif \
  -i /data/raw \
  -o /data/results \
  -r /data/ref_folder \
  -t <threads>
```

---

### 6. 🧾 Command Line Options (eboVar.sh)

| Flag        | Description                          |
|-------------|--------------------------------------|
| `-i`        | Input folder with raw reads          |
| `-o`        | Output directory                     |
| `-r`        | Reference genome folder              |
| `-t`        | Number of threads                    |

---

### 7. 📤 Output Structure

```
results/
├── fastqc/         # Quality control reports
├── trimmed/        # fastp output
├── alignment/      # BWA-aligned BAM files
├── variants/       # VCF output from bcftools
└── multiqc/        # Aggregated QC report
```

---

### 8. 🧪 Notes

- Add threads depending on your capacity
- Ensure your fastq format is {sample_name}_1.fastq.gz and {sample_name}_2.fastq.gz for read1 and read2 respectively

---

### 👩‍💻 Author Group

**Group Two**
- Raissa  
- Renatha  
- John  
- Liya



