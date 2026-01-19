# Variants_pipeline_ONT

Quick and dirty pipeline to extract **SNVs, SVs and CNAs** from **ONT tumor-only BAM/UBAM data**, mainly used on **GBM cell lines**.

The final output is a set of **summary plots** showing SNV / SV / CNV signals on a predefined list of genes used to profile GBM samples.

This is an **internal lab tool**, meant for fast exploratory analyses, not a polished or fully optimized workflow.

---

## What this pipeline does

Starting from an ONT BAM:

- SNV/INDEL calling with **ClairS-TO**
- Variant phasing + haplotagging with **Longphase**
- Annotation with **VEP**
- SV and CNA calling with **SAVANA (tumor-only mode)**
- Intersection with a list of genes of interest
- Automatic plotting of per-gene summaries

Everything is wrapped in a single SLURM bash script, plus a Python script for plotting.

---

## Where it was tested

The pipeline was developed and run on **LEONARDO (CINECA)**.

A non-trivial part was fitting the whole workflow inside the allowed **wall time and resource budget**, without recomputing expensive steps unnecessarily.

---

## Notes / gotchas (read this before running)

- **Longphase**
  Installing it directly was annoying, so it is run via **Singularity** using a `.sif` image.

- **SAVANA tumor-only**
  Tumor-only mode is a bit weird and the documentation is not super clear.

  Current behaviour:
  - Using a **phased / haplotagged BAM** → CNA works, SV VCF is empty
  - Using an **unphased BAM** → SVs are called correctly, CNAs are not

  For this reason:
  - CNA calling is done on the phased BAM
  - SV calling is done on the unphased BAM

  This seems to be the expected workaround for now. The reason is unclear. I will probably open an issue to their GitHub repo

- **Parallelization**
  Could definitely be improved, but for now it does the job.  
  Since tools change quickly, further optimization is postponed.

-**Instalaltions**
 As said, longphase-to was installed through (https://github.com/CCU-Bioinformatics-Lab/longphase-to?tab=readme-ov-file#installation) Docker the .sif image was obtained
 SAVANA (https://github.com/cortes-ciriano-lab/savana) 
 and ClairS-TO (https://github.com/HKU-BAL/ClairS-TO?tab=readme-ov-file#installation) were installed through conda (option 3 on ClairS's repo).
 Anyway, yaml files of the used environments are provided in this repo

---

## Usage

- Adjust all paths at the top of `scripts/run_pipeline.sh`
- The pipeline assumes:
  - SLURM
  - Conda environments already available
  - Singularity for Longphase
- Final plots are produced automatically at the end of the bash script via a Python script.

See `scripts/run_pipeline.sh` for the full workflow.

---

## Future

If/when this becomes more stable:
- refactor into Snakemake or Nextflow
- clean up environments and dependencies

For now, this is meant to be **fast, usable, and easy to tweak**.



