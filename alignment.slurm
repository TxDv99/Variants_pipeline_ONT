#!/bin/bash
#SBATCH --job-name=ALGN
#SBATCH --output=algn.%A_%a.out
#SBATCH --error=agln.%A_%a.err
#SBATCH --array=1-112
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=t.ducci@student.unisi.it
#SBATCH --account=IscrC_KREIOS_0
#SBATCH --partition=dcgp_usr_prod

#!/bin/bash
source "${CONDA_PREFIX}/etc/profile.d/conda.sh"
conda activate bioinfo

REFERENCE=/leonardo_work/IscrC_KREIOS/ENEO_setup/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
BAMS_PATH=/leonardo_work/IscrC_KREIOS/nanopore_cell_lines/TS576_11122025/
OUT_MAIN=/leonardo_scratch/large/userexternal/tducci00/nanopore_cell_lines/TS576_11122025/

BAMS=(ls ${BAMS_PATH}bam_pass/*.bam)
f=${BAMS[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename "$f" .bam)

samtools fastq "$f" \
| minimap2 -ax map-ont "$REFERENCE" - \
| samtools view -b \
| samtools sort -o "${OUT_MAIN}${SAMPLE_NAME}_aligned_sorted.bam"

conda deactivate

