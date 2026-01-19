#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100GB
#SBATCH -p boost_usr_prod
#SBATCH -q normal
#SBATCH -o complete_variants_pipeline.out
#SBATCH -e complete_variants_pipeline.err
#SBATCH -A IscrC_KREIOS
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=t.ducci@student.unisi.it
export OMP_PROC_BIND=true


BAM_FOLDER=/leonardo_scratch/large/userexternal/tducci00/nanopore_cell_lines/GSC11_11122025/
BAM=${BAM_FOLDER}GSC11_merged_sorted.bam
SAMPLE_NAME=$(basename "$BAM" _merged_sorted.bam)

CACHE=/leonardo_work/IscrC_KREIOS/ENEO_setup/homo_sapiens_vep_105_GRCh38
REFERENCE=/leonardo_work/IscrC_KREIOS/ENEO_setup/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
PON=/leonardo_work/IscrC_KREIOS/ENEO_setup/af-only-gnomad.hg38.vcf.gz
INTERESTING_GENES=/leonardo/home/userexternal/tducci00/genes_coordinates_final.tsv

OUTDIR_CLAIR=${BAM_FOLDER}clair_out/
OUTDIR_PHASING=${BAM_FOLDER}phasing/
OUTDIR_SAVANA=${BAM_FOLDER}savana/
OUTDIR_ANNOT=${BAM_FOLDER}annot/
SUMMARY_FOLDER=${BAM_FOLDER}summary/

PATH_TO_LONGPHASE_SIF=/leonardo/home/userexternal/tducci00/longphase-to/longphase-to_latest.sif
PATH_TO_PLOT_SCRIPT=/leonardo/home/userexternal/tducci00/Scripts/variants/build_reports.py

THREADS=24
BIND="/leonardo_scratch:/leonardo_scratch,/leonardo_work:/leonardo_work,/leonardo:/leonardo"
source "${CONDA_PREFIX}/etc/profile.d/conda.sh"

### clair ##
PLATFORM=ont_r10_dorado_sup_5khz_ssrs

mkdir -p "$OUTDIR_CLAIR"


conda activate clairs-to

"$HOME/ClairS-TO/run_clairs_to" \
  --tumor_bam_fn "${BAM}" \
  --ref_fn "${REFERENCE}" \
  --threads "${THREADS}" \
  --platform "${PLATFORM}" \
  --output_dir "${OUTDIR_CLAIR}"

conda deactivate
conda activate bioinfo

#phasing
bcftools concat -a ${OUTDIR_CLAIR}snv.vcf.gz ${OUTDIR_CLAIR}indel.vcf.gz -Oz -o ${OUTDIR_CLAIR}all_variants.vcf.gz
bcftools index ${OUTDIR_CLAIR}all_variants.vcf.gz

mkdir -p ${OUTDIR_PHASING}

conda deactivate
conda activate longphase_env

singularity exec -B ${BIND} ${PATH_TO_LONGPHASE_SIF} longphase-to phase \
-s ${OUTDIR_CLAIR}all_variants.vcf.gz \
-c clairs_to_ssrs \
-r ${REFERENCE} \
-t ${THREADS} \
--indels \
-b ${BAM} \
--pon-file=${PON} \
-o ${OUTDIR_PHASING}snv_indel_phased

rm ${OUTDIR_CLAIR}all_variants.vcf.gz ${OUTDIR_CLAIR}all_variants.vcf.gz.tbi

#haplotag#
singularity exec -B ${BIND} ${PATH_TO_LONGPHASE_SIF} longphase-to haplotag \
-r ${REFERENCE} \
-s ${OUTDIR_PHASING}snv_indel_phased.vcf \
-b ${BAM} \
-t ${THREADS} \
-o ${OUTDIR_PHASING}${SAMPLE_NAME}_phased

conda deactivate

conda activate bioinfo
samtools index ${OUTDIR_PHASING}${SAMPLE_NAME}_phased.bam
conda deactivate

###VEP####
mkdir -p ${OUTDIR_ANNOT}

conda activate vep

vep --input_file ${OUTDIR_PHASING}snv_indel_phased.vcf --output_file ${OUTDIR_ANNOT}snv_indels_phased_annotated.vcf \
            --fasta ${REFERENCE} --offline --cache --dir_cache ${CACHE} \
            --assembly GRCh38 \
            --pick --af --check_existing --coding_only --format vcf --vcf --symbol --terms SO --no_intergenic --tsl --fork ${THREADS}

conda deactivate

###SAVANA########
conda activate savana

mkdir -p ${OUTDIR_SAVANA}

#CNA#
savana to --tumour ${OUTDIR_PHASING}${SAMPLE_NAME}_phased.bam \
 --outdir ${OUTDIR_SAVANA}CNA \
 --ref ${REFERENCE} \
 --g1000_vcf 1000g_hg38 \
 --threads ${THREADS}
#SV#
savana to --tumour ${BAM} \
 --outdir ${OUTDIR_SAVANA}to \
 --ref ${REFERENCE} \
 --g1000_vcf 1000g_hg38 \
 --threads ${THREADS}

conda deactivate

#intersecting with genes, initially removing header and weird stuff (POS=0 in vcf 4.2 means mapping to strange contigs, telomeres for example)
conda activate bioinfo
mkdir -p ${SUMMARY_FOLDER}

awk '$2 != 0' ${OUTDIR_SAVANA}to/${SAMPLE_NAME}_merged_sorted.classified.somatic.vcf > ${SUMMARY_FOLDER}/${SAMPLE_NAME}_fixed.vcf
bedtools intersect -wo -a ${SUMMARY_FOLDER}/${SAMPLE_NAME}_fixed.vcf -b ${INTERESTING_GENES} > ${SUMMARY_FOLDER}/SV_overlapping_genes.vcf
rm ${SUMMARY_FOLDER}/${SAMPLE_NAME}_fixed.vcf

tail -n +2 ${OUTDIR_SAVANA}CNA/${SAMPLE_NAME}_phased_segmented_absolute_copy_number.tsv > ${SUMMARY_FOLDER}/${SAMPLE_NAME}_copy_number_no_header.tsv
bedtools intersect -wo -a ${SUMMARY_FOLDER}/${SAMPLE_NAME}_copy_number_no_header.tsv -b ${INTERESTING_GENES} > ${SUMMARY_FOLDER}/CNV_overlapping_genes.tsv
rm ${SUMMARY_FOLDER}/${SAMPLE_NAME}_copy_number_no_header.tsv

conda deactivate
###plotting results####
conda activate data_analysis

python3 ${PATH_TO_PLOT_SCRIPT} \
 --summary_folder ${SUMMARY_FOLDER} \
 --snv_file ${OUTDIR_ANNOT}snv_indels_phased_annotated.vcf \
 --genes_coordinates_file ${INTERESTING_GENES} \
 --sample_name ${SAMPLE_NAME}

conda deactivate

