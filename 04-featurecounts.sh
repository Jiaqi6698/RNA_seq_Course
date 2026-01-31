#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/data/users/jli2/rnaseq/log/out/%j.out
#SBATCH --error=/data/users/jli2/rnaseq/log/err/%j.err

set -euo pipefail

BASE=/data/users/jli2/rnaseq
BAM=$BASE/results/03_hisat2
OUT=$BASE/results/03_featureCounts
GTF=$BASE/references/Mus_musculus.GRCm39.109.gtf
T=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$OUT"
module load Subread

# assign paired-end RNA-seq reads to genes by counting 
# exon-overlapping reads with sufficient alignment quality, 
# using a non-strand-specific library setting and multi-threaded processing.
featureCounts -T "$T" -p -s 0 -Q 10 -t exon -g gene_id \
  -a "$GTF" \
  -o "$OUT/counts_matrix.txt" \
  "$BAM"/*.bam