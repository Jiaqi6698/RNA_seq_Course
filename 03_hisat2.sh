
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jiaqi.li1@students.unibe.ch
#SBATCH --output=/data/users/jli2/rnaseq/log/out/%j.out
#SBATCH --error=/data/users/jli2/rnaseq/log/err/%j.err

set -euo pipefail

BASE=/data/users/jli2/rnaseq
FQ=$BASE/dataset
OUT=$BASE/results/03_hisat2
IDX=$BASE/references/hisat2_index/mm39/mm39
SPL=$BASE/references/splicesites.txt
T=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$OUT"
module load HISAT2 SAMtools
# Loop over all paired-end FASTQ files (read 1)
for r1 in "$FQ"/*_1.fastq.gz; do
# Extract sample name from read 1 filename
  s=${r1##*/}; s=${s%_1.fastq.gz}
  # Define corresponding read 2 file
  r2=$FQ/${s}_2.fastq.gz
  # Check whether paired read 2 exists; skip sample if missing
  [[ -f "$r2" ]] || { echo "missing $r2" >&2; continue; }
# Align paired-end reads to the reference genome using HISAT2
  hisat2 -p "$T" --dta --known-splicesite-infile "$SPL" -x "$IDX" -1 "$r1" -2 "$r2" 2> "$OUT/$s.log" \
  # Convert SAM output to BAM format
  | samtools view -@ "$T" -b \
    # Sort BAM file by genomic coordinates
  | samtools sort -@ "$T" -o "$OUT/$s.bam"
# Index the sorted BAM file for downstream analysis
  samtools index "$OUT/$s.bam"
done