
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

for r1 in "$FQ"/*_1.fastq.gz; do
  s=${r1##*/}; s=${s%_1.fastq.gz}
  r2=$FQ/${s}_2.fastq.gz
  [[ -f "$r2" ]] || { echo "missing $r2" >&2; continue; }

  hisat2 -p "$T" --dta --known-splicesite-infile "$SPL" -x "$IDX" -1 "$r1" -2 "$r2" 2> "$OUT/$s.log" \
  | samtools view -@ "$T" -b \
  | samtools sort -@ "$T" -o "$OUT/$s.bam"

  samtools index "$OUT/$s.bam"
done