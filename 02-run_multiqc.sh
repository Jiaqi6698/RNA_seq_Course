#!/bin/bash

#SBATCH --partition=pibu_el8
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jiaqi.li1@students.unibe.ch
#SBATCH --output=/data/users/jli2/rnaseq/log/out/%j.out
#SBATCH --error=/data/users/jli2/rnaseq/log/err/%j.err

mkdir -p /data/users/jli2/rnaseq/results/multiqc

module load MultiQC

multiqc /data/users/jli2/rnaseq/results/fastqc \
    -o /data/users/jli2/rnaseq/results/multiqc \
    -n fastqc_report \
    --title "RNA-seq FastQC Quality Control Repor"
