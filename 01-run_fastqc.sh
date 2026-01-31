#!/bin/bash

#SBATCH --partition=pibu_el8
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=500MB
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jiaqi.li1@students.unibe.ch
#SBATCH --output=/data/users/jli2/rnaseq/log/out/%j.out
#SBATCH --output=/data/users/jli2/rnaseq/log/err/%j.err

mkdir -p /data/users/jli2/rnaseq/results/fastqc

module load FastQC/0.11.9-Java-11

fastqc /data/users/jli2/rnaseq/dataset/*.fastq.gz -o /data/users/jli2/rnaseq/results/fastqc