#!/bin/bash
#SBATCH --job-name=mouse_index
#SBATCH --output=m_index.o
#SBATCH --error=m_index.e
#SBATCH --time=24:00:00
#SBATCH --qos=high
#SBATCH --ntasks=16
#SBATCH --mem=128gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb

# STAR version=2.7.10b
module load conda
source activate mapping
STAR --runThreadN 16 \
 --runMode genomeGenerate \
 --genomeDir /path/to/output_folder/mouse_index \
 --genomeFastaFiles /path/to/fasta/file/mouse_pri_assembly_genome.fa \
 --sjdbGTFfile /path/to/genome/annotation/file/mouse_pri_assembly_annotation.gtf \
 --sjdbOverhang 99 \
