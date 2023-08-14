#!/bin/bash
#SBATCH --job-name=mbcL1
#SBATCH --output=mbcL1.o
#SBATCH --error=mbcL1.e
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
 --genomeDir /path/to/index/folder/mouse_index \
 --outSAMtype BAM SortedByCoordinate \
 --readFilesIn /path/to/fastq/file/R1/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R1_001.fastq /path/to/fastq/file/R2/SC3_v3_NextGem_DI_Neurons_5K_gex_S3_L001_R2_001.fastq \
 --outFileNamePrefix mbcL1_ \

 #mbcL1 = mouse brain cells L001 (Lane 1)
