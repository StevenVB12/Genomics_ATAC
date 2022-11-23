#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=144:00:00
#SBATCH --job-name=sra
#SBATCH --error=sra
#SBATCH --output=sra
#SBATCH --partition=rpapa
#SBATCH --ntasks=1

module load sratoolkit/2.10.3

prefetch --type fastq $(<SRA_accesions.txt)
