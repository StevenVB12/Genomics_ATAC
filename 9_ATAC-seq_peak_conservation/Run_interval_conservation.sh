#!/bin/bash
#SBATCH --mem-per-cpu=14gb
#SBATCH --time=48:00:00
#SBATCH --job-name=IDY
#SBATCH --error=IDY
#SBATCH --output=IDY
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load python2

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_1bpOverlap.bed -o 5thup_shared_PANpos_1bpOverlap_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_1bpOverlap.bed -o 5thup_eratoUnique_PANpos_1bpOverlap_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_1bpOverlap.bed -o 5thup_melpUnique_PANpos_1bpOverlap_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_1bpOverlap.bed -o D1up_shared_PANpos_1bpOverlap_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_1bpOverlap.bed -o D1up_eratoUnique_PANpos_1bpOverlap_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_1bpOverlap.bed -o D1up_melpUnique_PANpos_1bpOverlap_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_1bpOverlap.bed -o D2up_shared_PANpos_1bpOverlap_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_1bpOverlap.bed -o D2up_eratoUnique_PANpos_1bpOverlap_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_1bpOverlap.bed -o D2up_melpUnique_PANpos_1bpOverlap_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_1bpOverlap_strict.bed -o 5thup_shared_PANpos_1bpOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_1bpOverlap_strict.bed -o 5thup_eratoUnique_PANpos_1bpOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_1bpOverlap_strict.bed -o 5thup_melpUnique_PANpos_1bpOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_1bpOverlap_strict.bed -o D1up_shared_PANpos_1bpOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_1bpOverlap_strict.bed -o D1up_eratoUnique_PANpos_1bpOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_1bpOverlap_strict.bed -o D1up_melpUnique_PANpos_1bpOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_1bpOverlap_strict.bed -o D2up_shared_PANpos_1bpOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_1bpOverlap_strict.bed -o D2up_eratoUnique_PANpos_1bpOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_1bpOverlap_strict.bed -o D2up_melpUnique_PANpos_1bpOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_10percOverlap_strict.bed -o 5thup_shared_PANpos_10percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_10percOverlap_strict.bed -o D1up_shared_PANpos_10percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_10percOverlap_strict.bed -o D2up_shared_PANpos_10percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_30percOverlap_strict.bed -o 5thup_shared_PANpos_30percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_30percOverlap_strict.bed -o D1up_shared_PANpos_30percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_30percOverlap_strict.bed -o D2up_shared_PANpos_30percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_50percccOverlap_strict.bed -o 5thup_shared_PANpos_50percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_50percccOverlap_strict.bed -o D1up_shared_PANpos_50percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_50percccOverlap_strict.bed -o D2up_shared_PANpos_50percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_70perccccOverlap_strict.bed -o 5thup_shared_PANpos_70percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_70perccccOverlap_strict.bed -o D1up_shared_PANpos_70percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_70perccccOverlap_strict.bed -o D2up_shared_PANpos_70percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_10percOverlap_strict.bed -o 5thup_eratoUnique_PANpos_10percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_10percOverlap_strict.bed -o D1up_eratoUnique_PANpos_10percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_10percOverlap_strict.bed -o D2up_eratoUnique_PANpos_10percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_10percOverlap_strict.bed -o 5thup_melpUnique_PANpos_10percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_10percOverlap_strict.bed -o D1up_melpUnique_PANpos_10percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_10percOverlap_strict.bed -o D2up_melpUnique_PANpos_10percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_30percOverlap_strict.bed -o 5thup_eratoUnique_PANpos_30percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_30percOverlap_strict.bed -o D1up_eratoUnique_PANpos_30percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_30percOverlap_strict.bed -o D2up_eratoUnique_PANpos_30percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_30percOverlap_strict.bed -o 5thup_melpUnique_PANpos_30percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_30percOverlap_strict.bed -o D1up_melpUnique_PANpos_30percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_30percOverlap_strict.bed -o D2up_melpUnique_PANpos_30percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_50percccOverlap_strict.bed -o 5thup_eratoUnique_PANpos_50percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_50percccOverlap_strict.bed -o D1up_eratoUnique_PANpos_50percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_50percccOverlap_strict.bed -o D2up_eratoUnique_PANpos_50percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_50percccOverlap_strict.bed -o 5thup_melpUnique_PANpos_50percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_50percccOverlap_strict.bed -o D1up_melpUnique_PANpos_50percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_50percccOverlap_strict.bed -o D2up_melpUnique_PANpos_50percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_70perccccOverlap_strict.bed -o 5thup_eratoUnique_PANpos_70percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_70perccccOverlap_strict.bed -o D1up_eratoUnique_PANpos_70percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_70perccccOverlap_strict.bed -o D2up_eratoUnique_PANpos_70percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_70perccccOverlap_strict.bed -o 5thup_melpUnique_PANpos_70percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_70perccccOverlap_strict.bed -o D1up_melpUnique_PANpos_70percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_70perccccOverlap_strict.bed -o D2up_melpUnique_PANpos_70percOverlap_IDY_strict


#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_25percOverlap_strict.bed -o 5thup_shared_PANpos_25percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_25percOverlap_strict.bed -o D1up_shared_PANpos_25percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_25percOverlap_strict.bed -o D2up_shared_PANpos_25percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_75percOverlap_strict.bed -o 5thup_shared_PANpos_75percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_75percOverlap_strict.bed -o D1up_shared_PANpos_75percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_75percOverlap_strict.bed -o D2up_shared_PANpos_75percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_shared_PANpos_100percOverlap_strict.bed -o 5thup_shared_PANpos_100percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_shared_PANpos_100percOverlap_strict.bed -o D1up_shared_PANpos_100percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_shared_PANpos_100percOverlap_strict.bed -o D2up_shared_PANpos_100percOverlap_IDY_strict


#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_25percOverlap_strict.bed -o 5thup_eratoUnique_PANpos_25percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_25percOverlap_strict.bed -o D1up_eratoUnique_PANpos_25percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_25percOverlap_strict.bed -o D2up_eratoUnique_PANpos_25percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_25percOverlap_strict.bed -o 5thup_melpUnique_PANpos_25percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_25percOverlap_strict.bed -o D1up_melpUnique_PANpos_25percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_25percOverlap_strict.bed -o D2up_melpUnique_PANpos_25percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_75percOverlap_strict.bed -o 5thup_eratoUnique_PANpos_75percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_75percOverlap_strict.bed -o D1up_eratoUnique_PANpos_75percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_75percOverlap_strict.bed -o D2up_eratoUnique_PANpos_75percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_75percOverlap_strict.bed -o 5thup_melpUnique_PANpos_75percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_75percOverlap_strict.bed -o D1up_melpUnique_PANpos_75percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_75percOverlap_strict.bed -o D2up_melpUnique_PANpos_75percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_eratoUnique_PANpos_100percOverlap_strict.bed -o 5thup_eratoUnique_PANpos_100percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_eratoUnique_PANpos_100percOverlap_strict.bed -o D1up_eratoUnique_PANpos_100percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_eratoUnique_PANpos_100percOverlap_strict.bed -o D2up_eratoUnique_PANpos_100percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5thup_melpUnique_PANpos_100percOverlap_strict.bed -o 5thup_melpUnique_PANpos_100percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1up_melpUnique_PANpos_100percOverlap_strict.bed -o D1up_melpUnique_PANpos_100percOverlap_IDY_strict
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2up_melpUnique_PANpos_100percOverlap_strict.bed -o D2up_melpUnique_PANpos_100percOverlap_IDY_strict

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5th_p05FC1_erato_FWup_pan.bed -o 5th_p05FC1_erato_FWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5th_p05FC1_erato_HWup_pan.bed -o 5th_p05FC1_erato_HWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5th_p05FC1_melp_FWup_pan.bed -o 5th_p05FC1_melp_FWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b 5th_p05FC1_melp_HWup_pan.bed -o 5th_p05FC1_melp_HWup_pan_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1_p05FC1_erato_FWup_pan.bed -o D1_p05FC1_erato_FWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1_p05FC1_erato_HWup_pan.bed -o D1_p05FC1_erato_HWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1_p05FC1_melp_FWup_pan.bed -o D1_p05FC1_melp_FWup_pan_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D1_p05FC1_melp_HWup_pan.bed -o D1_p05FC1_melp_HWup_pan_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2_p05FC1_erato_FWup_pan.bed -o D2_p05FC1_erato_FWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2_p05FC1_erato_HWup_pan.bed -o D2_p05FC1_erato_HWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2_p05FC1_melp_FWup_pan.bed -o D2_p05FC1_melp_FWup_pan_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b D2_p05FC1_melp_HWup_pan.bed -o D2_p05FC1_melp_HWup_pan_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_5th_p05FC1_shared_FWup_panCoords.bed -o FWHW_5th_p05FC1_shared_FWup_panCoords_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_5th_p05FC1_shared_HWup_panCoords.bed -o FWHW_5th_p05FC1_shared_HWup_panCoords_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_D1_p05FC1_shared_FWup_panCoords.bed -o FWHW_D1_p05FC1_shared_FWup_panCoords_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_D1_p05FC1_shared_HWup_panCoords.bed -o FWHW_D1_p05FC1_shared_HWup_panCoords_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_D2_p05FC1_shared_FWup_panCoords.bed -o FWHW_D2_p05FC1_shared_FWup_panCoords_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_D2_p05FC1_shared_HWup_panCoords.bed -o FWHW_D2_p05FC1_shared_HWup_panCoords_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_shared_PAN.bed -o FWHW_shared_PAN_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b FWHW_unique_PAN.bed -o FWHW_unique_PAN_IDY

#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b Section_unique_PAN.bed -o Section_unique_PAN_IDY
#python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b Section_shared_PAN.bed -o Section_shared_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_5th_FWHW_FWup.bed -o sign_E_5th_FWHW_FWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_5th_FWHW_FWup.bed -o sign_M_5th_FWHW_FWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_5th_FWHW_FWup.bed -o sign_SHARED_5th_FWHW_FWup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_5th_FWHW_HWup.bed -o sign_E_5th_FWHW_HWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_5th_FWHW_HWup.bed -o sign_M_5th_FWHW_HWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_5th_FWHW_HWup.bed -o sign_SHARED_5th_FWHW_HWup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FWHW_FWup.bed -o sign_E_D1_FWHW_FWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_FWHW_FWup.bed -o sign_M_D1_FWHW_FWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D1_FWHW_FWup.bed -o sign_SHARED_D1_FWHW_FWup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FWHW_HWup.bed -o sign_E_D1_FWHW_HWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_FWHW_HWup.bed -o sign_M_D1_FWHW_HWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D1_FWHW_HWup.bed -o sign_SHARED_D1_FWHW_HWup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D2_FWHW_FWup.bed -o sign_E_D2_FWHW_FWup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D2_FWHW_HWup.bed -o sign_E_D2_FWHW_HWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D2_FWHW_HWup.bed -o sign_M_D2_FWHW_HWup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D2_FWHW_HWup.bed -o sign_SHARED_D2_FWHW_HWup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_5th_D1D2_5thup.bed -o sign_E_5th_D1D2_5thup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_5th_D1D2_5thup.bed -o sign_M_5th_D1D2_5thup_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_5th_D1D2_5thup.bed -o sign_SHARED_5th_D1D2_5thup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_5th_D1D2_5thdown.bed -o sign_E_5th_D1D2_5thdown_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_5th_D1D2_5thdown.bed -o sign_M_5th_D1D2_5thdown_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_5th_D1D2_5thdown.bed -o sign_SHARED_5th_D1D2_5thdown_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_5thD2_D1up.bed -o sign_E_D1_5thD2_D1up_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_5thD2_D1up.bed -o sign_M_D1_5thD2_D1up_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D1_5thD2_D1up.bed -o sign_SHARED_D1_5thD2_D1up_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_5thD2_D1down.bed -o sign_E_D1_5thD2_D1down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_5thD2_D1down.bed -o sign_M_D1_5thD2_D1down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D1_5thD2_D1down.bed -o sign_SHARED_D1_5thD2_D1down_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D2_5thD1_D2up.bed -o sign_E_D2_5thD1_D2up_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D2_5thD1_D2up.bed -o sign_M_D2_5thD1_D2up_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D2_5thD1_D2up.bed -o sign_SHARED_D2_5thD1_D2up_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D2_5thD1_D2down.bed -o sign_E_D2_5thD1_D2down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D2_5thD1_D2down.bed -o sign_M_D2_5thD1_D2down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_D2_5thD1_D2down.bed -o sign_SHARED_D2_5thD1_D2down_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_grad_pos.bed -o sign_E_E_D1_grad_pos_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_grad_pos.bed -o sign_M_D1_grad_pos_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_gradup.bed -o sign_SHARED_gradup_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_grad_neg.bed -o sign_E_E_D1_grad_neg_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_grad_neg.bed -o sign_M_D1_grad_neg_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_graddown.bed -o sign_SHARED_graddown_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FD_rest_up.bed -o sign_E_D1_FD_rest_up_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FD_rest_down.bed -o sign_E_D1_FD_rest_down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_FD_rest_down.bed -o sign_M_D1_FD_rest_down_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FM_rest_up.bed -o sign_E_D1_FM_rest_up_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FP_rest_up.bed -o sign_E_D1_FP_rest_up_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_FP_rest_up.bed -o sign_M_D1_FP_rest_up_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_FPdown.bed -o sign_SHARED_FPdown_PAN_IDY

python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_E_D1_FP_rest_down.bed -o sign_E_D1_FP_rest_down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_M_D1_FP_rest_down.bed -o sign_M_D1_FP_rest_down_PAN_IDY
python seq-seq-pan_bedfile_conservation.py -I /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/erato_melp_pan_noNewline.fasta -g 1,2 -b sign_SHARED_FPup.bed -o sign_SHARED_FPup_PAN_IDY


