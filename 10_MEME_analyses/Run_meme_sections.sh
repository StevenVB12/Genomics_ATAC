#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=12:00:00
#SBATCH --job-name=meme
#SBATCH --error=meme
#SBATCH --output=meme
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

indID=$((SLURM_ARRAY_TASK_ID -1))



module load gcc/6.3.0
module load memesuite/5.1.1
module load bedtools

#8
name=(D1_grad_neg_era D1_grad_neg_era_shared FM_pos_era D1_grad_pos_era FP_pos_era D1_grad_pos_era_shared FP_pos_era_shared FD_pos_era)

sed -i 's/ /	/g' /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.bed

bedtools getfasta -fi /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta -bed /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.bed > /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.fasta

meme-chip -o /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}") -time 300 -ccut 0 -order 1 -db /work/rpapa/share/SVB_Angelo/TFmotifs2.txt -bfile /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/background_meme_markov/dem_hydara_peaks_Genrich_all_merged.m1.mod -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.fasta

#5
#name=(FD_pos_melp D1_grad_neg_melp FM_pos_melp D1_grad_pos_melp FP_pos_melp)

#sed -i 's/ /	/g' /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.bed

#bedtools getfasta -fi /work/rpapa/share/REF/Hmel2/Hmel2.fa -bed /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.bed > /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.fasta

#meme-chip -o /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}") -time 300 -ccut 0 -order 1 -db /work/rpapa/share/SVB_Angelo/TFmotifs2.txt -bfile /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/background_meme_markov/dem_hydara_peaks_Genrich_all_merged.m1.mod -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/sections2/$(echo "${name[indID]}")_06122021.fasta
