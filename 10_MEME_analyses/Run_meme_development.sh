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




#12
#name=(development_D1_p05FC1_erato_D1down development_D2_p05FC1_erato_D2up development_5th_p05FC1_erato_5thdown development_D1_p05FC1_erato_D1up development_5th_p05FC1_erato_5thup development_D2_p05FC1_shared_D2down development_D1_p05FC1_shared_D1down development_D2_p05FC1_shared_D2up development_5th_p05FC1_shared_5thdown development_D1_p05FC1_shared_D1up development_5th_p05FC1_shared_5thup development_D2_p05FC1_erato_D2down)

#sed -i 's/ /	/g' /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").bed

#bedtools getfasta -fi /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta -bed /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").bed > /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").fasta

#meme-chip -o /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}") -time 300 -ccut 0 -order 1 -db /work/rpapa/share/SVB_Angelo/TFmotifs2.txt -bfile /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/background_meme_markov/dem_hydara_peaks_Genrich_all_merged.m1.mod -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").fasta




#6
#name=(development_D1_p05FC1_melp_D1down development_D2_p05FC1_melp_D2up development_5th_p05FC1_melp_5thdown development_D1_p05FC1_melp_D1up development_5th_p05FC1_melp_5thup development_D2_p05FC1_melp_D2down)

#sed -i 's/ /	/g' /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").bed

#bedtools getfasta -fi /work/rpapa/share/REF/Hmel2/Hmel2.fa -bed /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").bed > /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").fasta

#meme-chip -o /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}") -time 300 -ccut 0 -order 1 -db /work/rpapa/share/SVB_Angelo/TFmotifs2.txt -bfile /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/background_meme_markov/dem_hydara_peaks_Genrich_all_merged.m1.mod -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/development/$(echo "${name[indID]}").fasta






