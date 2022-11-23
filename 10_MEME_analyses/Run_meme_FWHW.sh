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
#name=(FWHW_5th_p05FC1_shared_FWup FWHW_5th_p05FC1_erato_FWup FWHW_5th_p05FC1_shared_HWup FWHW_D1_p05FC1_shared_FWup FWHW_5th_p05FC1_erato_HWup FWHW_D1_p05FC1_erato_FWup FWHW_D1_p05FC1_shared_HWup FWHW_D2_p05FC1_shared_FWup FWHW_D1_p05FC1_erato_HWup FWHW_D2_p05FC1_erato_FWup FWHW_D2_p05FC1_shared_HWup FWHW_D2_p05FC1_erato_HWup) 

#sed -i 's/ /	/g' /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").bed

#bedtools getfasta -fi /work/rpapa/share/REF/H_erato_dem/Herato_final.fasta -bed /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").bed > /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").fasta

#meme-chip -o /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}") -time 300 -ccut 0 -order 1 -db /work/rpapa/share/SVB_Angelo/TFmotifs2.txt -bfile /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/background_meme_markov/dem_hydara_peaks_Genrich_all_merged.m1.mod -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").fasta




#6
#name=(FWHW_D1_p05FC1_melp_HWup FWHW_D2_p05FC1_melp_FWup FWHW_D2_p05FC1_melp_HWup FWHW_5th_p05FC1_melp_FWup FWHW_5th_p05FC1_melp_HWup FWHW_D1_p05FC1_melp_FWup)

#sed -i 's/ /	/g' /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").bed

#bedtools getfasta -fi /work/rpapa/share/REF/Hmel2/Hmel2.fa -bed /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").bed > /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").fasta

#meme-chip -o /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}") -time 300 -ccut 0 -order 1 -db /work/rpapa/share/SVB_Angelo/TFmotifs2.txt -bfile /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/background_meme_markov/dem_hydara_peaks_Genrich_all_merged.m1.mod -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 10 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/MEME_out/FW_HW/$(echo "${name[indID]}").fasta







