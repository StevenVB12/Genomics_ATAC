#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=12:00:00
#SBATCH --job-name=homer
#SBATCH --error=homer
#SBATCH --output=homer
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

indID=$((SLURM_ARRAY_TASK_ID -1))

export PATH="$PATH:~/share/programs/homer/bin"


#name=(FWHW_5th_p05FC1_erato_FWup FWHW_5th_p05FC1_erato_HWup FWHW_D2_p05FC1_erato_FWup FWHW_D2_p05FC1_erato_HWup FWHW_D1_p05FC1_erato_FWup FWHW_D1_p05FC1_erato_HWup)

#~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").txt ~/share/REF/H_erato_dem/Herato_final.fasta -gff ~/share/REF/GFF/Heliconius_erato_demophoon_v1.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > $(echo "${name[indID]}")_HOMER.txt


#name=(FWHW_D1_p05FC1_melp_FWup FWHW_D1_p05FC1_melp_HWup FWHW_5th_p05FC1_melp_FWup FWHW_5th_p05FC1_melp_HWup FWHW_D2_p05FC1_melp_FWup FWHW_D2_p05FC1_melp_HWup)

#~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").txt /work/rpapa/share/REF/Hmel2/Hmel2.fa -gff ~/share/REF/GFF/Heliconius_melpomene_melpomene_Hmel2.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > $(echo "${name[indID]}")_HOMER.txt

name=(FWHW_5th_p05FC1_shared_FWup_eratoCoords FWHW_5th_p05FC1_shared_HWup_eratoCoords FWHW_D1_p05FC1_shared_FWup_eratoCoords FWHW_D1_p05FC1_shared_HWup_eratoCoords FWHW_D2_p05FC1_shared_FWup_eratoCoords FWHW_D2_p05FC1_shared_HWup_eratoCoords)

~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").bed ~/share/REF/H_erato_dem/Herato_final.fasta -gff ~/share/REF/GFF/Heliconius_erato_demophoon_v1.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > $(echo "${name[indID]}")_HOMER.txt
