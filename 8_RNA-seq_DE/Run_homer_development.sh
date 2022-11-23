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

#name=(development_5th_p05FC1_eratoUnique_5thdown development_5th_p05FC1_eratoUnique_5thup development_5th_p05FC1_shared_5thdown_eratoCoords development_5th_p05FC1_shared_5thup_eratoCoords development_D1_p05FC1_eratoUnique_D1down development_D1_p05FC1_eratoUnique_D1up development_D1_p05FC1_shared_D1down_eratoCoords development_D1_p05FC1_shared_D1up_eratoCoords development_D2_p05FC1_eratoUnique_D2down development_D2_p05FC1_eratoUnique_D2up development_D2_p05FC1_shared_D2down_eratoCoords development_D2_p05FC1_shared_D2up_eratoCoords)

#~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").bed ~/share/REF/H_erato_dem/Herato_final.fasta -gff ~/share/REF/GFF/Heliconius_erato_demophoon_v1.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > F$(echo "${name[indID]}")_HOMER.txt


name=(development_5th_p05FC1_melpUnique_5thdown development_5th_p05FC1_melpUnique_5thup development_5th_p05FC1_shared_5thdown_melpCoords development_5th_p05FC1_shared_5thup_melpCoords development_D1_p05FC1_melpUnique_D1down development_D1_p05FC1_melpUnique_D1up development_D1_p05FC1_shared_D1down_melpCoords development_D1_p05FC1_shared_D1up_melpCoords development_D2_p05FC1_melpUnique_D2down development_D2_p05FC1_melpUnique_D2up development_D2_p05FC1_shared_D2down_melpCoords development_D2_p05FC1_shared_D2up_melpCoords)

~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").bed /work/rpapa/share/REF/Hmel2/Hmel2.fa -gff ~/share/REF/GFF/Heliconius_melpomene_melpomene_Hmel2.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > $(echo "${name[indID]}")_HOMER.txt