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

#name=(sections_D1_FPup_rest_erato sections_D1_FMup_rest_erato sections_D1_FDup_rest_erato sections_D1_FPdown_rest_erato sections_D1_FMdown_rest_erato sections_D1_FDdown_rest_erato sections_D2_FPup_rest_erato sections_D2_FMup_rest_erato sections_D2_FDup_rest_erato sections_D2_FPdown_rest_erato sections_D2_FMdown_rest_erato sections_D2_FDdown_rest_erato)

#~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").bed ~/share/REF/H_erato_dem/Herato_final.fasta -gff ~/share/REF/GFF/Heliconius_erato_demophoon_v1.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > F$(echo "${name[indID]}")_HOMER.txt


name=(sections_D1_FPup_rest_melp sections_D1_FMup_rest_melp sections_D1_FDup_rest_melp sections_D1_FPdown_rest_melp sections_D1_FMdown_rest_melp sections_D1_FDdown_rest_melp sections_D2_FPup_rest_melp sections_D2_FMup_rest_melp sections_D2_FDup_rest_melp sections_D2_FPdown_rest_melp sections_D2_FMdown_rest_melp sections_D2_FDdown_rest_melp)

~/share/programs/homer/bin/annotatePeaks.pl $(echo "${name[indID]}").bed /work/rpapa/share/REF/Hmel2/Hmel2.fa -gff ~/share/REF/GFF/Heliconius_melpomene_melpomene_Hmel2.gff3 -annStats $(echo "${name[indID]}")_blabla.txt > $(echo "${name[indID]}")_HOMER.txt