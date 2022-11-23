#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=48:00:00
#SBATCH --job-name=htsq
#SBATCH --error=htsq
#SBATCH --output=htsq
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=8

echo "test"

module load python3/3.6.10 

indID=$((SLURM_ARRAY_TASK_ID -1))

# erato 39
#name=(A4_HPo2f B2_FD2f B4_HA2f B8_FP1f C3_FD1f C4_FM2f C6_HA1f C9_FP2f D10_FP1f D6_HA2f D9_HPo2f E3_FM2f E6_FP2f F21_HPo1f F8_FM1f G21_FM1f G5_FM1f G6_FD2f H4_HA1f Heh017PAf Heh033DPf Heh033EDf Heh034IPf Heh034JDf Heh17DBf Heh25CPf Heh25DDf Heh36PCf Heh46HWHf Hep40FWCf Hep40HWDf Hep41FWC_Crawlerf Hep41HWD_Crawlerf Hep47HWHf I8_FD1f J4_FD1f J8_FP2f N3_FM2f S3_HPo1f)  

#htseq-count -c htseq/$(echo "${name[indID]}")_erato_nonStranded.counts -n 8 -s no --idattr Parent BAM/$(echo "${name[indID]}")_erato.filtered.sorted.nd.bam ~/share/REF/GFF/Heliconius_erato_demophoon_v1.gff3




# melp 50
name=(13Ff_ 13Gf_ 13Hf_ 13If_ 13Jf_ 14Af_ 14Bf_ 14Cf_ 14Df_ 14Ef_ 14Ff_ 14Gf_ 14Hf_ 14Jf_ 15Af_ 15Bf_ 15Cf_ 15Df_ 15Ef_ 16Af_ 16Bf_ 16Df_ 16Ef_ 16Hf_ 17Ff_ 17Gf_ 17Hf_ 17If_ 17Jf_ 2213Ff_ 2213Gf_ 2213Jf_ 24Df_ Hmm49DAf_ Hmm49PBf_ Hmm52EPf_ Hmm54GPf_ Hmm54HDf_ Hmm55HWJf_ Hmm56HWBf_ Hmm59DDf_ Hmm59PCf_ Hmm61DFf_ Hmm61PEf_ Hmm62DHf_ Hmm62PGf_ Hmr37HWBf_ Hmr44HWHf_ Hmr45FWEf_ Hmr45HWFf_)

htseq-count -c htseq/$(echo "${name[indID]}")melp_nonStranded.counts -n 8 -s no --idattr Parent BAM/$(echo "${name[indID]}")melp.filtered.sorted.nd.bam ~/share/REF/GFF/Heliconius_melpomene_melpomene_Hmel2.gff3

