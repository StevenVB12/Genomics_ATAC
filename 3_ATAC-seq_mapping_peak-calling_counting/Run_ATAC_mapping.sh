#!/bin/bash
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=8:00:00
#SBATCH --job-name=trimm
#SBATCH --error=trimm
#SBATCH --output=trimm
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=8

module load bowtie2
module load samtools
module load jdk

ID=$((SLURM_ARRAY_TASK_ID -1))

#83
samples=(BR10_Demophoon_Brain BR14_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain E4-Head LB_3 LB_45 LB_50 LI2_demophoon_FD LI21_demophoon_FD LI27_Demophoon_FD LB_2 LB_44 LB_49 LI2_demophoon_FM LI21_demophoon_FM LI27_Demophoon_FP E3_FW E4-FW FW-pboy LB_41 LI7_demophoon_FW LB_1 LB_43 LB_48 LI2_demophoon_FP LI21_demophoon_FP LI27_Demophoon_FM LB_4 LB_46 LB_51 LI2_demophoon_HA LI21_demophoon_HA LI27_Demophoon_HA E3_HW LB_42 LI7_demophoon_HW LB_47 LB_5 LB_52 LI2_demophoon_HP LI20_demophoon_HP LI21_demophoon_HP LI27_Demophoon_HP BR2_hydara_brain BR3_hydara_brain BR6_Hydara_Brain LB_24 LB_33 LB_38 LB_65 LI1_hydara_FD LI4_hydara_FD LB_23 LB_32 LB_37 LB_64 LI1_hydara_FM LI4_hydara_FM LB_20 LB_27 LB_29 LB_22 LB_31 LB_36 LB_63 LI1_hydara_FP LI4_hydara_FP LB_25 LB_34 LB_39 LB_66 LI1_hydara_HA LI4_hydara_HA LB_21 LB_28 LB_30 LB_26 LB_35 LB_40 LB_67 LI1_hydara_HP LI4_hydara_HP)

#62 samples H. melpomene
#samples=(BR11_Rosina_Brain BR5_Rosina_Brain M4-Head LB_72 LB_77 LB_82 LI13_rosina_FD LI19_rosina_FD LB_71 LB_76 LB_81 LI13_rosina_FM LI19_rosina_FM2 LB_68 LI14_rosina_FW M4-FW LB_70 LB_80 LI13_rosina_FP LI19_rosina_FP LB_78 LB_83 LI13_rosina_HA LI19_rosina_HA LB_69 LI14_rosina_HW M4-HW LB_79 LB_84 LI13_rosina_HP LI19_rosina_HP LI22_melpomene_FD LI23_melpomene_FD LI25_melpomene_FD LI28_melpomene_FD LI6_melpomene_FD LI22_melpomene_FM LI23_melpomene_FM LI25_melpomene_FM LI28_melpomene_FM LI6_melpomene_FM LI12_melpomene_FW LI15_melpomene_FW LI8_melpomene_FW LI22_melpomene_FP LI23_melpomene_FP LI25_melpomene_FP LI28_melpomene_FP LI6_melpomene_FP LI22_melpomene_HA LI23_melpomene_HA LI25_melpomene_HA LI28_melpomene_HA LI6_melpomene_HA LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI22_melpomene_HP LI23_melpomene_HP LI25_melpomene_HP LI28_melpomene_HP LI6_melpomene_HP)


REF=REF/H_erato_dem/Herato_final_bowtie
#REF=REF/Hmel2/Hmel2_bowtie

REFNAME=Herato
#REFNAME=Hmel2

SIZES=/work/rpapa/share/REF/H_erato_dem/Herato_final.fasta.sizes
#SIZES=/work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes

#### Run bowtie

gunzip reads_ATACseq_trimmomatic/$(echo "${samples[ID]}")_trim_R1.fastq.gz
gunzip reads_ATACseq_trimmomatic/$(echo "${samples[ID]}")_trim_R2.fastq.gz

bowtie2 -t -k 2 -p 8 --local -x $REF \
-1 reads_ATACseq_trimmomatic/$(echo "${samples[ID]}")_trim_R1.fastq \
-2 reads_ATACseq_trimmomatic/$(echo "${samples[ID]}")_trim_R2.fastq |\
samtools view -bS - > BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.bam

gzip reads_ATACseq_trimmomatic/$(echo "${samples[ID]}")_trim_R1.fastq
gzip reads_ATACseq_trimmomatic/$(echo "${samples[ID]}")_trim_R2.fastq


#### filter bam files

samtools view -f 0x02 -q 20 -b BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.bam > BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.bam

samtools sort BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.bam BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted

java -jar /work/rpapa/sbelleghem/Programs/picard-tools-2.5.0/picard.jar MarkDuplicates \
I=BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.bam \
O=BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.nd.bam \
Remove_Duplicates=true  M=BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_dup_metrics.txt ASSUME_SORTED=true

rm BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.bam
rm BAM_ATACseq_trimmomatic/$(echo "${samples[ID]}")_$REFNAME.trim.filtered.sorted.bam

