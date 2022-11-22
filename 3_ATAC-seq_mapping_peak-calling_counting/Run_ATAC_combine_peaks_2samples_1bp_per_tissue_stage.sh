######
# Generate peak sets in H. erato for each tissue/time 1bp overlap in at least 2 samples
######

FILEL=(E4-FW FW-pboy LI7_demophoon_FW LB_41 LB_20 LB_27 LB_29 E3_FW)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > 5thFW_ERA_1bp.bed

sort -k1,1 -k2,2n 5thFW_ERA_1bp.bed > 5thFW_ERA_1bp.sort.bed
bedtools merge -i 5thFW_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > 5thFW_ERA_1bp.m1.sort.bed

FILEL=(E3_HW LB_42 LI7_demophoon_HW LB_21 LB_28 LB_30)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > 5thHW_ERA_1bp.bed

sort -k1,1 -k2,2n 5thHW_ERA_1bp.bed > 5thHW_ERA_1bp.sort.bed
bedtools merge -i 5thHW_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > 5thHW_ERA_1bp.m1.sort.bed

FILEL=(LB_45 LI21_demophoon_FD LI27_Demophoon_FD LB_24 LB_65 LI4_hydara_FD)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1FD_ERA_1bp.bed

sort -k1,1 -k2,2n D1FD_ERA_1bp.bed > D1FD_ERA_1bp.sort.bed
bedtools merge -i D1FD_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D1FD_ERA_1bp.m1.sort.bed

FILEL=(LB_44 LI21_demophoon_FM LI27_Demophoon_FM LB_23 LB_64 LI4_hydara_FM)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1FM_ERA_1bp.bed

sort -k1,1 -k2,2n D1FM_ERA_1bp.bed > D1FM_ERA_1bp.sort.bed
bedtools merge -i D1FM_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D1FM_ERA_1bp.m1.sort.bed

FILEL=(LB_43 LI21_demophoon_FP LI27_Demophoon_FP LB_26 LB_22 LB_63 LI4_hydara_FP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1FP_ERA_1bp.bed

sort -k1,1 -k2,2n D1FP_ERA_1bp.bed > D1FP_ERA_1bp.sort.bed
bedtools merge -i D1FP_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D1FP_ERA_1bp.m1.sort.bed

FILEL=(LB_46 LI21_demophoon_HA LI27_Demophoon_HA LB_47 LI20_demophoon_HP LI21_demophoon_HP LI27_Demophoon_HP LB_25 LB_66 LI4_hydara_HA LB_26 LB_67 LI4_hydara_HP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1HW_ERA_1bp.bed

sort -k1,1 -k2,2n D1HW_ERA_1bp.bed > D1HW_ERA_1bp.sort.bed
bedtools merge -i D1HW_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D1HW_ERA_1bp.m1.sort.bed

FILEL=(LB_3 LB_50 LI2_demophoon_FD LB_33 LB_38 LI1_hydara_FD)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2FD_ERA_1bp.bed

sort -k1,1 -k2,2n D2FD_ERA_1bp.bed > D2FD_ERA_1bp.sort.bed
bedtools merge -i D2FD_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D2FD_ERA_1bp.m1.sort.bed

FILEL=(LB_2 LB_49 LI2_demophoon_FM LB_32 LB_37 LI1_hydara_FM)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2FM_ERA_1bp.bed

sort -k1,1 -k2,2n D2FM_ERA_1bp.bed > D2FM_ERA_1bp.sort.bed
bedtools merge -i D2FM_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D2FM_ERA_1bp.m1.sort.bed

FILEL=(LB_1 LB_48 LI2_demophoon_FP LB_31 LB_36 LI1_hydara_FP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2FP_ERA_1bp.bed

sort -k1,1 -k2,2n D2FP_ERA_1bp.bed > D2FP_ERA_1bp.sort.bed
bedtools merge -i D2FP_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D2FP_ERA_1bp.m1.sort.bed

FILEL=(LB_4 LB_51 LI2_demophoon_HA LB_5 LB_52 LI2_demophoon_HP LB_34 LB_39 LI1_hydara_HA LB_35 LB_40 LI1_hydara_HP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2HW_ERA_1bp.bed

sort -k1,1 -k2,2n D2HW_ERA_1bp.bed > D2HW_ERA_1bp.sort.bed
bedtools merge -i D2HW_ERA_1bp.sort.bed -c 1 -o count | awk '$4>1' | cut -d$'\t' -f 1-3 > D2HW_ERA_1bp.m1.sort.bed



######
# Generate peak sets in H. melpomene for each tissue/time 1bp overlap in 2 samples
######

FILEL=(LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI14_rosina_HW M4-HW LB_68)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > 5thFW_MELP_1bp.bed

sort -k1,1 -k2,2n 5thFW_MELP_1bp.bed > 5thFW_MELP_1bp.sort.bed
bedtools merge -i 5thFW_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > 5thFW_MELP_1bp.m1.sort.bed

FILEL=(LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI14_rosina_HW M4-HW LB_69)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > 5thHW_MELP_1bp.bed

sort -k1,1 -k2,2n 5thHW_MELP_1bp.bed > 5thHW_MELP_1bp.sort.bed
bedtools merge -i 5thHW_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > 5thHW_MELP_1bp.m1.sort.bed

FILEL=(LI22_melpomene_FD LI23_melpomene_FD LI6_melpomene_FD LB_72 LB_77 LI13_rosina_FD)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1FD_MELP_1bp.bed

sort -k1,1 -k2,2n D1FD_MELP_1bp.bed > D1FD_MELP_1bp.sort.bed
bedtools merge -i D1FD_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FD_MELP_1bp.m1.sort.bed

FILEL=(LI22_melpomene_FM LI23_melpomene_FM LI6_melpomene_FM LB_71 LB_76 LI13_rosina_FM)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1FM_MELP_1bp.bed

sort -k1,1 -k2,2n D1FM_MELP_1bp.bed > D1FM_MELP_1bp.sort.bed
bedtools merge -i D1FM_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FM_MELP_1bp.m1.sort.bed

FILEL=(LI23_melpomene_FP LI6_melpomene_FP LB_70 LI13_rosina_FP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1FP_MELP_1bp.bed

sort -k1,1 -k2,2n D1FP_MELP_1bp.bed > D1FP_MELP_1bp.sort.bed
bedtools merge -i D1FP_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FP_MELP_1bp.m1.sort.bed

FILEL=(LI22_melpomene_HA LI23_melpomene_HA LI6_melpomene_HA LI22_melpomene_HP LI23_melpomene_HP LI6_melpomene_HP LB_78 LI13_rosina_HA LB_79 LI13_rosina_HP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D1HW_MELP_1bp.bed

sort -k1,1 -k2,2n D1HW_MELP_1bp.bed > D1HW_MELP_1bp.sort.bed
bedtools merge -i D1HW_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1HW_MELP_1bp.m1.sort.bed

FILEL=(LI25_melpomene_FD LI28_melpomene_FD LB_82 LI19_rosina_FD)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2FD_MELP_1bp.bed

sort -k1,1 -k2,2n D2FD_MELP_1bp.bed > D2FD_MELP_1bp.sort.bed
bedtools merge -i D2FD_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FD_MELP_1bp.m1.sort.bed

FILEL=(LI25_melpomene_FM LI28_melpomene_FM LB_81 LI19_rosina_FM2)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2FM_MELP_1bp.bed

sort -k1,1 -k2,2n D2FM_MELP_1bp.bed > D2FM_MELP_1bp.sort.bed
bedtools merge -i D2FM_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FM_MELP_1bp.m1.sort.bed

FILEL=(LI25_melpomene_FP LI28_melpomene_FP LB_80 LI19_rosina_FP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2FP_MELP_1bp.bed

sort -k1,1 -k2,2n D2FP_MELP_1bp.bed > D2FP_MELP_1bp.sort.bed
bedtools merge -i D2FP_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FP_MELP_1bp.m1.sort.bed

FILEL=(LI25_melpomene_HA LI28_melpomene_HA LI25_melpomene_HP LI28_melpomene_HP LB_83 LI19_rosina_HA LB_84 LI19_rosina_HP)
ALL_LIST=""
for FILE in ${FILEL[*]}
do
ALL_LIST="$ALL_LIST $FILE"_peaks.narrowPeak""
done
eval command=\$$(echo ALL_LIST)
cat $(echo $command) > D2HW_MELP_1bp.bed

sort -k1,1 -k2,2n D2HW_MELP_1bp.bed > D2HW_MELP_1bp.sort.bed
bedtools merge -i D2HW_MELP_1bp.sort.bed -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2HW_MELP_1bp.m1.sort.bed

