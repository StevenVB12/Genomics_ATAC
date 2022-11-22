######
# Generate peak sets in H. erato for each tissue/time 50% reciprocal overlap in at least 2 samples
######

FILEL="E4-FW FW-pboy LI7_demophoon_FW LB_41 LB_20 LB_27 LB_29 E3_FW"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.5thFW_50Overlap.common
fi
done
done

cat *5thFW_50Overlap.common > 5thFW_50Overlap.all.common
sort -k1,1 -k2,2n 5thFW_50Overlap.all.common > 5thFW_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' 5thFW_50Overlap.all.common.sort > 5thFW_50Overlap.all.common.sort.clean
bedtools merge -i 5thFW_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > 5thFW_50Overlap.all.common.sort.clean.m1.bed


FILEL="E3_HW LB_42 LI7_demophoon_HW LB_21 LB_28 LB_30"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.5thHW_50Overlap.common
fi
done
done

cat *5thHW_50Overlap.common > 5thHW_50Overlap.all.common
sort -k1,1 -k2,2n 5thHW_50Overlap.all.common > 5thHW_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' 5thHW_50Overlap.all.common.sort > 5thHW_50Overlap.all.common.sort.clean
bedtools merge -i 5thHW_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > 5thHW_50Overlap.all.common.sort.clean.m1.bed



FILEL="LB_45 LI21_demophoon_FD LI27_Demophoon_FD LB_24 LB_65 LI4_hydara_FD"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1FD_50Overlap.common
fi
done
done

cat *D1FD_50Overlap.common > D1FD_50Overlap.all.common
sort -k1,1 -k2,2n D1FD_50Overlap.all.common > D1FD_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1FD_50Overlap.all.common.sort > D1FD_50Overlap.all.common.sort.clean
bedtools merge -i D1FD_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FD_50Overlap.all.common.sort.clean.m1.bed


FILEL="LB_44 LI21_demophoon_FM LI27_Demophoon_FM LB_23 LB_64 LI4_hydara_FM"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1FM_50Overlap.common
fi
done
done

cat *D1FM_50Overlap.common > D1FM_50Overlap.all.common
sort -k1,1 -k2,2n D1FM_50Overlap.all.common > D1FM_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1FM_50Overlap.all.common.sort > D1FM_50Overlap.all.common.sort.clean
bedtools merge -i D1FM_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FM_50Overlap.all.common.sort.clean.m1.bed


FILEL="LB_43 LI21_demophoon_FP LI27_Demophoon_FP LB_26 LB_22 LB_63 LI4_hydara_FP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1FP_50Overlap.common
fi
done
done

cat *D1FP_50Overlap.common > D1FP_50Overlap.all.common
sort -k1,1 -k2,2n D1FP_50Overlap.all.common > D1FP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1FP_50Overlap.all.common.sort > D1FP_50Overlap.all.common.sort.clean
bedtools merge -i D1FP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FP_50Overlap.all.common.sort.clean.m1.bed



FILEL="LB_46 LI21_demophoon_HA LI27_Demophoon_HA LB_47 LI20_demophoon_HP LI21_demophoon_HP LI27_Demophoon_HP LB_25 LB_66 LI4_hydara_HA LB_26 LB_67 LI4_hydara_HP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1HW_50Overlap.common
fi
done
done

cat *D1HW_50Overlap.common > D1HW_50Overlap.all.common
sort -k1,1 -k2,2n D1HW_50Overlap.all.common > D1HW_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1HW_50Overlap.all.common.sort > D1HW_50Overlap.all.common.sort.clean
bedtools merge -i D1HW_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1HW_50Overlap.all.common.sort.clean.m1.bed


FILEL="LB_3 LB_50 LI2_demophoon_FD LB_33 LB_38 LI1_hydara_FD"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2FD_50Overlap.common
fi
done
done

cat *D2FD_50Overlap.common > D2FD_50Overlap.all.common
sort -k1,1 -k2,2n D2FD_50Overlap.all.common > D2FD_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2FD_50Overlap.all.common.sort > D2FD_50Overlap.all.common.sort.clean
bedtools merge -i D2FD_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FD_50Overlap.all.common.sort.clean.m1.bed


FILEL="LB_2 LB_49 LI2_demophoon_FM LB_32 LB_37 LI1_hydara_FM"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2FM_50Overlap.common
fi
done
done

cat *D2FM_50Overlap.common > D2FM_50Overlap.all.common
sort -k1,1 -k2,2n D2FM_50Overlap.all.common > D2FM_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2FM_50Overlap.all.common.sort > D2FM_50Overlap.all.common.sort.clean
bedtools merge -i D2FM_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FM_50Overlap.all.common.sort.clean.m1.bed



FILEL="LB_1 LB_48 LI2_demophoon_FP LB_31 LB_36 LI1_hydara_FP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2FP_50Overlap.common
fi
done
done

cat *D2FP_50Overlap.common > D2FP_50Overlap.all.common
sort -k1,1 -k2,2n D2FP_50Overlap.all.common > D2FP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2FP_50Overlap.all.common.sort > D2FP_50Overlap.all.common.sort.clean
bedtools merge -i D2FP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FP_50Overlap.all.common.sort.clean.m1.bed



FILEL="LB_4 LB_51 LI2_demophoon_HA LB_5 LB_52 LI2_demophoon_HP LB_34 LB_39 LI1_hydara_HA LB_35 LB_40 LI1_hydara_HP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2HW_50Overlap.common
fi
done
done

cat *D2HW_50Overlap.common > D2HW_50Overlap.all.common
sort -k1,1 -k2,2n D2HW_50Overlap.all.common > D2HW_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2HW_50Overlap.all.common.sort > D2HW_50Overlap.all.common.sort.clean
bedtools merge -i D2HW_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2HW_50Overlap.all.common.sort.clean.m1.bed



######
# Generate peak sets in H. erato for each tissue/time 50% reciprocal overlap in at least 2 samples
######

FILEL="LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI14_rosina_HW M4-HW LB_68"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.5thFW_MELP_50Overlap.common
fi
done
done

cat *5thFW_MELP_50Overlap.common > 5thFW_MELP_50Overlap.all.common
sort -k1,1 -k2,2n 5thFW_MELP_50Overlap.all.common > 5thFW_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' 5thFW_MELP_50Overlap.all.common.sort > 5thFW_MELP_50Overlap.all.common.sort.clean
bedtools merge -i 5thFW_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > 5thFW_MELP_50Overlap.all.common.sort.clean.m1.bed


FILEL="LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI14_rosina_HW M4-HW LB_69"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.5thHW_MELP_50Overlap.common
fi
done
done

cat *5thHW_MELP_50Overlap.common > 5thHW_MELP_50Overlap.all.common
sort -k1,1 -k2,2n 5thHW_MELP_50Overlap.all.common > 5thHW_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' 5thHW_MELP_50Overlap.all.common.sort > 5thHW_MELP_50Overlap.all.common.sort.clean
bedtools merge -i 5thHW_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > 5thHW_MELP_50Overlap.all.common.sort.clean.m1.bed


FILEL="LI22_melpomene_FD LI23_melpomene_FD LI6_melpomene_FD LB_72 LB_77 LI13_rosina_FD"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1FD_MELP_50Overlap.common
fi
done
done

cat *D1FD_MELP_50Overlap.common > D1FD_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D1FD_MELP_50Overlap.all.common > D1FD_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1FD_MELP_50Overlap.all.common.sort > D1FD_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D1FD_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FD_MELP_50Overlap.all.common.sort.clean.m1.bed



FILEL="LI22_melpomene_FM LI23_melpomene_FM LI6_melpomene_FM LB_71 LB_76 LI13_rosina_FM"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1FM_MELP_50Overlap.common
fi
done
done

cat *D1FM_MELP_50Overlap.common > D1FM_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D1FM_MELP_50Overlap.all.common > D1FM_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1FM_MELP_50Overlap.all.common.sort > D1FM_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D1FM_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FM_MELP_50Overlap.all.common.sort.clean.m1.bed




FILEL="LI23_melpomene_FP LI6_melpomene_FP LB_70 LI13_rosina_FP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1FP_MELP_50Overlap.common
fi
done
done

cat *D1FP_MELP_50Overlap.common > D1FP_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D1FP_MELP_50Overlap.all.common > D1FP_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1FP_MELP_50Overlap.all.common.sort > D1FP_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D1FP_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1FP_MELP_50Overlap.all.common.sort.clean.m1.bed



FILEL="LI22_melpomene_HA LI23_melpomene_HA LI6_melpomene_HA LI22_melpomene_HP LI23_melpomene_HP LI6_melpomene_HP LB_78 LI13_rosina_HA LB_79 LI13_rosina_HP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D1HW_MELP_50Overlap.common
fi
done
done

cat *D1HW_MELP_50Overlap.common > D1HW_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D1HW_MELP_50Overlap.all.common > D1HW_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D1HW_MELP_50Overlap.all.common.sort > D1HW_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D1HW_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D1HW_MELP_50Overlap.all.common.sort.clean.m1.bed



FILEL="LI25_melpomene_FD LI28_melpomene_FD LB_82 LI19_rosina_FD"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2FD_MELP_50Overlap.common
fi
done
done

cat *D2FD_MELP_50Overlap.common > D2FD_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D2FD_MELP_50Overlap.all.common > D2FD_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2FD_MELP_50Overlap.all.common.sort > D2FD_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D2FD_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FD_MELP_50Overlap.all.common.sort.clean.m1.bed


FILEL="LI25_melpomene_FM LI28_melpomene_FM LB_81 LI19_rosina_FM2"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2FM_MELP_50Overlap.common
fi
done
done

cat *D2FM_MELP_50Overlap.common > D2FM_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D2FM_MELP_50Overlap.all.common > D2FM_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2FM_MELP_50Overlap.all.common.sort > D2FM_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D2FM_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FM_MELP_50Overlap.all.common.sort.clean.m1.bed


FILEL="LI25_melpomene_FP LI28_melpomene_FP LB_80 LI19_rosina_FP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2FP_MELP_50Overlap.common
fi
done
done

cat *D2FP_MELP_50Overlap.common > D2FP_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D2FP_MELP_50Overlap.all.common > D2FP_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2FP_MELP_50Overlap.all.common.sort > D2FP_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D2FP_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2FP_MELP_50Overlap.all.common.sort.clean.m1.bed



FILEL="LI25_melpomene_HA LI28_melpomene_HA LI25_melpomene_HP LI28_melpomene_HP LB_83 LI19_rosina_HA LB_84 LI19_rosina_HP"

for file1 in ${FILEL}
do 
for file2 in ${FILEL}
do 
if [ $file1 != $file2 ] 
then intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $file1.$file2.D2HW_MELP_50Overlap.common
fi
done
done

cat *D2HW_MELP_50Overlap.common > D2HW_MELP_50Overlap.all.common
sort -k1,1 -k2,2n D2HW_MELP_50Overlap.all.common > D2HW_MELP_50Overlap.all.common.sort
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' D2HW_MELP_50Overlap.all.common.sort > D2HW_MELP_50Overlap.all.common.sort.clean
bedtools merge -i D2HW_MELP_50Overlap.all.common.sort.clean -c 1 -o count | awk '$4>0' | cut -d$'\t' -f 1-3 > D2HW_MELP_50Overlap.all.common.sort.clean.m1.bed
