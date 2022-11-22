######
# Generate peak sets in H. erato for each tissue/time 1bp overlap in all samples
######

FILEL=(E4-FW FW-pboy LI7_demophoon_FW LB_41 LB_20 LB_27 LB_29 E3_FW)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.5thFWERA.1bp.common
            echo $FILENAME.5thFWERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.5thFWERA.1bp.common
            intersectBed -a $FILENAME.5thFWERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.5thFWERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(E3_HW LB_42 LI7_demophoon_HW LB_21 LB_28 LB_30)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.5thHWERA.1bp.common
            echo $FILENAME.5thHWERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.5thHWERA.1bp.common
            intersectBed -a $FILENAME.5thHWERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.5thHWERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_45 LI21_demophoon_FD LI27_Demophoon_FD LB_24 LB_65 LI4_hydara_FD)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1FDERA.1bp.common
            echo $FILENAME.D1FDERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FDERA.1bp.common
            intersectBed -a $FILENAME.D1FDERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1FDERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_44 LI21_demophoon_FM LI27_Demophoon_FM LB_23 LB_64 LI4_hydara_FM)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1FMERA.1bp.common
            echo $FILENAME.D1FMERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FMERA.1bp.common
            intersectBed -a $FILENAME.D1FMERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1FMERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_43 LI21_demophoon_FP LI27_Demophoon_FP LB_26 LB_22 LB_63 LI4_hydara_FP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1FPERA.1bp.common
            echo $FILENAME.D1FPERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FPERA.1bp.common
            intersectBed -a $FILENAME.D1FPERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1FPERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_46 LI21_demophoon_HA LI27_Demophoon_HA LB_47 LI20_demophoon_HP LI21_demophoon_HP LI27_Demophoon_HP LB_25 LB_66 LI4_hydara_HA LB_26 LB_67 LI4_hydara_HP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1HWERA.1bp.common
            echo $FILENAME.D1HWERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1HWERA.1bp.common
            intersectBed -a $FILENAME.D1HWERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1HWERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_3 LB_50 LI2_demophoon_FD LB_33 LB_38 LI1_hydara_FD)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2FDERA.1bp.common
            echo $FILENAME.D2FDERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2FDERA.1bp.common
            intersectBed -a $FILENAME.D2FDERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2FDERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_2 LB_49 LI2_demophoon_FM LB_32 LB_37 LI1_hydara_FM)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2FMERA.1bp.common
            echo $FILENAME.D2FMERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2FMERA.1bp.common
            intersectBed -a $FILENAME.D2FMERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2FMERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_1 LB_48 LI2_demophoon_FP LB_31 LB_36 LI1_hydara_FP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2FPERA.1bp.common
            echo $FILENAME.D2FPERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2FPERA.1bp.common
            intersectBed -a $FILENAME.D2FPERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2FPERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LB_4 LB_51 LI2_demophoon_HA LB_5 LB_52 LI2_demophoon_HP LB_34 LB_39 LI1_hydara_HA LB_35 LB_40 LI1_hydara_HP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2HWERA.1bp.common
            echo $FILENAME.D2HWERA.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2HWERA.1bp.common
            intersectBed -a $FILENAME.D2HWERA.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2HWERA.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done




len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $FILENAME.D1FPERA.common
            echo $FILENAME.D1FPERA.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FPERA.common
            intersectBed -a $FILENAME.D1FPERA.common -b $file2"_peaks.narrowPeak" -f 0.5 -r > $FILENAME2.D1FPERA.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


######
# Generate peak sets in H. melpomene for each tissue/time 1bp overlap in all samples
######


FILEL=(LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI14_rosina_HW M4-HW LB_68)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.5thFWMELP.1bp.common
            echo $FILENAME.5thFWMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.5thFWMELP.1bp.common
            intersectBed -a $FILENAME.5thFWMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.5thFWMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done



FILEL=(LI12_melpomene_HW LI15_melpomene_HW LI8_melpomene_HW LI14_rosina_HW M4-HW LB_69)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.5thHWMELP.1bp.common
            echo $FILENAME.5thHWMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.5thHWMELP.1bp.common
            intersectBed -a $FILENAME.5thHWMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.5thHWMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done



FILEL=(LI22_melpomene_FD LI23_melpomene_FD LI6_melpomene_FD LB_72 LB_77 LI13_rosina_FD)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1FDMELP.1bp.common
            echo $FILENAME.D1FDMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FDMELP.1bp.common
            intersectBed -a $FILENAME.D1FDMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1FDMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done



FILEL=(LI22_melpomene_FM LI23_melpomene_FM LI6_melpomene_FM LB_71 LB_76 LI13_rosina_FM)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1FMMELP.1bp.common
            echo $FILENAME.D1FMMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FMMELP.1bp.common
            intersectBed -a $FILENAME.D1FMMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1FMMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done



FILEL=(LI23_melpomene_FP LI6_melpomene_FP LB_70 LI13_rosina_FP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1FPMELP.1bp.common
            echo $FILENAME.D1FPMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1FPMELP.1bp.common
            intersectBed -a $FILENAME.D1FPMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1FPMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LI22_melpomene_HA LI23_melpomene_HA LI6_melpomene_HA LI22_melpomene_HP LI23_melpomene_HP LI6_melpomene_HP LB_78 LI13_rosina_HA LB_79 LI13_rosina_HP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D1HWMELP.1bp.common
            echo $FILENAME.D1HWMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D1HWMELP.1bp.common
            intersectBed -a $FILENAME.D1HWMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D1HWMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LI25_melpomene_FD LI28_melpomene_FD LB_82 LI19_rosina_FD)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2FDMELP.1bp.common
            echo $FILENAME.D2FDMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2FDMELP.1bp.common
            intersectBed -a $FILENAME.D2FDMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2FDMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LI25_melpomene_FM LI28_melpomene_FM LB_81 LI19_rosina_FM2)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2FMMELP.1bp.common
            echo $FILENAME.D2FMMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2FMMELP.1bp.common
            intersectBed -a $FILENAME.D2FMMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2FMMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LI25_melpomene_FP LI28_melpomene_FP LB_80 LI19_rosina_FP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2FPMELP.1bp.common
            echo $FILENAME.D2FPMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2FPMELP.1bp.common
            intersectBed -a $FILENAME.D2FPMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2FPMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


FILEL=(LI25_melpomene_HA LI28_melpomene_HA LI25_melpomene_HP LI28_melpomene_HP LB_83 LI19_rosina_HA LB_84 LI19_rosina_HP)

len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" > $FILENAME.D2HWMELP.1bp.common
            echo $FILENAME.D2HWMELP.1bp.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2HWMELP.1bp.common
            intersectBed -a $FILENAME.D2HWMELP.1bp.common -b $file2"_peaks.narrowPeak" > $FILENAME2.D2HWMELP.1bp.common
            FILENAME=$FILENAME.$((i+1))
      fi
done




len=${#FILEL[@]}-1

FILENAME=F0.1

for (( i=0; i<$len; i++ ))
do
   file1=${FILEL[$i]}
   file2=${FILEL[$i+1]}

   if [ $i == 0 ]
         then
            intersectBed -a $file1"_peaks.narrowPeak" -b $file2"_peaks.narrowPeak" -f 0.5 -r > $FILENAME.D2HWMELP.common
            echo $FILENAME.D2HWMELP.common
         else
            FILENAME2=$FILENAME.$((i+1))
            echo $FILENAME.D2HWMELP.common
            intersectBed -a $FILENAME.D2HWMELP.common -b $file2"_peaks.narrowPeak" -f 0.5 -r > $FILENAME2.D2HWMELP.common
            FILENAME=$FILENAME.$((i+1))
      fi
done


