#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=4:00:00
#SBATCH --job-name=bed
#SBATCH --error=bed
#SBATCH --output=bed
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1


module load bedtools
module load samtools

#83
sampleList1=(BR10_Demophoon_Brain BR14_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain E4-Head LB_3 LB_45 LB_50 LI2_demophoon_FD LI21_demophoon_FD LI27_Demophoon_FD LB_2 LB_44 LB_49 LI2_demophoon_FM LI21_demophoon_FM LI27_Demophoon_FP E3_FW E4-FW FW-pboy LB_41 LI7_demophoon_FW LB_1 LB_43 LB_48 LI2_demophoon_FP LI21_demophoon_FP LI27_Demophoon_FM LB_4 LB_46 LB_51 LI2_demophoon_HA LI21_demophoon_HA LI27_Demophoon_HA E3_HW LB_42 LI7_demophoon_HW LB_47 LB_5 LB_52 LI2_demophoon_HP LI20_demophoon_HP LI21_demophoon_HP LI27_Demophoon_HP BR2_hydara_brain BR3_hydara_brain BR6_Hydara_Brain LB_24 LB_33 LB_38 LB_65 LI1_hydara_FD LI4_hydara_FD LB_23 LB_32 LB_37 LB_64 LI1_hydara_FM LI4_hydara_FM LB_20 LB_27 LB_29 LB_22 LB_31 LB_36 LB_63 LI1_hydara_FP LI4_hydara_FP LB_25 LB_34 LB_39 LB_66 LI1_hydara_HA LI4_hydara_HA LB_21 LB_28 LB_30 LB_26 LB_35 LB_40 LB_67 LI1_hydara_HP LI4_hydara_HP)

## make sample list command
ALL_LIST=""
for FILE in ${sampleList1[*]}
do
ALL_LIST="$ALL_LIST /work/rpapa/share/BAM_ATACseq_trimmomatic/${FILE}_Herato.trim.filtered.sorted.nd.bam"
done
echo $ALL_LIST

eval command=\$$(echo ALL_LIST)

#195
indID=$((SLURM_ARRAY_TASK_ID -1))
scafs=(Herato0101 Herato0201 Herato0202 Herato0203 Herato0204 Herato0205 Herato0206 Herato0207 Herato0208 Herato0209 Herato0210 Herato0211 Herato0212 Herato0213 Herato0214 Herato0215 Herato0301 Herato0302 Herato0303 Herato0304 Herato0305 Herato0306 Herato0307 Herato0308 Herato0309 Herato0310 Herato0401 Herato0402 Herato0403 Herato0404 Herato0405 Herato0406 Herato0407 Herato0408 Herato0409 Herato0410 Herato0411 Herato0412 Herato0413 Herato0414 Herato0415 Herato0416 Herato0417 Herato0418 Herato0419 Herato0501 Herato0502 Herato0503 Herato0504 Herato0505 Herato0506 Herato0507 Herato0508 Herato0509 Herato0510 Herato0511 Herato0601 Herato0602 Herato0603 Herato0604 Herato0605 Herato0606 Herato0607 Herato0608 Herato0609 Herato0701 Herato0801 Herato0802 Herato0803 Herato0804 Herato0805 Herato0806 Herato0807 Herato0808 Herato0809 Herato0810 Herato0811 Herato0812 Herato0813 Herato0814 Herato0815 Herato0816 Herato0817 Herato0818 Herato0819 Herato0820 Herato0821 Herato0901 Herato0902 Herato0903 Herato0904 Herato1001 Herato1002 Herato1003 Herato1004 Herato1005 Herato1006 Herato1007 Herato1101 Herato1102 Herato1103 Herato1104 Herato1105 Herato1106 Herato1107 Herato1108 Herato1109 Herato1110 Herato1111 Herato1112 Herato1113 Herato1114 Herato1115 Herato1116 Herato1201 Herato1202 Herato1301 Herato1401 Herato1402 Herato1403 Herato1404 Herato1405 Herato1406 Herato1407 Herato1408 Herato1409 Herato1410 Herato1411 Herato1501 Herato1502 Herato1503 Herato1504 Herato1505 Herato1506 Herato1507 Herato1508 Herato1509 Herato1510 Herato1511 Herato1512 Herato1513 Herato1514 Herato1515 Herato1516 Herato1517 Herato1518 Herato1519 Herato1520 Herato1521 Herato1522 Herato1523 Herato1524 Herato1601 Herato1602 Herato1603 Herato1604 Herato1605 Herato1701 Herato1702 Herato1703 Herato1704 Herato1705 Herato1706 Herato1707 Herato1708 Herato1709 Herato1710 Herato1711 Herato1712 Herato1713 Herato1714 Herato1715 Herato1716 Herato1717 Herato1718 Herato1719 Herato1801 Herato1802 Herato1803 Herato1804 Herato1805 Herato1806 Herato1807 Herato1901 Herato1902 Herato1903 Herato1904 Herato1905 Herato1906 Herato1907 Herato1908 Herato1909 Herato1910 Herato2001 Herato2101)
# count reads in peaks
bedtools multicov -bams $(echo $command) -bed MACS2_all_combined/bed_scafs_erato/dem_hyd_MACS2_peaks.$(echo "${scafs[indID]}").merged.sort.bed > MACS2_all_combined/bed_counts_erato_dem_hyd/dem_hyd_peaks_MACS2.$(echo "${scafs[indID]}").counts

