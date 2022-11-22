#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=8:00:00
#SBATCH --job-name=scale
#SBATCH --error=scale
#SBATCH --output=scale
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1


indID=$((SLURM_ARRAY_TASK_ID -1))

#83
#samples=(BR10_Demophoon_Brain BR14_Demophoon_Brain BR12_Demophoon_Brain BR15_Demophoon_Brain E4-Head LB_3 LB_45 LB_50 LI2_demophoon_FD LI21_demophoon_FD LI27_Demophoon_FD LB_2 LB_44 LB_49 LI2_demophoon_FM LI21_demophoon_FM LI27_Demophoon_FP E3_FW E4-FW FW-pboy LB_41 LI7_demophoon_FW LB_1 LB_43 LB_48 LI2_demophoon_FP LI21_demophoon_FP LI27_Demophoon_FM LB_4 LB_46 LB_51 LI2_demophoon_HA LI21_demophoon_HA LI27_Demophoon_HA E3_HW LB_42 LI7_demophoon_HW LB_47 LB_5 LB_52 LI2_demophoon_HP LI20_demophoon_HP LI21_demophoon_HP LI27_Demophoon_HP BR2_hydara_brain BR3_hydara_brain BR6_Hydara_Brain LB_24 LB_33 LB_38 LB_65 LI1_hydara_FD LI4_hydara_FD LB_23 LB_32 LB_37 LB_64 LI1_hydara_FM LI4_hydara_FM LB_20 LB_27 LB_29 LB_22 LB_31 LB_36 LB_63 LI1_hydara_FP LI4_hydara_FP LB_25 LB_34 LB_39 LB_66 LI1_hydara_HA LI4_hydara_HA LB_21 LB_28 LB_30 LB_26 LB_35 LB_40 LB_67 LI1_hydara_HP LI4_hydara_HP)
#scale=(1.827449083 0.314679357 1.140134058 0.993308977 0.528384968 0.766024856 0.950165287 0.591433052 1.51478408 0.730195322 2.210744815 0.705723713 0.740203786 1.472849785 1.118127352 0.645861735 1.077976515 2.303345967 1.016188922 0.408890344 1.059150245 1.252479246 0.84177101 0.800494165 0.810974971 1.623929731 0.624383429 1.366700108 2.348512899 0.934940952 0.692317208 1.786358234 0.591443093 2.354779038 0.992913794 0.927246174 1.344028224 1.19177926 1.894605452 0.511703774 2.272828089 0.497169626 0.62612469 1.289661708 0.553677245 0.85892251 0.631878559 0.528164722 1.486235032 1.24978839 0.827187864 1.637900268 1.207989044 0.987601487 0.681264836 0.892173381 0.834348586 0.869769488 0.861203369 0.334190063 1.454375434 1.908115471 0.591222726 0.6844118 1.403119961 1.060781505 0.786715512 1.074699795 1.010399115 0.639424371 1.198196393 0.761098448 4.31597668 0.961801384 0.91798891 2.214649707 0.970554925 0.554721793 0.882950252 1.936802963 0.714278807 1.308109887 1.048325009)

#60
#samples=(BR11_Rosina_Brain BR5_Rosina_Brain M4-Head LB_72 LB_77 LB_82 LI13_rosina_FD LI19_rosina_FD LB_71 LB_76 LB_81 LI13_rosina_FM LI19_rosina_FM2 LB_68 LI14_rosina_FW M4-FW LB_70 LB_80 LI13_rosina_FP LI19_rosina_FP LB_78 LB_83 LI13_rosina_HA LI19_rosina_HA LB_69 LI14_rosina_HW M4-HW LB_79 LB_84 LI13_rosina_HP LI19_rosina_HP LI22_melpomene_FD LI23_melpomene_FD LI25_melpomene_FD LI28_melpomene_FD LI6_melpomene_FD LI22_melpomene_FM LI23_melpomene_FM LI25_melpomene_FM LI28_melpomene_FM LI6_melpomene_FM LI12_melpomene_FW LI15_melpomene_FW LI8_melpomene_FW LI23_melpomene_FP LI25_melpomene_FP LI28_melpomene_FP LI6_melpomene_FP LI22_melpomene_HA LI23_melpomene_HA LI25_melpomene_HA LI28_melpomene_HA LI6_melpomene_HA LI12_melpomene_HW LI8_melpomene_HW LI22_melpomene_HP LI23_melpomene_HP LI25_melpomene_HP LI28_melpomene_HP LI6_melpomene_HP)
#scale=(2.580556658 3.278080306 0.42508338 2.124077436 1.555910025 1.639068744 0.731944057 2.45102513 2.017488131 1.7463467 0.952374321 0.892778494 4.725132783 0.949031132 1.222057172 1.11712052 1.400251222 1.264178374 0.447203845 0.868883025 1.02286447 1.432263886 0.652184459 0.815850352 0.892511769 1.094311474 0.59636932 2.883359238 0.868132491 0.563576857 0.979013323 0.680154144 0.801245068 0.687675277 1.020552472 0.733450686 0.643839022 0.605590839 0.721764599 1.031611967 0.508385411 1.279063655 0.871185578 1.46227805 0.69301579 0.895407549 1.277431112 0.498741119 0.679758296 0.697494275 0.915126369 1.90495131 0.738719833 0.994727366 0.836299777 0.784027235 0.598834918 0.620406046 1.019517655 0.698017713)


#module load bedtools
#bedtools genomecov -ibam /work/rpapa/share/BAM_ATACseq_trimmomatic/$(echo "${samples[indID]}")_Hmel2.trim.filtered.sorted.nd.bam -bg -scale $(echo "${scale[indID]}") -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes > /work/rpapa/share/BAM_ATACseq_trimmomatic_scaled_bg/$(echo "${samples[indID]}")_Hmel2_normalised.bg

module load wiggletools
cd /work/rpapa/share/BAM_ATACseq_trimmomatic_scaled_bg

nohup wiggletools write_bg brain_5th_H_erato_normalized_mean.bg mean BR10_Demophoon_Brain_Herato_normalised.bg BR14_Demophoon_Brain_Herato_normalised.bg BR12_Demophoon_Brain_Herato_normalised.bg BR15_Demophoon_Brain_Herato_normalised.bg E4-Head_Herato_normalised.bg &
nohup wiggletools write_bg FW_5th_H_erato_normalized_mean.bg mean E3_FW_Herato_normalised.bg E4-FW_Herato_normalised.bg FW-pboy_Herato_normalised.bg LB_41_Herato_normalised.bg LI7_demophoon_FW_Herato_normalised.bg &
nohup wiggletools write_bg HW_5th_H_erato_normalized_mean.bg mean E3_HW_Herato_normalised.bg LB_42_Herato_normalised.bg LI7_demophoon_HW_Herato_normalised.bg &

nohup wiggletools write_bg FD_D1_H_erato_normalized_mean.bg mean LB_45_Herato_normalised.bg LI21_demophoon_FD_Herato_normalised.bg LI27_Demophoon_FD_Herato_normalised.bg &
nohup wiggletools write_bg FM_D1_H_erato_normalized_mean.bg mean LB_44_Herato_normalised.bg LI21_demophoon_FM_Herato_normalised.bg LI27_Demophoon_FM_Herato_normalised.bg &
nohup wiggletools write_bg FP_D1_H_erato_normalized_mean.bg mean LI27_Demophoon_FP_Herato_normalised.bg LB_43_Herato_normalised.bg LI21_demophoon_FP_Herato_normalised.bg &
nohup wiggletools write_bg HW_D1_H_erato_normalized_mean.bg mean LB_46_Herato_normalised.bg LI21_demophoon_HA_Herato_normalised.bg LI27_Demophoon_HA_Herato_normalised.bg LB_47_Herato_normalised.bg LI21_demophoon_HP_Herato_normalised.bg LI27_Demophoon_HP_Herato_normalised.bg &

nohup wiggletools write_bg FD_D2_H_erato_normalized_mean.bg mean LB_3_Herato_normalised.bg LB_50_Herato_normalised.bg LI2_demophoon_FD_Herato_normalised.bg &
nohup wiggletools write_bg FM_D2_H_erato_normalized_mean.bg mean LB_49_Herato_normalised.bg LI2_demophoon_FM_Herato_normalised.bg LB_2_Herato_normalised.bg &
nohup wiggletools write_bg FP_D2_H_erato_normalized_mean.bg mean LB_48_Herato_normalised.bg LI2_demophoon_FP_Herato_normalised.bg LB_1_Herato_normalised.bg &
nohup wiggletools write_bg HW_D2_H_erato_normalized_mean.bg mean LB_4_Herato_normalised.bg LB_51_Herato_normalised.bg LI2_demophoon_HA_Herato_normalised.bg LB_5_Herato_normalised.bg LB_52_Herato_normalised.bg LI2_demophoon_HP_Herato_normalised.bg &


nohup wiggletools write_bg brain_5th_H_melp_normalized_mean.bg mean BR11_Rosina_Brain_Hmel2_normalised.bg BR5_Rosina_Brain_Hmel2_normalised.bg M4-Head_Hmel2_normalised.bg &
nohup wiggletools write_bg FW_5th_H_melp_normalized_mean.bg mean LB_68_Hmel2_normalised.bg LI14_rosina_FW_Hmel2_normalised.bg M4-FW_Hmel2_normalised.bg &
nohup wiggletools write_bg HW_5th_H_melp_normalized_mean.bg mean LB_69_Hmel2_normalised.bg LI14_rosina_HW_Hmel2_normalised.bg M4-HW_Hmel2_normalised.bg &

nohup wiggletools write_bg FD_D1_H_melp_normalized_mean.bg mean LB_72_Hmel2_normalised.bg LB_77_Hmel2_normalised.bg LI13_rosina_FD_Hmel2_normalised.bg &
nohup wiggletools write_bg FM_D1_H_melp_normalized_mean.bg mean LB_71_Hmel2_normalised.bg LB_76_Hmel2_normalised.bg LI13_rosina_FM_Hmel2_normalised.bg &
nohup wiggletools write_bg FP_D1_H_melp_normalized_mean.bg mean LB_70_Hmel2_normalised.bg LI13_rosina_FP_Hmel2_normalised.bg &
nohup wiggletools write_bg HW_D1_H_melp_normalized_mean.bg mean LB_78_Hmel2_normalised.bg LI13_rosina_HA_Hmel2_normalised.bg LB_79_Hmel2_normalised.bg LI13_rosina_HP_Hmel2_normalised.bg &

nohup wiggletools write_bg FD_D2_H_melp_normalized_mean.bg mean LB_82_Hmel2_normalised.bg LI19_rosina_FD_Hmel2_normalised.bg &
nohup wiggletools write_bg FM_D2_H_melp_normalized_mean.bg mean LB_81_Hmel2_normalised.bg LI19_rosina_FM2_Hmel2_normalised.bg &
nohup wiggletools write_bg FP_D2_H_melp_normalized_mean.bg mean LB_80_Hmel2_normalised.bg LI19_rosina_FP_Hmel2_normalised.bg &
nohup wiggletools write_bg HW_D2_H_melp_normalized_mean.bg mean LB_83_Hmel2_normalised.bg LI19_rosina_HA_Hmel2_normalised.bg LB_84_Hmel2_normalised.bg LI19_rosina_HP_Hmel2_normalised.bg &



nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b brain_5th_H_erato_normalized_mean.bg -c 4 -o mean > brain_5th_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FW_5th_H_erato_normalized_mean.bg -c 4 -o mean > FW_5th_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b HW_5th_H_erato_normalized_mean.bg -c 4 -o mean > HW_5th_H_erato_normalized_mean.w30s0bin.bg &

nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FD_D1_H_erato_normalized_mean.bg -c 4 -o mean > FD_D1_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FM_D1_H_erato_normalized_mean.bg -c 4 -o mean > FM_D1_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FP_D1_H_erato_normalized_mean.bg -c 4 -o mean > FP_D1_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b HW_D1_H_erato_normalized_mean.bg -c 4 -o mean > HW_D1_H_erato_normalized_mean.w30s0bin.bg &

nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FD_D2_H_erato_normalized_mean.bg -c 4 -o mean > FD_D2_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FM_D2_H_erato_normalized_mean.bg -c 4 -o mean > FM_D2_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b FP_D2_H_erato_normalized_mean.bg -c 4 -o mean > FP_D2_H_erato_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g ~/share/REF/H_erato_dem/Herato_final.fasta.sizes -w 30 | bedtools map -a - -b HW_D2_H_erato_normalized_mean.bg -c 4 -o mean > HW_D2_H_erato_normalized_mean.w30s0bin.bg &


nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b brain_5th_H_melp_normalized_mean.bg -c 4 -o mean > brain_5th_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FW_5th_H_melp_normalized_mean.bg -c 4 -o mean > FW_5th_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b HW_5th_H_melp_normalized_mean.bg -c 4 -o mean > HW_5th_H_melp_normalized_mean.w30s0bin.bg &

nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FD_D1_H_melp_normalized_mean.bg -c 4 -o mean > FD_D1_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FM_D1_H_melp_normalized_mean.bg -c 4 -o mean > FM_D1_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FP_D1_H_melp_normalized_mean.bg -c 4 -o mean > FP_D1_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b HW_D1_H_melp_normalized_mean.bg -c 4 -o mean > HW_D1_H_melp_normalized_mean.w30s0bin.bg &

nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FD_D2_H_melp_normalized_mean.bg -c 4 -o mean > FD_D2_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FM_D2_H_melp_normalized_mean.bg -c 4 -o mean > FM_D2_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b FP_D2_H_melp_normalized_mean.bg -c 4 -o mean > FP_D2_H_melp_normalized_mean.w30s0bin.bg &
nohup bedtools makewindows -g /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes -w 30 | bedtools map -a - -b HW_D2_H_melp_normalized_mean.bg -c 4 -o mean > HW_D2_H_melp_normalized_mean.w30s0bin.bg &


nohup grep -v '\.$' brain_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > brain_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FW_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FW_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' HW_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > HW_5th_H_erato_normalized_mean.w30s0bin.sort.bg &

nohup grep -v '\.$' FD_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FD_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FM_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FM_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FP_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FP_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' HW_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > HW_5th_H_erato_normalized_mean.w30s0bin.sort.bg &

nohup grep -v '\.$' FD_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FD_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FM_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FM_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FP_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FP_5th_H_erato_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' HW_5th_H_erato_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > HW_5th_H_erato_normalized_mean.w30s0bin.sort.bg &


nohup grep -v '\.$' brain_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > brain_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FW_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FW_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' HW_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > HW_5th_H_melp_normalized_mean.w30s0bin.sort.bg &

nohup grep -v '\.$' FD_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FD_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FM_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FM_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FP_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FP_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' HW_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > HW_5th_H_melp_normalized_mean.w30s0bin.sort.bg &

nohup grep -v '\.$' FD_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FD_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FM_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FM_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' FP_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > FP_5th_H_melp_normalized_mean.w30s0bin.sort.bg &
nohup grep -v '\.$' HW_5th_H_melp_normalized_mean.w30s0bin.bg | sort -g -k1,1 -k2,2n > HW_5th_H_melp_normalized_mean.w30s0bin.sort.bg &


nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig brain_5th_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes brain_5th_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FW_5th_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FW_5th_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig HW_5th_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes HW_5th_H_erato_normalized_mean.w30s0bin.bw &

nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FD_D1_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FD_D1_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FM_D1_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FM_D1_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FP_D1_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FP_D1_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig HW_D1_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes HW_D1_H_erato_normalized_mean.w30s0bin.bw &

nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FD_D2_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FD_D2_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FM_D2_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FM_D2_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FP_D2_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes FP_D2_H_erato_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig HW_D2_H_erato_normalized_mean.w30s0bin.sort.bg ~/share/REF/H_erato_dem/Herato_final.fasta.sizes HW_D2_H_erato_normalized_mean.w30s0bin.bw &


nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig brain_5th_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes brain_5th_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FW_5th_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FW_5th_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig HW_5th_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes HW_5th_H_melp_normalized_mean.w30s0bin.bw &

nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FD_D1_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FD_D1_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FM_D1_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FM_D1_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FP_D1_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FP_D1_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig HW_D1_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes HW_D1_H_melp_normalized_mean.w30s0bin.bw &

nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FD_D2_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FD_D2_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FM_D2_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FM_D2_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig FP_D2_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes FP_D2_H_melp_normalized_mean.w30s0bin.bw &
nohup /work/rpapa/sbelleghem/scripts/bedGraphToBigWig HW_D2_H_melp_normalized_mean.w30s0bin.sort.bg /work/rpapa/share/REF/Hmel2/Hmel2.fa.sizes HW_D2_H_melp_normalized_mean.w30s0bin.bw &




