#!/bin/bash
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=1:00:00
#SBATCH --job-name=bed
#SBATCH --error=bed
#SBATCH --output=bed
#SBATCH --partition=rpapaplus
#SBATCH --ntasks=1

module load seq-seq-pan/current
module load python2 

ID=$((SLURM_ARRAY_TASK_ID -1))

samples1=(5thFW_ERA_1bp.m1.sort.bed \
5thHW_ERA_1bp.m1.sort.bed \
D1FD_ERA_1bp.m1.sort.bed \
D1FM_ERA_1bp.m1.sort.bed \
D1FP_ERA_1bp.m1.sort.bed \
D1HW_ERA_1bp.m1.sort.bed \
D2FD_ERA_1bp.m1.sort.bed \
D2FM_ERA_1bp.m1.sort.bed \
D2FP_ERA_1bp.m1.sort.bed \
D2HW_ERA_1bp.m1.sort.bed \
F0.1.2.3.4.5.5thHWERA.common \
F0.1.2.3.4.5.6.7.5thFWERA.common \
F0.1.2.3.4.5.6.D1FPERA.common \
F0.1.2.3.4.5.D1FDERA.common \
F0.1.2.3.4.5.D1FMERA.common \
F0.1.2.3.4.5.6.7.8.9.10.11.12.D1HWERA.common \
F0.1.2.3.4.5.D2FDERA.common \
F0.1.2.3.4.5.D2FMERA.common \
F0.1.2.3.4.5.D2FPERA.common \
F0.1.2.3.4.5.6.7.8.9.10.11.D2HWERA.common)


samples2=(5thFW_MELP_1bp.m1.sort.bed \
5thHW_MELP_1bp.m1.sort.bed \
D1FD_MELP_1bp.m1.sort.bed \
D1FM_MELP_1bp.m1.sort.bed \
D1FP_MELP_1bp.m1.sort.bed \
D1HW_MELP_1bp.m1.sort.bed \
D2FD_MELP_1bp.m1.sort.bed \
D2FM_MELP_1bp.m1.sort.bed \
D2FP_MELP_1bp.m1.sort.bed \
D2HW_MELP_1bp.m1.sort.bed \
F0.1.2.3.4.5.5thFWMELP.common \
F0.1.2.3.4.5.5thHWMELP.common \
F0.1.2.3.4.5.D1FDMELP.common \
F0.1.2.3.4.5.D1FMMELP.common \
F0.1.2.3.D1FPMELP.common \
F0.1.2.3.4.5.6.7.8.9.D1HWMELP.common \
F0.1.2.3.D2FDMELP.common \
F0.1.2.3.D2FMMELP.common \
F0.1.2.3.D2FPMELP.common \
F0.1.2.3.4.5.6.7.D2HWMELP.common)

cd ~/share/MACS2_out_common 

# make bedfile from counts file
cut -d$'\t' -f 1-3 $(echo "${samples[ID]}") > $(echo "${samples[ID]}").BED

# transform coordinates to genome coordinates
python /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/se-seq-pan_bedgraph_chrompos.py -I $(echo "${samples[ID]}").BED -g 1 -p /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/genome_pos_H_e_dem.txt -o start_end_MACS2/$(echo "${samples[ID]}")
#python /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/se-seq-pan_bedgraph_chrompos.py -I $(echo "${samples[ID]}").BED -g 2 -p /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/genome_pos_H_m_melp.txt -o start_end_MACS2/$(echo "${samples[ID]}")

# map coords
nohup seq-seq-pan map -c /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/SeqSeqPan_erato_melp_consensus.fasta -p ./ -n start_end_MACS2/$(echo "${samples[ID]}").BED_start_PAN -i start_end_MACS2/$(echo "${samples[ID]}").BED_start.txt
nohup seq-seq-pan map -c /work/rpapa/sbelleghem/SeqSeqPan_erato_melp/SeqSeqPan_erato_melp_consensus.fasta -p ./ -n start_end_MACS2/$(echo "${samples[ID]}").BED_end_PAN -i start_end_MACS2/$(echo "${samples[ID]}").BED_end.txt

#combine
lineN=$(< start_end_MACS2/$(echo "${samples[ID]}").BED_start.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(cut -d$'\t' -f 2 start_end_MACS2/$(echo "${samples[ID]}").BED_start_PAN.txt) <(cut -d$'\t' -f 2 start_end_MACS2/$(echo "${samples[ID]}").BED_end_PAN.txt) > start_end_MACS2/$(echo "${samples[ID]}").BED_start_end_PAN.txt

# sort rows
lineN=$(< start_end_MACS2/$(echo "${samples[ID]}").BED_start.txt -l)
pr -mts <(yes "pan" | head -n $lineN) <(awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s\t", a[i] ); printf( "\n" ); }' <(cut -d$'\t' -f 2,3 start_end_MACS2/$(echo "${samples[ID]}").BED_start_end_PAN.txt)) > start_end_MACS2/$(echo "${samples[ID]}").BED_start_end_PAN.rowsort.txt
sed -i 's/[[:space:]]*$//' start_end_MACS2/$(echo "${samples[ID]}").BED_start_end_PAN.rowsort.txt

awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' start_end_MACS2/$(echo "${samples[ID]}").BED_start_end_PAN.rowsort.txt > start_end_MACS2/$(echo "${samples[ID]}").BED_start_end_PAN.rowsort.clean.txt

for e in {0..19}
do intersectBed -a $(echo "${samples1[$e]}").BED_start_end_PAN.rowsort.clean.txt -b $(echo "${samples2[$e]}").BED_start_end_PAN.rowsort.clean.txt -f 0.5 -r | wc -l
done


cat 5thFW_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
5thHW_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FD_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FM_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FP_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1HW_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FD_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FM_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FP_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2HW_ERA_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt > ALL_ERA_1bp_pan.bed

cat F0.1.2.3.4.5.5thHWERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.5thFWERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.D1FPERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FDERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FMERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.8.9.10.11.12.D1HWERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D2FDERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D2FMERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D2FPERA.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.8.9.10.11.D2HWERA.common.BED_start_end_PAN.rowsort.clean.txt > ALL_ERA_50_pan.bed

cat 5thFW_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
5thHW_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FD_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FM_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FP_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D1HW_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FD_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FM_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FP_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt \
D2HW_MELP_1bp.m1.sort.bed.BED_start_end_PAN.rowsort.clean.txt > ALL_MELP_1bp_pan.bed

cat F0.1.2.3.4.5.5thFWMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.5thHWMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FDMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FMMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D1FPMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.8.9.D1HWMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D2FDMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D2FMMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D2FPMELP.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.D2HWMELP.common.BED_start_end_PAN.rowsort.clean.txt > ALL_MELP_50_pan.bed

sort -k1,1 -k2,2n ALL_ERA_1bp_pan.bed > ALL_ERA_1bp_pan.sort.bed
sort -k1,1 -k2,2n ALL_ERA_50_pan.bed > ALL_ERA_50_pan.sort.bed
sort -k1,1 -k2,2n ALL_MELP_1bp_pan.bed > ALL_MELP_1bp_pan.sort.bed
sort -k1,1 -k2,2n ALL_MELP_50_pan.bed > ALL_MELP_50_pan.sort.bed

bedtools merge -i ALL_ERA_1bp_pan.sort.bed > ALL_ERA_1bp_pan.merged.sort.bed
bedtools merge -i ALL_ERA_50_pan.sort.bed > ALL_ERA_50_pan.merged.sort.bed
bedtools merge -i ALL_MELP_1bp_pan.sort.bed > ALL_MELP_1bp_pan.merged.sort.bed
bedtools merge -i ALL_MELP_50_pan.sort.bed > ALL_MELP_50_pan.merged.sort.bed

awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_ERA_1bp_pan.merged.sort.bed > ALL_ERA_1bp_pan.merged.sort.clean.bed
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_ERA_50_pan.merged.sort.bed > ALL_ERA_50_pan.merged.sort.clean.bed
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_MELP_1bp_pan.merged.sort.bed > ALL_MELP_1bp_pan.merged.sort.clean.bed
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_MELP_50_pan.merged.sort.bed > ALL_MELP_50_pan.merged.sort.clean.bed

wc -l ALL_ERA_1bp_pan.merged.sort.clean.bed
wc -l ALL_ERA_50_pan.merged.sort.clean.bed
wc -l ALL_MELP_1bp_pan.merged.sort.clean.bed
wc -l ALL_MELP_50_pan.merged.sort.clean.bed

intersectBed -a ALL_ERA_1bp_pan.merged.sort.clean.bed -b ALL_MELP_1bp_pan.merged.sort.clean.bed | wc -l
intersectBed -a ALL_ERA_1bp_pan.merged.sort.clean.bed -b ALL_MELP_1bp_pan.merged.sort.clean.bed -f 0.5 -r | wc -l
intersectBed -a ALL_ERA_50_pan.merged.sort.clean.bed -b ALL_MELP_50_pan.merged.sort.clean.bed | wc -l
intersectBed -a ALL_ERA_50_pan.merged.sort.clean.bed -b ALL_MELP_50_pan.merged.sort.clean.bed -f 0.5 -r | wc -l



cat 5thFW_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
5thHW_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FD_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FM_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1HW_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FD_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FM_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2HW_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt > ALL_ERA_2samp50bp_pan.bed

cat F0.1.2.3.4.5.5thHWERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.5thFWERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.8.9.10.11.12.D1HWERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.8.9.10.11.D2HWERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.D1FPERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FDERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FMERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D2FDERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D2FMERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D2FPERA.1bp.common.BED_start_end_PAN.rowsort.clean.txt > ALL_ERA_50bpsamp1bp_pan.bed


cat 5thFW_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
5thHW_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FD_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FM_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1FP_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D1HW_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FD_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FM_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2FP_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt \
D2HW_MELP_50Overlap.all.common.sort.clean.m1.bed.BED_start_end_PAN.rowsort.clean.txt > ALL_MELP_2samp50bp_pan.bed

cat F0.1.2.3.4.5.5thFWMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.5thHWMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.8.9.D1HWMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.6.7.D2HWMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FDMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.4.5.D1FMMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D1FPMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D2FDMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D2FMMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt \
F0.1.2.3.D2FPMELP.1bp.common.BED_start_end_PAN.rowsort.clean.txt > ALL_MELP_50bpsamp1bp_pan.bed

sort -k1,1 -k2,2n ALL_ERA_2samp50bp_pan.bed > ALL_ERA_2samp50bp_pan.sort.bed
sort -k1,1 -k2,2n ALL_ERA_50bpsamp1bp_pan.bed > ALL_ERA_50bpsamp1bp_pan.sort.bed
sort -k1,1 -k2,2n ALL_MELP_2samp50bp_pan.bed > ALL_MELP_2samp50bp_pan.sort.bed
sort -k1,1 -k2,2n ALL_MELP_50bpsamp1bp_pan.bed > ALL_MELP_50bpsamp1bp_pan.sort.bed

bedtools merge -i ALL_ERA_2samp50bp_pan.sort.bed > ALL_ERA_2samp50bp_pan.merged.sort.bed
bedtools merge -i ALL_ERA_50bpsamp1bp_pan.sort.bed > ALL_ERA_50bpsamp1bp_pan.merged.sort.bed
bedtools merge -i ALL_MELP_2samp50bp_pan.sort.bed > ALL_MELP_2samp50bp_pan.merged.sort.bed
bedtools merge -i ALL_MELP_50bpsamp1bp_pan.sort.bed > ALL_MELP_50bpsamp1bp_pan.merged.sort.bed

awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_ERA_2samp50bp_pan.merged.sort.bed > ALL_ERA_2samp50bp_pan.merged.sort.clean.bed
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_ERA_50bpsamp1bp_pan.merged.sort.bed > ALL_ERA_50bpsamp1bp_pan.merged.sort.clean.bed
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_MELP_2samp50bp_pan.merged.sort.bed > ALL_MELP_2samp50bp_pan.merged.sort.clean.bed
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' ALL_MELP_50bpsamp1bp_pan.merged.sort.bed > ALL_MELP_50bpsamp1bp_pan.merged.sort.clean.bed

wc -l ALL_ERA_2samp50bp_pan.merged.sort.clean.bed
wc -l ALL_ERA_50bpsamp1bp_pan.merged.sort.clean.bed
wc -l ALL_MELP_2samp50bp_pan.merged.sort.clean.bed
wc -l ALL_MELP_50bpsamp1bp_pan.merged.sort.clean.bed

intersectBed -a ALL_ERA_2samp50bp_pan.merged.sort.clean.bed -b ALL_MELP_2samp50bp_pan.merged.sort.clean.bed | wc -l
intersectBed -a ALL_ERA_2samp50bp_pan.merged.sort.clean.bed -b ALL_MELP_2samp50bp_pan.merged.sort.clean.bed -f 0.5 -r | wc -l
intersectBed -a ALL_ERA_50bpsamp1bp_pan.merged.sort.clean.bed -b ALL_MELP_50bpsamp1bp_pan.merged.sort.clean.bed | wc -l
intersectBed -a ALL_ERA_50bpsamp1bp_pan.merged.sort.clean.bed -b ALL_MELP_50bpsamp1bp_pan.merged.sort.clean.bed -f 0.5 -r | wc -l

for e in {0..19}
do intersectBed -a $(echo "${samples1[$e]}").BED_start_end_PAN.rowsort.clean.txt -b $(echo "${samples2[$e]}").BED_start_end_PAN.rowsort.clean.txt -f 0.5 -r | wc -l
done
