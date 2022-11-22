# make bedfile from counts file
cut -d$'\t' -f 1-3 start_end_MACS2/dem_hyd.MACS2.peaks.all.merged.sort.counts > start_end_MACS2/erato_peaks.bed
cut -d$'\t' -f 1-3 start_end_MACS2/melp_ros.MACS2.peaks.all.merged.sort.counts > start_end_MACS2/melp_peaks.bed

# transform coordinates to genome coordinates
python se-seq-pan_bedgraph_chrompos.py -I start_end_MACS2/erato_peaks.bed -g 1 -p genome_pos_H_e_dem.txt -o start_end_MACS2/erato_peaks
python se-seq-pan_bedgraph_chrompos.py -I start_end_MACS2/melp_peaks.bed -g 2 -p genome_pos_H_m_melp.txt -o start_end_MACS2/melp_peaks

# map coords
module load seq-seq-pan/current
nohup seq-seq-pan map -c SeqSeqPan_erato_melp_consensus.fasta -p ./ -n start_end_MACS2/erato_peaks_start_pan -i start_end_MACS2/erato_peaks_start.txt &
nohup seq-seq-pan map -c SeqSeqPan_erato_melp_consensus.fasta -p ./ -n start_end_MACS2/erato_peaks_end_pan -i start_end_MACS2/erato_peaks_end.txt &

nohup seq-seq-pan map -c SeqSeqPan_erato_melp_consensus.fasta -p ./ -n start_end_MACS2/melp_peaks_start_pan -i start_end_MACS2/melp_peaks_start.txt &
nohup seq-seq-pan map -c SeqSeqPan_erato_melp_consensus.fasta -p ./ -n start_end_MACS2/melp_peaks_end_pan -i start_end_MACS2/melp_peaks_end.txt &

#combine
lineN=$(< start_end_MACS2/erato_peaks_start.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(cut -d$'\t' -f 2 start_end_MACS2/erato_peaks_start_pan.txt) <(cut -d$'\t' -f 2 start_end_MACS2/erato_peaks_end_pan.txt) > start_end_MACS2/erato_peaks_start_end_pan.txt

lineN=$(< start_end_MACS2/melp_peaks_start.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(cut -d$'\t' -f 2 start_end_MACS2/melp_peaks_start_pan.txt) <(cut -d$'\t' -f 2 start_end_MACS2/melp_peaks_end_pan.txt) > start_end_MACS2/melp_peaks_start_end_pan.txt

# sort rows
lineN=$(< start_end_MACS2/erato_peaks_start_end_pan.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s\t", a[i] ); printf( "\n" ); }' <(cut -d$'\t' -f 2,3 start_end_MACS2/erato_peaks_start_end_pan.txt)) > start_end_MACS2/erato_peaks_start_end_pan.rowsort.txt
sed -i 's/[[:space:]]*$//' start_end_MACS2/erato_peaks_start_end_pan.rowsort.txt

lineN=$(< start_end_MACS2/melp_peaks_start_end_pan.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(awk ' {split( $0, a, " " ); asort( a ); for( i = 1; i <= length(a); i++ ) printf( "%s\t", a[i] ); printf( "\n" ); }' <(cut -d$'\t' -f 2,3 start_end_MACS2/melp_peaks_start_end_pan.txt)) > start_end_MACS2/melp_peaks_start_end_pan.rowsort.txt
sed -i 's/[[:space:]]*$//' start_end_MACS2/melp_peaks_start_end_pan.rowsort.txt


awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' erato_peaks_start_end_pan.rowsort.txt > erato_peaks_start_end_pan.rowsort.clean.txt
awk '{ if ($2 < $3 && $3-$2 < 10000 && $2!=$3) print $0}' melp_peaks_start_end_pan.rowsort.txt > melp_peaks_start_end_pan.rowsort.clean.txt

bedtools intersect -r -f 0.5 -wao -a erato_peaks_start_end_pan.rowsort.clean.txt -b melp_peaks_start_end_pan.rowsort.clean.txt > match_melp_to_erato.0.5.txt
bedtools intersect -r -f 0.5 -wao -b erato_peaks_start_end_pan.rowsort.clean.txt -a melp_peaks_start_end_pan.rowsort.clean.txt > match_erato_to_melp.0.5.txt

bedtools intersect -b erato_peaks_start_end_pan.rowsort.clean.txt -a melp_peaks_start_end_pan.rowsort.clean.txt > intersect_erato_melp.txt
bedtools intersect -f 0.5 -r -b erato_peaks_start_end_pan.rowsort.clean.txt -a melp_peaks_start_end_pan.rowsort.clean.txt > intersect_erato_melp.0.5.txt

# unique peaks
grep '\-1' match_melp_to_erato.0.5.txt > unique_peaks_erato.0.5.txt
grep '\-1' match_erato_to_melp.0.5.txt > unique_peaks_melp.0.5.txt

cut -d$'\t' -f 1,2,3 unique_peaks_erato.txt > unique_peaks_erato.bed
cut -d$'\t' -f 1,2,3 unique_peaks_melp.txt > unique_peaks_melp.bed

# intersect unique peaks with missing
bedtools intersect -f 0.5 -r -a start_end_MACS2/unique_peaks_melp.bed -b 1_missing_intervals.bed > start_end_MACS2/unique_peaks_melp_missing_in_erato_0.5.txt
bedtools intersect -f 0.5 -r -a start_end_MACS2/unique_peaks_erato.bed -b 2_missing_intervals.bed > start_end_MACS2/unique_peaks_erato_missing_in_melp_.0.5.txt

# merge erato and melp peaks
cat erato_peaks_start_end_pan.rowsort.clean.txt melp_peaks_start_end_pan.rowsort.clean.txt > erato_melp_peaks_start_end_pan.rowsort.clean.bed
sort -k1,1 -k2,2n erato_melp_peaks_start_end_pan.rowsort.clean.bed > erato_melp_peaks_start_end_pan.rowsort.clean.sort.bed
bedtools merge -i erato_melp_peaks_start_end_pan.rowsort.clean.sort.bed > erato_melp_peaks_start_end_pan.rowsort.clean.sort.merged.bed
