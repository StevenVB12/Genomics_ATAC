# Remove first line
sed -i -n -e '2,$p' 1_blocks_intervals.txt
sed -i -n -e '2,$p' 2_blocks_intervals.txt

# Remove space from end of line
sed -i 's/[[:space:]]*$//' 1_blocks_intervals.txt
sed -i 's/[[:space:]]*$//' 2_blocks_intervals.txt

# Add ‘pan’ as scaffold in the .bed file
lineN=$(< 1_blocks_intervals.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(cut -d$'\t' -f 1,2 1_blocks_intervals.txt) > 1_blocks_intervals.corr.bed

lineN=$(< 2_blocks_intervals.txt wc -l)
pr -mts <(yes "pan" | head -n $lineN) <(cut -d$'\t' -f 1,2 2_blocks_intervals.txt) > 2_blocks_intervals.corr.bed

# Subtract blocks to get unique sequences
bedtools subtract -sorted -a 1_blocks_intervals.corr.bed -b 2_blocks_intervals.corr.bed > blocks_unique_1.txt
bedtools subtract -sorted -a 2_blocks_intervals.corr.bed -b 1_blocks_intervals.corr.bed > blocks_unique_2.txt


