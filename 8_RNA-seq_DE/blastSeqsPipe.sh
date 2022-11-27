#!/usr/bin/bash

in_file="$1"
fmt=\"6 qseqid sseqid slen sstart send evalue pident sseq\"
while read -r database sample_id; do
    blastn -db "../BLAST_DBs/"$database -query "./DNA-seqs/"$sample_id".exons.bed.fa" -outfmt 6\ qseqid\ sseqid\ slen\ sstart\ send\ evalue\ pident\ sseq -evalue 1e-10 -max_target_seqs 1 -num_threads 12 -out "intermediate_files/"$sample_id".ATAC-v-"$database".blast_out"
    python REtogenome_blast.py "intermediate_files/"$sample_id".ATAC-v-"$database".blast_out" "intermediate_files/"$sample_id".ATAC-v-"$database".seqs.fa"
    blastn -db "chrLibs/"$sample_id -query "intermediate_files/"$sample_id".ATAC-v-"$database".seqs.fa" -outfmt 6\ qseqid\ sseqid\ slen\ sstart\ send\ evalue\ pident\ sseq -evalue 1e-10 -max_target_seqs 1 -num_threads 12 -out "intermediate_files/"$sample_id".ATAC-v-"$database".recip.blast_out"
    python2 genometoRE_blast.py "intermediate_files/"$sample_id".ATAC-v-"$database".recip.blast_out" "output/"$sample_id".ATAC-v-"$database".txt"
done < "$in_file"
