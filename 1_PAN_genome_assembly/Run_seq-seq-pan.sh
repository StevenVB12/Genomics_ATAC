seq-seq-pan-wga --config genomefile=genome_list.txt outfilename=SeqSeqPan_erato_melp

## genome_list.txt includes the names of each genome fasta file. Fasta files are available through http://lepbase.org/
# Herato_final.fasta
# Hmel2.fa

## Remove newline from xmfa file for downstream processing
perl -pe 'chomp if /^[ATCGNSBDHVMRWYKatcgnsbdhvmrwyk-]/' SeqSeqPan_erato_melp.xmfa | sed 's/\=/\n\=/g' | sed 's/>/\n>/g' | sed '/^$/d'  > SeqSeqPan_erato_melp_noNewline.xmfa
