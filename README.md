# Genomics_ATAC

### 1. PAN genome assembly

Uses the seq-seq-pan software (https://gitlab.com/rki_bioinformatics/seq-seq-pan) to generate an .xmfa multiple genome alignment file. 

```
Run_seq-seq-pan.sh
```

### 2. PAN genome processing

Python script to create .bed files with sequence blocks present in pan genome for each genome.
```
python seq-seq-pan_blocks_intervals.py -I SeqSeqPan_erato_melp_noNewline.xmfa -g 1,2
```

Bash commands for intersecting block .bed with bedtools.
```
Intersect_bed.sh
```

Transform .xmfa file to fasta.
```
python seq-seq-pan_toFasta.py -I SeqSeqPan_erato_melp_noNewline.xmfa -g 1,2
```

### 3. ATAC-seq read mapping, peak calling and counting. Note that this bash code is written for a slurm array job and would still require to add slurm job submission parameters.

Clean reads using Trimmomatic.
```
Run_ATAC_trimmomatic.sh
```

Read mapping.
```
Run_ATAC_mapping.sh
```

Peak calling MACS2.
```
Run_ATAC_MACS2.sh
```

