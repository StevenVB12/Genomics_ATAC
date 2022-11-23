# Genomics_ATAC

### 1. PAN genome assembly

Uses the seq-seq-pan software (https://gitlab.com/rki_bioinformatics/seq-seq-pan) to generate an .xmfa multiple genome alignment file. Genomes and gff files are available through lepbase.org.

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

### 3. ATAC-seq read mapping, peak calling and counting. 

Note that this bash code is written for a slurm array job and would still require to add slurm job submission parameters.

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

Combine MACS peaks to creat peak reference set
```
# Retain peaks present in 2 samples with at least 1 bp overlap
Run_ATAC_combine_peaks_2samples_1bp_per_tissue_stage.sh

# Retain peaks present in all samples with at least 1 bp overlap
Run_ATAC_combine_peaks_ALLsamples_1bp_per_tissue_stage.sh

# Retain peaks present in 2 samples with at least 50% reciprocal overlap
Run_ATAC_combine_peaks_2samples_50perc_per_tissue_stage.sh

# Retain peaks present in all samples with at least 50% reciprocal overlap
Run_ATAC_combine_peaks_ALLsamples_50perc_per_tissue_stage.sh
```

### 4. ATAC-seq QC

Count number of reads in .bam file.
```
samtools view -c SAMPLE.bam
```

Calculate fraction of reads in peaks (FRIP).
```
Run_ATAC_Fraction_reads_in_peaks.sh
```

Calculate Transcription Start Site enrichment score.
```
TSS_enrichment.R
```

### 5. Counting reads in peaks

Counting.
```
# Count read number in reference peak set in H. erato
Run_ATAC_count_H_erato.sh

# Count read number in reference peak set in H. melpomene
Run_ATAC_count_H_melpomene.sh
```

Map  ATAC-seq peak intervals to PAN genome.
```
map_peak_intervals_to_pan.sh
```

### 6. Map bedgraphs and other features to PAN genome

Python code to transform scaffold positions to genome positions for mapping with seq-seq-pan map.
```
seq-seq-pan_bedgraph_chrompos.py
```

Map MACS2 peaks sets to PAN genome coordinates and intersect.
```
Run_ATAC_map_peaks_to_PAN.sh
```

Create bedgraph files from .bam files.
```
Run_bedgraphs.sh
```

Scale and average bedgraphs for tissue/time.
```
Run_bedgraphs_scale.sh
```
 
Map bedgraphs to PAN genome assembly.
```
# Extract positions from bedgraph files for mapping with seq-seq-pan map
Run_bedgraphs_map_to_PAN_preprocess.sh

# Map position to pan genome
Run_bedgraphs_map_to_PAN_seqseqpan_mapping.sh

# Combine start and end positions of intervals mapped to pan genome
Run_bedgraphs_map_to_PAN_postprocessing.sh
```

### 7. Differential accessibility analyses

Calculate size factors (used when scaling and combining bedgraphs).
```
ATAC_sizefactor.R
```

R code for differential accessibility (DA) analyses.
Includes:
- Code to match H. erato and H. melpomene peak counts
- code for DA between developmental time points, wings and sections
- Code for PCA
- Merging of DA peaks with DNA sequence conservation
- Foldchange correlation analysis
- Code to output DA peak sets
```
ATAC_DA_erato_melp_development.R
ATAC_DA_erato_melp_FWHW.R
ATAC_DA_erato_melp_sections.R
```

### 8. Differential expression analyses

Download data
```
Download_RNAseq_data.sh
```

Map RNA-seq reads and count gene expression.
```
Run_mapping.sh
Run_count.sh
```

Create counts table from individual sample mappings.
```
Create_counts_tables.R
```

Differential expression forewing versus hindwing.
```
diff_expression_analysis_FWHW.R
```

Volcano plots with highlighted genes near DA peaks.
```
diff_expression_volcano_FWHW.R
```

Identify genes with shared expression patterns forewing versus hindwing.
```
shared_genes_with_peaks.R
```

Differential expression forewing sections.
```
diff_expression_analysis_sections.R
```


Identify genes close to DA ATAC-seq peak.
```
Run_homer_development.sh
Run_homer_FWHW.sh
Run_homer_sections.sh
```

Correlate ATAC-seq accessibility with gene expression.
```
# Over development
shared_unique_development_expression.R

# between forewing and hindwing
shared_unique_FWHW_expression.R
```

### 9. ATAC-seq peak conservation

Python script to calculate conservation of intervals between pair of genomes.
```
seq-seq-pan_bedfile_conservation.py
```

Calculate conservation for different interval sets.
```
Run_interval_conservation.sh
```

Summaries and visualisation in R.
```
calculate_IDY_mean.R
```

### 10. TF enrichment

Run meme-chip to find TF enrichment patterns between tissues and time points.
```
Run_meme_development.sh
Run_meme_FWHW.sh
Run_meme_sections.sh
```

Parse meme-chip html outputs and extract TFs and enrichment values.
```
parse_MEME.py
```

### 11. Color pattern analysis of optix CRE mutants

R script to extract color and compare with wild type phentypes.
```
patternize_optix4.R
```

### 12. Visualize a PAN genome segment with ATAC-seq data

<i>Ubx</i>
```
Plot_PAN_ATAC_Ubx.R
```

<i>optix</i>
```
Plot_PAN_ATAC_optix.R
```