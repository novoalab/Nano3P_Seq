
# Nano3P-seq
Analysis of Nano3P-seq nanopore libraries (direct cDNA first strand sequencing with template switching)

## Table of Contents
- [General command line steps used to analyze Nano3P-seq datasets](#General-command-line-steps-used-to-analyze-Nano3P-seq-datasets)
- [1. Base-calling, demultiplexing and mapping](#1-Base-calling-demultiplexing-and-mapping)
 * [2. Filtering mapped reads based on annotations and assigning reads to gene biotype](#2-Filtering-mapped-reads-based-on-annotations-and-assigning-reads-to-gene-biotype)
 * [3. Assigning reads to transcripts/isoforms](#3-Assigning-reads-to-transcripts/isoforms)
 * [4. Per-read poly(A) tail length estimations](#4-Per-read-polyA-tail-length-estimations)
 * [5. Post-processing of polyA tail length estimations](#5-Post-processing-of-polyA-tail-length-estimations)
 * [6. Visualizing the results](#6-Visualizing-the-results)
- [Software versions used](#Software-versions-used) 
- [Citation](#Citation) 

## General command line steps used to analyze Nano3P-seq datasets

### 1. Base-calling, demultiplexing and mapping
Base-calling and demultiplexing with Guppy:
```bash
guppy_3_6_1 --device cuda:0 -c dna_r9.4.1_450bps_hac.cfg  --fast5_out -ri input_fast5 -s output_fastq
```
Mapping with minimap2:
```bash
# Mapping to transcriptome
minimap2 -ax map-ont --MD reference.fasta input.fastq | samtools view -hSb -F 3844 - > output.sam
samtools view  -f 0x10 -bq 59 output.sam | samtools sort - output.sorted && samtools index output.sorted.bam

# Mapping to genome
minimap2 -ax splice -k14 -uf --MD $ref input.fastq | samtools view -hSb -F 3844 - >  output.sam
samtools sort output.sam output.sorted && rm output.sam && samtools index output.sorted.bam

```

### 2. Filtering mapped reads based on annotations and assigning reads to gene biotype

#### 2.1. Convert the BAM into BED
```bash
# File needed:
## BAM mapped to genome and rRNAs and sorted

# Convert BAM to BED
bedtools bamtobed -i data.bam > data.bed
```

#### 2.2. Extract Reads starts 

```bash
# File needed : 
## BED converted from BAM mapped to genome and rRNAs and sorted
Rscript --vanilla executable_R_scripts/extract_readstarts.R data.bed <outputFile.readstarts.bed>
```



#### 2.3. Intersect the read starts with Small RNAs transcript ends

```bash
# File needed :
## Reads Starts file of the BAM mapped to Genome and rRNAs and sorted
## BAM Mapped to Genome and rRNAs and sorted
## A file containing the Transcript Ends of SmallRNAs 


bedtools intersect -a data.readstarts.bed -b SmallRNA_TranscriptEnds.bed -wa -wb > data_smallRNAs.bed


#Remove duplicate reads and export read IDs
awk '!seen[$4]++' data_smallRNAs.bed | cut -f4 > data_smallRNAs.reads


### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=data.bam \
       O=data.smallRNAs.bam\
       READ_LIST_FILE=data_smallRNAs.reads \
       FILTER=includeReadList

#Sort the bam
samtools sort data.smallRNAs.bam data.smallRNAs.sorted
#Index the bam
samtools index data.smallRNAs.sorted.bam



### EXTRACT THE REST TO INTERSECT AGAIN
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=data.bam \
       O=data.restRNAs.bam\
       READ_LIST_FILE=data_smallRNAs.reads \
       FILTER=excludeReadList

#Sort the bam
samtools sort data.restRNAs.bam data.restRNAs.sorted
#Index the bam
samtools index data.restRNAs.sorted.bam
```


#### 2.4. Convert the BAM into BED 

```bash
bedtools bamtobed -i data.restRNAs.sorted.bam > data.restRNAs.sorted.bed
```

#### 2.5. Extract Reads starts of the Rest Reads (non-smallRNA reads)

```bash
# File needed : 
## BED converted from BAM mapped to genome and rRNAs and sorted

Rscript --vanilla executable_R_scripts/extract_readstarts.R data.restRNAs.sorted.bed <outputFile.readstarts.bed>

```



#### 2.6. Intersect Read Starts with Transcript Ends for the Rest RNAs (non-smallRNA reads)
```bash
# File needed :
## Reads Starts file of the BAM mapped to Genome and rRNAs and sorted (Rest RNAs)
## BAM Mapped to Genome and rRNAs and sorted
## A file containing the Transcript Ends of All RNAs 


bedtools intersect -a data.restRNAs.readstarts.bed -b All_TranscriptEnds.bed -wa -wb > data.restRNAs.complete.bed



#Remove duplicate reads and extract Read IDs
awk '!seen[$4]++' data.restRNAs.complete.bed | cut -f4 > data.restRNAs.complete.reads


### EXTRACT THE BAM FOR Rest Complete Reads
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=data.bam \
       O=data.restRNAs.complete.bam\
       READ_LIST_FILE=data.restRNAs.complete.reads \
       FILTER=includeReadList

#Sort the BAM file
samtools sort data.restRNAs.complete.bam data.restRNAs.complete.sorted
#Index the BAM file
samtools index data.restRNAs.complete.sorted.bam

#Convert BAM to BED
bedtools bamtobed -i data.restRNAs.complete.sorted.bam > data.restRNAs.complete.sorted.bed

#Remove duplicate reads and extract Read IDs
awk '!seen[$4]++' data.restRNAs.complete.sorted.bed | cut -f4 > data.restRNAs.complete.sorted.reads

```


#### 2.7. Create the BAM file for the miRNAs
```bash
#Files needed
## miRNA Gene BED file extracted from annotation file
## BAM mapped to genome and rRNAs and sorted
## Read ID file from the restRNAs (non-smallRNAs)
## Read ID file from the small RNAs


# Intersect BED with miRNA Genes BED
bedtools intersect -abam data.bam -b miRNA_Gene.bed -wa -wb -bed -split -S > data.miRNAs.bed

# Remove duplicate reads and extract read IDs
awk '!seen[$4]++' data.miRNAs.bed | cut -f4 > data.miRNAs.reads

# Remove the reads from miRNAs that are overlapping with restRNAs and smallRNAs
diff data.restRNAs.complete.sorted.reads data.miRNAs.reads |grep ">"|cut -c 3- > data.miRNAs.RestExcluded.reads


diff data_smallRNAs.reads data.miRNAs.RestExcluded.reads |grep ">"|cut -c 3- > data.miRNAs.RestandSmallExcluded.reads


# Now we processed the miRNA reads file, exluding overlaps with genes
#cDNA964321_all.genome38_sequin_overlapping_miRNAs_RestSmallExcluded.reads

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=data.bam \
       O=data.miRNAs.bam\
       READ_LIST_FILE=data.miRNAs.RestandSmallExcluded.reads\
       FILTER=includeReadList

samtools sort data.miRNAs.bam data.miRNAs.sorted
samtools index data.miRNAs.sorted.bam
```

#### 2.8. Merge all the BED files with their annotations

```bash

#Intersect BAM file of Rest RNAs (non-smallRNAs) with the Rest Exon Annotation BED File
bedtools intersect -abam data.restRNAs.complete.sorted.bam -b Rest_EXON.bed -wa -wb -bed -split -S | awk '!seen[$4]++'>  data.restRNAs_FINAL.bed

#Intersect BAM file of Small RNAs with the Small RNA Gene Annotation BED File
bedtools intersect -abam data.smallRNAs.sorted.bam -b SmallRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > data.smallRNAs_FINAL.bed

#Intersect BAM file of miRNAs with the miRNA Gene Annotation BED File
bedtools intersect -abam data.miRNAs.sorted.bam -b miRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > data.miRNAs_FINAL.bed

cat  data.restRNAs_FINAL.bed data.smallRNAs_FINAL.bed data.miRNAs_FINAL.bed > data.allRNAs.bed
```



### 3. Assigning reads to transcripts/isoforms
```bash
#Prerequisite : isoquant.py tool
python isoquant.py --genedb gtf_file --complete_genedb --bam data.bam --data_type nanopore -o OUTPUT_FOLDER
```

### 4. Per-read poly(A) tail length estimations
```R
#Prerequisite : tailfindR tool
tails <- find_tails(fast5_dir ='fast5_location',
save_dir = './',
csv_filename = 'Tails.csv' ,
num_cores = 10)

```
### 5. Post-processing of polyA tail length estimations 
```bash
Rscript --vanilla executable_R_scripts/process_tail.R Tails.csv <outputFile.tails.processed.tsv>
```

### 6. Visualizing the results

#### a) Scatterplots of poly(A) tail length estimations across biological replicates
```bash
Rscript --vanilla executable_R_scripts/scatter_tails_replicates.R Rep1.tails Rep2.tails Rep1.bed Rep2.bed label
```
#### b) Dotplots of poly(A) tail length estimations across time points 
```bash
Rscript --vanilla executable_R_scripts/dotplot_timepoints.R tails 2hpf.bed 4hpf.bed 6hpf.bed
```
#### c) Line plots of transcript abundances across time points
```
Rscript --vanilla executable_R_scripts/line_plot.R Tails 2hpf.bed 4hpf.bed 6hpf.bed gene_list.txt
```

#### d) Getting your poly(A) tail ends to be visible in IGV

```bash
#Porechop tool to remove the adapter sequences
porechop -i input.fastq> trimmed.fastq

# Remap the reads to the same reference
# Load the BAM file into IGV
# Go to View>Preferences>Alignments>Click on "Show soft-clipped bases"
# Zoom into the 3'end of the gene
```

## Software versions used

* Guppy version 3.6.1
* minimap2 version 2.17
* samtools version 0.1.19
* R version 3.6.0
* TailfindR v1.2
* picard.jar v2.25.0
* bedtools v2.29.1
* Isoquant v1.3
* porechop v0.2.4

## Citation
If you find this work useful, please cite: 


Begik O, Liu H, Delgado-Tejedor A, Kontur C, Giraldez AJ, Beaudoin JD, Mattick JS and Novoa EM. Nano3P-seq: transcriptome-wide analysis of gene expression and tail dynamics using end-capture nanopore sequencing. bioRxiv 2021. doi: https://doi.org/10.1101/2021.09.22.461331. 
