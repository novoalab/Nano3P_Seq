
# Nano3P-seq
Analysis of Nano3P-seq nanopore libraries (direct cDNA first strand sequencing with template switching)

## General command line steps used to analyze Nano3P-seq

### 1. Base-calling, demultiplexing, mapping
Base-calling and demultiplexing with Guppy:
```bash
guppy_3_6_1 --device cuda:0 -c dna_r9.4.1_450bps_hac.cfg  --fast5_out -ri input_fast5 -s output_fastq
```
Mapping with minimap2:
```

```

### 2. Filtering mapped reads based on annotations
```

```
### 3. Assigning reads to gene biotype
```

```
### 4. Assigning reads to transcripts/isoforms
```

```
### 5. Per-read poly(A) tail length estimations
```

```
### 6. Post-processing of polyA tail length estimations 
```

```
### 7. Analysis of tail composition
```

```
### 8. Visualizing the results

#### a) Scatterplots of poly(A) tail length estimations across biological replicates
```

```
#### b) Dotplots of poly(A) tail length estimations across time points 
```

```
#### c) Line plots of transcript abundances across time points
```

```
### 9. BONUS: Getting your poly(A) tail ends to be visible in IGV



## Software versions used

* Guppy version 3.6.1
* minimap2 version XX
* R version XX
* Nanopolish version XX
* TailfindR version XX

## Citation
If you find this work useful, please cite: 

Begik O, Liu H, Delgado-Tejedor A, Kontur C, Giraldez AJ, Beaudoin JD, Mattick JS and Novoa EM. Nano3P-seq: transcriptome-wide analysis of gene expression and tail dynamics using end-capture nanopore sequencing. bioRxiv 2021. doi: XXX. 



>>>>>>> 6a36f7669be491a93cc5802bf97616426f3e69e8
