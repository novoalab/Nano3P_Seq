<<<<<<< HEAD
# Nano3P-seq study

This github repo contains the scripts used in the study below : 

Beğik O, Liu H, Delgado-Tejedor A, Kontur C, Giraldez AJ, Beaudoin JD, Mattick JS, Novoa EM (2021) Nano3P-seq: transcriptome-wide analysis of gene expression and tail dynamics using end-capture nanopore sequencing. 


## Each analyses is seperated into different folders


- [Analyses of curlcakes](https://github.com/novoalab/Nano3P_Seq/tree/master/curlcakes)


- [Analyses of yeast total RNA](https://github.com/novoalab/Nano3P_Seq/tree/master/yeast)


- [Analyses of mouse nuclear/mitochondrial RNA](https://github.com/novoalab/Nano3P_Seq/tree/master/mouse)


- [Analyses of zebrafish RNA](https://github.com/novoalab/Nano3P_Seq/tree/master/zebrafish)

=======
# Nano3P-seq
Analysis of Nano3P-seq nanopore libraries (direct cDNA first strand sequencing with template switching)

## General command line steps used to analyze Nano3P-seq

### 1. Base-calling, demultiplexing, mapping
Base-calling and demultiplexing with Guppy:
```

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

* Guppy version XX
* minimap2 version XX
* R version XX
* Nanopolish version XX
* TailfindR version XX

## Citation
If you find this work useful, please cite: 

Begik O, Liu H, Delgado-Tejedor A, Kontur C, Giraldez AJ, Beaudoin JD, Mattick JS and Novoa EM. Nano3P-seq: transcriptome-wide analysis of gene expression and tail dynamics using end-capture nanopore sequencing. bioRxiv 2021. doi: XXX. 



>>>>>>> 6a36f7669be491a93cc5802bf97616426f3e69e8
