
# Creating BED files from GTF annotation file for Nano3p-seq analyses 


## Create a BED file that contains Exon cordinates

```R
Rscript --vanilla exons.R <GTF_file> <Label>

#Example
Rscript --vanilla exons.R zebrafish_chr1.gtf zebrafish
```

## Create a BED file that contains Gene cordinates

```R
Rscript --vanilla genes.R <GTF_file> <Label>

#Example
Rscript --vanilla genes.R zebrafish_chr1.gtf zebrafish
```

## Create a BED file that contains all transcript end coordinates

```R
Rscript --vanilla transcript_ends.R <GTF_file> <Label>

#Example
Rscript --vanilla transcript_ends.R zebrafish_chr1.gtf zebrafish
```



## Create a BED file that contains miRNA gene cordinates

```R
Rscript --vanilla mirna_genes.R <GTF_file> <Label>

#Example
Rscript --vanilla mirna_genes.R zebrafish_chr1.gtf zebrafish
```

## Create a BED file that contains smallRNA gene cordinates

```R
Rscript --vanilla smallrna_genes.R <GTF_file> <Label>

#Example
Rscript --vanilla smallrna_genes.R zebrafish_chr1.gtf zebrafish
```

## Create a BED file that contains smallRNA transcript end coordinates

```R
Rscript --vanilla smallrna_transcriptends.R <GTF_file> <Label>

#Example
Rscript --vanilla smallrna_transcriptends.R zebrafish_chr1.gtf zebrafish
```

