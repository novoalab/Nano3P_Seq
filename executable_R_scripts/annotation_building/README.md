
# Creating BED files from GTF annotation file for detailed analyses 


## Create a BED file that contains Exon cordinates

```R
Rscript --vanilla exons.R <GTF_file> <Label>

#Example
Rscript --vanilla exons.R Saccer64.gtf yeast
```

## Create a BED file that contains Gene cordinates

```R
Rscript --vanilla genes.R <GTF_file> <Label>

#Example
Rscript --vanilla genes.R Saccer64.gtf yeast
```

## Create a BED file that contains all transcript end coordinates

```R
Rscript --vanilla transcript_ends.R <GTF_file> <Label>

#Example
Rscript --vanilla transcript_ends.R Saccer64.gtf yeast
```



## Create a BED file that contains miRNA gene cordinates

```R
Rscript --vanilla mirna_genes.R <GTF_file> <Label>

#Example
Rscript --vanilla mirna_genes.R Saccer64.gtf yeast
```

## Create a BED file that contains smallRNA gene cordinates

```R
Rscript --vanilla smallrna_genes.R <GTF_file> <Label>

#Example
Rscript --vanilla smallrna_genes.R Saccer64.gtf yeast
```

## Create a BED file that contains smallRNA transcript end coordinates

```R
Rscript --vanilla smallrna_transcriptends.R <GTF_file> <Label>

#Example
Rscript --vanilla smallrna_transcriptends.R Saccer64.gtf yeast
```

