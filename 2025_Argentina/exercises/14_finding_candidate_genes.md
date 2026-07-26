## Finding candidate genes
(Note, this tutorial is an adapted version of the speciation genomics tutorial by Mark Ravinet and Joana Meier: https://speciationgenomics.github.io/candidate_genes)

Now that we have SNPs that we think might be under selection, we can take the next step to identify what the genes lying close those SNPs are. We will do this entirely in `R`.

### Setting up the R environment

The first thing we need to do is clear our `R` environment and load the packages we need. Like so:

```r
# clear environment
rm(list = ls())
# load packages
library(tidyverse)
```
That's it! Just the tidyverse for this section.

```r
# read in the selection scan data
house_bac <- read_tsv("house_bac_xpEHH.tsv")
```

Now we're ready to proceed

### Reading in a gff file

In order to identify genes, we need a `gff` file - which is a gene annotation (or feature) file. You can learn more about the `gff` format [here](https://useast.ensembl.org/info/website/upload/gff.html). For this tutorial, we will use a subset of the house sparrow genome annotation produced by [Elgvin et al (2017)](https://advances.sciencemag.org/content/3/6/e1602996.full). The version we use here is just for chromosome 8.

The gff file is already in the workshop folder, if you are working on your private computer, you can get it [here](https://github.com/speciationgenomics/data/raw/refs/heads/master/house_sparrow_chr8.gff).

```r
# read in the gff
gff <- read_tsv("~/workshop_materials/28_ehh/data/house_sparrow_chr8.gff", col_names = FALSE)
# subset and clear up the gff - add names
colnames(gff) <- c("chr", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute")
```

Before we use the `gff`, we will subset it further so it only includes genes and also rearrange it to sort the order. Finally, we will use `mutate` to add a new column which contains the midpoint of the genes.

```r
# select genes only
new_gff <- gff %>% filter(feature == "gene")
# arrange the gff by position
new_gff <- new_gff %>% arrange(start, end)
# make a gene mid point variable
new_gff <- new_gff %>% mutate(mid = start + (end-start)/2)
```

Next we will identify our peak of selection and look for the genes close to it.

### Identifying genes close to a region of selection

First of all, let's replot our selection scan and remind ourselves where our peak is.

```r
# plot selection scan again
ggplot(house_bac, aes(position, logpvalue)) + geom_point()
# Let's zoom in to the largest peak
ggplot(house_bac, aes(position, logpvalue)) + geom_point() + xlim(1.5E7,2.5E7)
```

How can we find out the identify of the highest peak here? We can do this with a few simple `dplyr` commands.

```r
# identify the highest peak of selection
hits <- house_bac %>% arrange(desc(logpvalue)) %>% top_n(3)
```

For the rest of this tutorial, we will focus on the highest point - the first position in the `hits` `data.frame`.

```r
# find the nearest genes to our highest hit
x <- hits$position[1]
```

Next we will alter the `gff` to include a new column that shows the **absolute distance** from our top selection hit.

```r
# find hits closest to genes
new_gff %>% mutate(hit_dist = abs(mid - x)) %>% arrange(hit_dist)
```

Now we can use this column to identify the genes that occur within **250 Kb** of our target.

```r
# find hits within 250 Kb
gene_hits <- new_gff %>% mutate(hit_dist = abs(mid - x)) %>% arrange(hit_dist) %>% filter(hit_dist < 250000)
```

Last of all let's find out what these genes are...

```r
# what are these genes?
gene_hits <- gene_hits %>% select(chr, start,end, attribute,hit_dist)
# separate out the attribute column
gene_hits %>% pull(attribute)
```

Let's label the nearest genes on a zoomed in peak region:
```r
# as the attributes are very long, let's extract just the key name
new_gff <- new_gff %>% mutate(gene=gsub(".*Note=Similar to (.+):.*","\\1",new_gff$attribute))

# Plot the peak region zoomed in (with base R for a change)
par(mfrow=c(1,1))
plot(house_bac$position,house_bac$logpvalue,cex=0.5,pch=19,
     xlab="position",ylab="log p value",xlim=c(1.8e7,2.1e7), ylim=c(0,12))
rect(xleft = new_gff$start,xright = new_gff$end,ybottom = 9.5,ytop=10)
text(new_gff$mid,10.2,labels = new_gff$gene,srt=45,cex=0.8,adj = 0)
```


So here we have [COL11A1](https://www.genecards.org/cgi-bin/carddisp.pl?gene=COL11A1) and [AMY2A](https://www.genecards.org/cgi-bin/carddisp.pl?gene=AMY2A).
