## Detecting selective sweeps with extended haplotype statistics
(Note, this tutorial is an adapted version of the speciation genomics tutorial by Mark Ravinet and Joana Meier: https://speciationgenomics.github.io/haplotypes)

As another measure of selection, we will use tests that rely on extended haplotype lengths. During a selective sweep, a variant rises to high frequency so rapidly that linkage disequilibrium with neighbouring polymorphisms is not disrupted by recombination, giving rise to long haplotypes. In regions of low recombination, all haplotypes are expected to be longer than in regions of high recombination. Therefore, it is important to compare the haplotypes against other haplotypes at the same genomic region. This can either be within a population, whereby an ongoing sweep would lead to a single haplotype being very long compared to the other haplotypes, or between populations whereby the population that experienced a sweep has longer haplotypes than the population which was not affected by selection in that genomic region. Haplotypes are best assessed with phased dataset.

To compute detect regions with extended haplotype lengths, we will use the excellent `rehh` R package. For more information see the [vignette](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html). `rehh` is well maintained, continually updated and has [a very informative tutorial](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html) which we recommend you also check out.

To run `rehh` and perform our analyses, we need to run things in `R`. You can either download the data from [our github](https://github.com/speciationgenomics) and run it locally on your own machine, or we can use our `RStudio` server. Once `R` is running, we are ready to go!

### Setting up the R environment
The first thing we need to do is clear our `R` environment, load the packages we need and move into the right directory to have access to our files. Like so:

```r
# clear environment
rm(list = ls())
# load packages
library(rehh)
library(tidyverse)

# Move into the working directory:
setwd("~/workshop_materials/28_ehh/data")
```

### Reading in data from vcfs

Conveniently, `rehh` is able to read in (and also filter) your data from a vcf. However, it is quite tricky to split up individuals so instead we will read in a vcf for each population. We read in data using the `data2haplohh` function:

```r
# read in data for each species
# house
house_hh <- data2haplohh(hap_file = "house_chr8.vcf.gz",
                   polarize_vcf = FALSE)
# bactrianus
bac_hh <- data2haplohh(hap_file = "bac_chr8.vcf.gz",
                         polarize_vcf = FALSE)
```

This will take a few moments but once the commands are run, we will have read in our data. It is really important that we set `polarize_vcf` to `FALSE` because we have not used an outgroup genome to set our alleles as derived or ancestral. **NB. Some plots and functions in `rehh` will still refer to ancestral and derived even if you use this option**. Instead `rehh` will use the minor and major alleles (in terms of frequency) from our data.

Next, we will filter our data on a **minor allele frequency** or **MAF**. This is really simple in `rehh` with the `subset` function:

```r
# filter on MAF - here 0.05
house_hh_f <- subset(house_hh, min_maf = 0.05)
bac_hh_f <- subset(bac_hh, min_maf = 0.05)
```

This will remove some sites - and we'll be ready to run our haplotype scans.

### Performing a haplotype genome scan - *iHS*

Before we can calculate the statistics we are interested in - *iHS* and *xpEHH* - we need to calculate *iES* statistics. Luckily this is really easy using the `rehh` function `scan_hh`.

```r
# perform scans
house_scan <- scan_hh(house_hh_f, polarized = FALSE)
bac_scan <- scan_hh(bac_hh_f, polarized = FALSE)
```

Note that we once again set the `polarized` argument to `FALSE`. Next we can use the output of this scan to calculate *iHS*. We do this with the `ihh2ihs` function.

```r
# perform iHS on house
house_ihs <- ihh2ihs(house_scan, freqbin = 1)
```
ihh2ihs computes the log ratio of iHH of two focal alleles as described in Voight et al. (2006).

 The standardization is performed within each bins separately because of the frequency-dependence of expected iHS values under neutrality. An implicit assumption of this approach is that each bin is dominated by neutral markers. However, as we do not have polarised data, (i.e. we do not know which allele is ancestral or derived), we cannot use that information and set `freqbin = 1`. If we did, `rehh` can apply weights to different bins of allele frequencies in order to test whether there is a significant deviation in the *iHS* statistic.

house_ihs has two elements:
ihs: a data frame with markers in rows and the columns for chromosome name, marker position, iHS and, if standardized, p-value in a negative log10 scale. Optionally, allele frequencies are included.
frequency.class: a data frame with bins in rows and columns for the number of markers, mean uniHS, standard deviation uniHS, lower quantile uniHS, upper quantile uniHS.

Having calculated the statistic, let's plot the results. We can either plot the statistic itself like so:

```r
ggplot(house_ihs$ihs, aes(POSITION, IHS)) + geom_point()
```

Or we can plot the log *P*-value to test for outliers.

```r
# plot
ggplot(house_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()
```

Here a log *P*-value of 6 is equivalent to something like P = 10<sup>-6</sup> - which is quite a conservative threshold for an outlier!

### Performing a haplotype genome scan - *xpEHH*

Next we will calculate *xpEHH* which is the cross-population *EHH* test. This is essentially a test for the probability that if we randomly sampled haplotypes from different populations, we would get different haplotypes. Again, `rehh` makes this simple with the `ies2xpehh` function.

```r
# perform xp-ehh
house_bac <- ies2xpehh(bac_scan, house_scan,
                       popname1 = "bactrianus", popname2 = "house",
                       include_freq = T)
```

Here we provide the names of our previous *iES* scans (`bac_scan` and `house_scan`). We can also provide the function with the names of our populations and finally, if we set `include_freq` to `TRUE`, we get the frequencies of alleles in our output, which might be useful if we want to see how selection is acting on a particular position.

Next, we can plot the *xpEHH* values, like so:

```r
# plot
ggplot(house_bac, aes(POSITION, XPEHH_bactrianus_house)) + geom_point()
```

In this plot, highly negative values suggest selection in population 2 (house in this case) whereas positive values indicate selection in population 1. Alternatively, like with *iHS*, we could plot the log *P* values.

```r
ggplot(house_bac, aes(POSITION, LOGPVALUE)) + geom_point()
```

### Examining haplotype structure around a target of selection

One other nice feature of `rehh` is that we can examine haplotype structure around SNPs we think might be under selection. Before we do that, we need to identify the SNP in our dataset with the strongest evidence of being an *xpEHH* outlier, but it could also be a variant you found with another approach such as GWAS or another selection or introgression test.

```r
# find the highest hit
hit <- house_bac %>% arrange(desc(LOGPVALUE)) %>% top_n(1)
# get SNP position
x <- hit$POSITION
```

Here we also set the position of our putative selection SNP as the object `x`. This is because we need to identify where it occurs in our haplotype objects - unfortunately we cannot use the position for this. In the code below, we find the marker id for both our datasets.

```r
marker_id_h <- which(house_hh_f@positions == x)
marker_id_b <- which(bac_hh_f@positions == x)
```

Now we are ready to plot the bifurcation of haplotypes around our site of selection. We do this like so:

```r
house_furcation <- calc_furcation(house_hh_f, mrk = marker_id_h)
bac_furcation <- calc_furcation(bac_hh_f, mrk = marker_id_b)
```

We can also plot both of these to have a look at them:

```r
plot(house_furcation)
plot(bac_furcation)
```

Note that the two colours represent the two alleles. It says that the blue allele is ancestral but this is actually the reference allele. If the reference is an outgroup, it may indeed often be the ancestral allele but for our purpose, we just treat them as two different alleles.

Calculating the furcation pattern also makes it possible to calculate the haplotype length around our signature of selection.

```r
house_haplen <- calc_haplen(house_furcation)
bac_haplen <- calc_haplen(bac_furcation)
```

With the haplotype length calculated, we can now plot this to see how haplotype structure differs between our two populations.

```r
par(mfrow=c(1,2))
plot(house_haplen)
plot(bac_haplen)
```

Here we can see the blue haplotype is much larger around this target and is also more numerous in the European house sparrow.

### Writing out the data for later used

Finally, before we move on to the last tutorial, we are going to write out the data. We'll also make the column names all smaller letters, to make downstream scripting a bit easier.

Use the following code to achieve this:

```r
# write out house bactrianus xpEHH
house_bac <- tbl_df(house_bac)
colnames(house_bac) <- tolower(colnames(house_bac))
write_tsv(house_bac, "house_bac_xpEHH.tsv")
```

In the next tutorial, we will use `R` to identify genes that are close to our outlier SNPs.
