## Filtering vcf

In the last session, we learned how to call variants and handle VCFs. In this session, we are going to focus on how to filter VCFs. This might seem like a relatively straightforward task but it is actually exceptionally important and something you should spend a lot of time thinking carefully about.

### How many unfiltered variants?

One thing we didn't check yet is how many variants we actually have. Each line in the main output of a vcf represents a single call so we can use the following code to work it out:

```shell
bcftools view -H sara_sapho.vcf.gz| wc -l
```
We have close to 48,393 variants in our full VCF. At present, we have applied no filters at all. This is intentional - we want to see what happens when filters are applied. However, it is also a good idea to perform an initial analysis, to get an idea of how to set filters. However as we have just seen, it takes time to perform operations on a large VCF.

For this reason, it is a good idea to subsample our variant calls and get an idea of the general distribution of a few key attributes of the data.

### Generating statistics from a VCF

In order to generate statistics from our VCF and also actually later apply filters, we are going to use `vcftools`, [a very useful and fast program for handling vcf files](https://vcftools.github.io/examples.html).

Determining how to set filters on a dataset is a bit of a nightmare - it is something newcomers (and actually experienced people too) really struggle with. There also isn't always a clear answer - if you can justify your actions then there are often multiple solutions to how you set your filters. What is important is that you get to know your data and from that determine some sensible thresholds.

Luckily, `vcftools` makes it possible to easily calculate these statistics. In this section, we will analyse our VCF in order to get a sensible idea of how to set such filtering thresholds. The main areas we will consider are:

* **Depth:** You should always include a minimum depth filter and ideally also a maximum depth one too. Minimum depth cutoffs will remove false positive calls and will ensure higher quality calls too. A maximum cut off is important because regions with very, very high read depths are likely repetitive ones mapping to multiple parts of the genome.
* **Quality** Genotype quality is also an important filter - essentially you should not trust any genotype with a Phred score below 20 which suggests a less than 99% accuracy.
* **Missing data** How much missing data are you willing to tolerate? It will depend on the study but typically any site with >25% missing data should be dropped.

#### Setting up

Before we calculate our stats, lets make a little effort to make our commands simpler and also to ensure the output is written to the right place. First we need to make a directory for our results.

```shell
mkdir ~/vcftools
```
Next we will declare to variables to save us some typing below.

```shell
VCF=~/vcf_real/sara_sapho.vcf.gz
OUT=~/vcftools/sara_sapho
```
#### Calculate mean depth per individual

Next we calculate the mean depth of coverage per individual.

```shell
vcftools --gzvcf $VCF --depth --out $OUT
```
#### Calculate mean depth per site

Similarly, we also estimate the mean depth of coverage for each site.

```shell
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
```

#### Calculate proportion of missing data per individual

Another individual level statistic - we calculate the proportion of missing data per sample.

```shell
vcftools --gzvcf $VCF --missing-indv --out $OUT
```
#### Calculate proportion of missing data per site

And more missing data, just this time per site rather than per individual.

```shell
vcftools --gzvcf $VCF --missing-site --out $OUT
```
---

With the statistics calculated, take a moment to have a quick look at the output in the `~/vcftools/` directory. We will now need to download our output data onto our local machines in order to work with R.

```shell
scp -i c28.pem user28@IP:~/vcftools ./
```

Examining statistics in R
------------------------

### Setting up the R environment

First load the `tidyverse` and `ggplot2` packages and ensure you have moved the `vcftools` output into the working directory you are operating in. You may want to set up an RStudio Project to manage this analysis. See [here](https://speciationgenomics.github.io/more_advanced_R/) for a guide on how to do this.

``` r
# load tidyverse package
library(tidyverse)
library(ggplot2)
```

Variant based statistics
------------------------

The first thing we will do is look at the statistics we generated for each of the variants in our subset VCF - quality, depth, missingness and allele frequency.

### Variant mean depth

Next we will examine the mean depth for each of our variants. This is essentially the number of reads that have mapped to this position. The output we generated with `vcftools` is the mean of the read depth across all individuals - it is for both alleles at a position and is not partitioned between the reference and the alternative. First we read in the data.

``` r
var_depth <- read_delim("sara_sapho.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)
```

Again take a moment to look at the data - `mean_depth` is our column of interest but note that you can also get a an idea of the variance in depth among individuals from the `var_depth` column. Once again, we will use `ggplot` to look at the distribution of read depths.

``` r
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() + xlim(0, 100)
```


Hmm - this plot is a bit misleading because clearly, there are very few variants with extremely high coverage indeed. Let's take a closer at the mean depth:

``` r
summary(var_depth$mean_depth)
```

    ##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
    ##    3.50   16.13      17.39      17.46   18.29     1269.41 

Since we all took different subsets, these values will likely differ slightly but clearly in this case most variants have a depth of 17x whereas there are some extreme outliers.

This gives a better idea of the distribution. We could set our minimum coverage at the 5 and 95% quantiles but we should keep in mind that the more reads that cover a site, the higher confidence our basecall is. 10x is a good rule of thumb as a minimum cutoff for read depth, although if we wanted to be conservative, we could go with 15x.

What is more important here is that we set a good **maximum depth** cufoff. As the outliers show, some regions clearly have extremely high coverage and this likely reflects mapping/assembly errors and also paralogous or repetitive regions. We want to exclude these as they will bias our analyses. Usually a good rule of thumb is something the mean depth x 2 - so in this case we could set our maximum depth at 34x.

**So we will set our minimum depth to 10x and our maximum depth to 34x.**

### Variant missingness

Next up we will look at the proportion of missingness at each variant. This is a measure of how many individuals **lack** a genotype at a call site. Again, we read in the data with `read_delim`.

``` r
var_miss <- read_delim("sara_sapho.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)
```

Then we plot the data. One thing to keep in mind here is that different datasets will likely have different missingness profiles. RAD-sequencing data for example is likely to have a slightly higher mean missingnes than whole genome resequencing data because it is a random sample of RAD sites from each individual genome - meaning it is very unlikely all individuals will share exactly the same loci (although you would hope the majority share a subset).

``` r
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
```

Our sara/sapho data has clearly most individuals have a call at almost every site. Indeed if we look at the summary of the data we can see this even more clearly.

``` r
summary(var_miss$fmiss)
```

    ##    Min. 1st Qu.  Median     Mean    3rd Qu.    Max.
    ## 0.000 0.000   0.00877   0.016250 0.026316 0.096491

Most sites have almost no missing data. Although clearly, there are sum (as the max value shows). This means we can be quite conservative when we set our missing data threshold. We will remove all sites where **over 10% of individuals are missing a genotype**. One thing to note here is that `vcftools` inverts the direction of missigness, so our 10% threshold means **we will tolerate 90% missingness** (yes this is confusing and counterintuitive... but that's the way it is!). Typically missingness of 75-95% is used.

Individual based statistics
---------------------------

As well as a our per variant statistics we generated earlier, we also calculated some individual metrics too. We can look at the distribution of these to get an idea whether some of our individuals have not sequenced or mapped as well as others. This is good practice to do with a new dataset. A lot of these statistics can be compared to other measures generated from the data (i.e. principal components as a measure of population structure) to see if they drive any apparent patterns in the data.

### Mean depth per individual

First we will look at the distribution of mean depth among individuals. We read the data in with `read_delim`:

``` r
ind_depth <- read_delim("sara_sapho.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)
```

Then we plot the distribution as a histogram using `ggplot` and `geom_hist`.

``` r
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
```

The mean depth is around 17. However, there is one individual with a very low mean depth (approximately 8).

### Proportion of missing data per individual

Next we will look at the proportion of missing data per individual. We read in the data below:

``` r
ind_miss  <- read_delim("sara_sapho.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
```

This is very similar to the missing data per site. Here we will focus on the `fmiss` column - i.e. the proportion of missing data.

``` r
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()
```

This shows that the proportion of missing data per individual is very small. For most individuals, it ranges between 0.01 and 0.16, so we can confidently say that our individuals were well sequenced. However, there is one individual with a higher proportion of missing data (0.28).

Considering these two results, you should decide whether to remove individuals that have either low mean depth and/or a high amount of missing data.

### Applying filters to the VCF

Now we have an idea of how to set out thresholds, we will do just that. First of all, we will set some simple variables in order to make our filtering command more straightforward.

```shell
VCF_IN=~/vcf_real/sara_sapho.vcf.gz
VCF_OUT=~/vcf_real/sara_sapho_filtered.vcf.gz
```
Then next we will set our chosen filters like so:

```shell
# set filters
MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=30
```
Finally we run the following `vcftools` command on the data to produce a filtered vcf. We will investigate the options as the filtering is running.

```shell
# move to the vcf directory
cd vcf_real
# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indv D5252__Hvenez \
--remove-indels --max-missing $MISS \
--minQ $QUAL \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT
```

What have we done here?

* `--gvcf` - input path -- denotes a gzipped vcf file
* `--remove-indv` - name of the individuals to be removed (Based on the results of mean depth and missing data)
* `--remove-indels` - remove all indels (SNPs only)
* `--max-missing` - set minimum missing data. A little counterintuitive - 0 is totally missing, 1 is none missing. Here 0.9 means we will tolerate 10% missing data.
* `--minQ` - this is just the minimum quality score required for a site to pass our filtering threshold. Here we set it to 30.
* `--min-meanDP` - the minimum mean depth for a site.
* `--max-meanDP` - the maximum mean depth for a site.
* `--recode` - recode the output - necessary to output a vcf
* `--stdout` - pipe the vcf out to the stdout (easier for file handling)

Now, how many variants remain? There are two ways to examine this - look at the vcftools log or the same way we tried before.

```shell
cat out.log
bcftools view -H sara_sapho_filtered.vcf.gz | wc -l
```

You can see we have substantially filtered our dataset!

---
