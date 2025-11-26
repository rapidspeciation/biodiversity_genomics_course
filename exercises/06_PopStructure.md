## Principal components analysis (PCA)

Here, we will investigate population structure using [principal components analysis](https://en.wikipedia.org/wiki/Principal_component_analysis). Examining population structure can give us a great deal of insight into the history and origin of populations. Model-free methods for examining population structure and ancestry, such as principal components analysis are extremely popular in population genomic research. This is because it is typically simple to apply and relatively easy to interpret.

Essentially, PCA aims to identify the main axes of variation in a dataset with each axis being independent of the next (i.e. there should be no correlation between them). The first component summarizes the major axis variation and the second the next largest and so on, until cumulatively all the available variation is explained. In the context of genetic data, PCA summarizes the major axes of variation in allele frequencies and then produces the coordinates of individuals along these axes.

To perform a PCA on our Heliconius data, we will use plink - specifically [version 2.0](https://www.cog-genomics.org/plink/2.0/) (although be aware [older](https://www.cog-genomics.org/plink/1.9/) versions are available). Note that plink was written with human data in mind. As a result, we need to provide a bit of extra info to get it to work on our non-human dataset.


### Linkage pruning

One of the major assumptions of PCA is that the data we use is indpendent - i.e. there are no spurious correlations among the genomic sites. This is not the case for most genomic data as allele frequencies are correlated due to physical linkage and linkage disequilibrium. We also want to get a genome-wide picture and avoid that e.g. a large inversion strongly affects the results. So as a first step, we need to prune our dataset of variants that are strongly linked (=highly correlated alleles).

There are two options when it comes to keeping only independent sites.

(1) **Filter based on linkage disequilibrium (LD) directly calculated from the data.** It makes sense to use this approach when the dataset is composed by a randomly mating population (i.e. a single species). 

plink --indep-pairwise


(2) **Filter based on physical distance**. When two or more species are present in your dataset, it only makes sense to calculate LD for each species separately. You may want prune the datase by calculating LD for one of the species (provided large enough samples sizes). Alternatively, you may chose to calculate LD-decay (i.e. how LD values decrease with physical distance) and then prune based on physical distance between SNPs (e.g. average minimum distance between SNPs when LD ≤ 0.2). Sometimes, this is already known, as in the case of *Heliconius* in which LD ≤ 0.2 after ~10-kb.
(Note: it might be worth comparing LD-decay among the different species in the dataset, provided large enough sample sizes)

plink --bp-space


First things first, we will make a directory called 06_PopulationStructure where we will run our analysis

```shell
# move to your directory
cd /home/genomics/scratch/<yourname>

# make a new directory for this excercise
mkdir pca

# move into it
cd pca
```

Since we are only interested in investigating population structure within the melpomene-timareta-cydno clade, we will exclude outgroup (*H. numata*) individuals from our dataset. This can be done using the `--keep` option in plink.
```shell
# To avoid having lots of copies of this vcf file, you can directly specify the full path to the file without copying it to your folder
VCF="/home/genomics/scratch/data/martin2019/wgenome.martin2019.biallelic.mac2.vcf.gz"

# create a file listing individuals in ingroup (melpomene, timareta, cydno), excluding the outgroup "Hnum". This file will be used in plink to keep only these individuals

bcftools query -l $VCF | grep -v Hnum > mel_tim_cyd.keep

```

PCA requires only variable SNPs that are bi-allelic (2 alternative alleles). Let's keep only the ingroup individuals and apply a filter to get variable sites where the allele is found at least twice `--min-alleles 2 --max-alleles 2 --mac 2`. We will also remove sites where more than 5 individuals have missing data `--geno 5`. Keep only one site per 10000 bp.

```shell
plink2 \
--vcf $VCF \
--threads 2 \
--allow-extra-chr \
--keep mel_tim_cyd.keep \
--geno 0 \
--min-alleles 2 \
--max-alleles 2 \
--mac 2 \
--bp-space 10000 \
--export vcf bgz \
--out wgenome.martin2019.ingroup.mac2.prune10kb
```

Since in plink, filtering commands are processed in a pre-defined [order](https://www.cog-genomics.org/plink/2.0/order) we we can run all these filters with a single command.

Here an explanation for all parameters we used:
- `--vcf` - specified the location of our VCF file.
- `--threads` - number of compute threads to use.
- `--allow-extra-chr` - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome.
- `--keep` - Filter out all samples not named in a file.
- `--geno` - Filter out variants with missing call rates exceeding the provided value.
- `--min-alleles` - Filter out variants with fewer than given # of alleles
- `-max-alleles` - Filter out variants with more than given # of alleles
- `--mac` - Filter out variants with minor allele count lower than #
- `--bp-space` - Remove variants so each pair is no closer than the given distance.
- `--export vcf bgz` - Make a bgz-compressed vcf file

If you wish to prune the dataset based on LD directly estimated from your dataset, you can find a series of commands to do so at the end of this tutorial.


### Perform a PCA

Next we rerun plink with a few additional arguments to get it to conduct a PCA. First, we need to produce a file with allele frequencies per SNP, necessary for the PCA analyses (note this is not necessary in earlier versions of plink)

First, we generate a file with allele frequencies for all sites
```shell
plink2 \
--vcf wgenome.martin2019.ingroup.mac2.prune10kb.vcf.gz \
--threads 2 \
--allow-extra-chr \
--set-missing-var-ids @:# \
--freq \
--out wgenome.martin2019.ingroup.mac2.prune10kb
```

Now we can create our PCA.
```shell
# create pca
plink2 \
--vcf wgenome.martin2019.ingroup.mac2.prune10kb.vcf.gz \
--threads 8 \
--allow-extra-chr \
--set-missing-var-ids @:# \
--read-freq wgenome.martin2019.ingroup.mac2.prune10kb.afreq \
--pca \
--out wgenome.martin2019.ingroup.mac2.prune10kb
```

This is very similar to our previous command. What did we do here?

- `--freq` - this just lets plink know we want to extract only these positions from our VCF - in other words, the analysis will only be conducted on these.
- `--read-freq` - this is necessary to write out some additional files for another type of population structure analysis - a model based approach with admixture.
- `--set-missing-var-ids`- also necessary to set a variant ID for our SNPs. Human and model organisms often have annotated SNP names and so `plink` will look for these. We do not have them so instead we set ours to default to `chromosome:position` which can be achieved in `plink` by setting the option `@:#` - [see here](https://www.cog-genomics.org/plink/1.9/data#set_missing_var_ids) for more info.
- `--pca` - fairly self explanatory, this tells plink to calculate a principal components analysis.


Once the command is run, we will see a series of new files. We will break these down too:

PCA output:

- wgenome.martin2019.ingroup.mac2.prune10kb.**eigenval** - the eigenvalues from our analysis for each PC axis
- wgenome.martin2019.ingroup.mac2.prune10kb.**eigenvec**- the eigenvectors from our analysis for each individual



### Plotting the PCA output

Let's now download the relevant files to our local computers to plot the PCA in R on your own computer. In a new terminal that is not connected to the server, write (changing <yourname> to the name of your users folder)
Make sure to go into the folder where you want to download the files into.

```shell
# define a variable with your users folder name
NAME=joana_m
scp genomics@toko.uncu.edu.ar:/home/genomics/scratch/users/$NAME/pca/wgenome.martin2019.ingroup.mac2.prune10kb.eigenvec ./
scp genomics@toko.uncu.edu.ar:/home/genomics/scratch/users/$NAME/pca/wgenome.martin2019.ingroup.mac2.prune10kb.eigenval ./
scp genomics@toko.uncu.edu.ar:/home/genomics/scratch/data/martin2019/martin2019_info.txt ./
```

#### Setting up the R environment
First load the `tidyverse` package and ensure you have moved the plink output into the working directory you are operating in. You may want to set up an RStudio Project to manage this analysis. See [here](https://speciationgenomics.github.io/more_advanced_R/) for a guide on how to do this.

```r
# load tidyverse package
library(tidyverse)

# set the folder where you downloaded the files into as your working directory
setwd("C:/Users/jm66/Dropbox/teaching/BiodiversityGenomics/Argentina2025/PCA/")

# read in data
pca <- read_table2("wgenome.martin2019.ingroup.mac2.prune10kb.eigenvec", col_names = TRUE)
eigenval <- scan("wgenome.martin2019.ingroup.mac2.prune10kb.eigenval")
info <- read_table2("martin2019_info.txt")
```

#### Cleaning up the data
Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. We will also give our pca data.frame proper column names.

```r
# sort out the pca data

# set names
names(pca)[1] <- "ind"

# add the species information
pca <- as_tibble(merge(pca, info, by="ind"))
```

#### Plotting the data
Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).
```r
# first convert to percentage variance explained
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

# With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.

# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
ylab("Percentage variance explained") + theme_classic()
```
Cumulatively, they explain 100% of the variance but PC1, PC2 and possible PC3 together explain about 54% of the variance. We could calculate this with the cumsum function, like so:

```r
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
```

Next we move on to actually plotting our PCA. Given the work we did earlier to get our data into shape, this doesn't take much effort at all.

```r
# plot pca
ggplot(pca, aes(PC1, PC2, col = pop)) + geom_point(size = 3) +
coord_equal() + theme_classic() +
xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +
ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
```

Note that this R code block also includes arguments to display the percentage of variance explained on each axis. Here we only plot PC1 and PC2. Given that PC3 also shows a high percentage of variance explained, it could be worth it to also plot PC1 against PC3.



### Supplementary: Perform LD-based prunning (if working with individuals from the same species)
Linkage pruning, here in a window of 50 sites will remove sites that show strong linkage, i.e. correlation in allele frequency. So if multiple sites basically contain the same information, we only want to keep one of them.

```shell
# prune snps based on linkage disequilibrium (LD)
plink2 \
--vcf $VCF \
--threads 8 \
--allow-extra-chr \
--bad-ld \
--set-missing-var-ids @:# \
--keep mel_tim_cyd.keep \
--min-alleles 2 \
--max-alleles 2 \
--mac 2 \
--indep-pairwise 50 10 0.2 \
--out wgenome.martin2019.ingroup.mac2

# extract LD-pruned sites
plink2 \
--vcf $VCF \
--threads 8 \
--allow-extra-chr \
--set-missing-var-ids @:# \
--keep mel_tim_cyd.keep \
--min-alleles 2 \
--max-alleles 2 \
--mac 2 \
--extract wgenome.martin2019.ingroup.mac2.prune.in \
--export vcf id-paste=iid \
--out wgenome.martin2019.ingroup.mac2.ld_prune
```

As well as being versatile, plink is very fast. It will quickly produce a linkage analysis for all our data and write plenty of information to the screen. When complete, it will write out two files wgenome.martin2019.ingroup.mac2.prune.in and wgenome.martin2019.ingroup.mac2.prune.out. The first of these is a list of sites which fell below our linkage threshold - i.e. those we should retain. The other file is the opposite of this. In the next step, we will produce a PCA from these linkage-pruned sites.
