## Principal components analysis (PCA)

Here, we will investigate population structure using [principal components analysis](https://en.wikipedia.org/wiki/Principal_component_analysis). Examining population structure can give us a great deal of insight into the history and origin of populations. Model-free methods for examining population structure and ancestry, such as principal components analysis are extremely popular in population genomic research. This is because it is typically simple to apply and relatively easy to interpret.

Essentially, PCA aims to identify the main axes of variation in a dataset with each axis being independent of the next (i.e. there should be no correlation between them). The first component summarizes the major axis variation and the second the next largest and so on, until cumulatively all the available variation is explained. In the context of genetic data, PCA summarizes the major axes of variation in allele frequencies and then produces the coordinates of individuals along these axes.

To perform a PCA on our Heliconius data, we will use plink - specifically [version 2.0](https://www.cog-genomics.org/plink/2.0/) (although be aware [older](https://www.cog-genomics.org/plink/1.9/) versions are available). Note that plink was written with human data in mind. As a result, we need to provide a bit of extra info to get it to work on our non-human dataset.


### Linkage pruning

One of the major assumptions of PCA is that the data we use is indpendent - i.e. there are no spurious correlations among the genomic sites. This is not the case for most genomic data as allele frequencies are correlated due to physical linkage and linkage disequilibrium. We also want to get a genome-wide picture and avoid that e.g. a large inversion strongly affects the results. So as a first step, we need to prune our dataset of variants that are strongly linked (=highly correlated alleles).

There are two options when it comes to keeping only independent sites.

(1) **Filter based on linkage disequilibrium (LD) directly calculated from the data.** It makes sense to use this approach when the dataset is composed by a randomly mating population (i.e. a single species). 

```shell
plink --indep-pairwise
```

(2) **Filter based on physical distance**. When two or more species are present in your dataset, it only makes sense to calculate LD for each species separately. You may want prune the datase by calculating LD for one of the species (provided large enough samples sizes). Alternatively, you may chose to calculate LD-decay (i.e. how LD values decrease with physical distance) and then prune based on physical distance between SNPs (e.g. average minimum distance between SNPs when LD ≤ 0.2). Sometimes, this is already known, as in the case of *Heliconius* in which LD ≤ 0.2 after ~10-kb. It might be worth comparing LD-decay in your dataset to specify a threshhold. Here a [tutorial](https://speciationgenomics.github.io/ld_decay/) on how to do that. If you do not have a genome assembly, and work with e.g. de novo assembled RAD loci or UCEs, you can use [STACKS](https://catchenlab.life.illinois.edu/stacks/comp/populations.php) to select one SNP per RAD/UCE locus.

```shell
plink --bp-space
```

First things first, we will make a directory called popStructure where we will run our analysis

```shell
# make a directory for this exercise
mkdir ~/popStructure
# move into it
cd ~/popStructure
# activate conda
module load envs/anaconda3
```

Since we are only interested in investigating population structure within the melpomene-timareta-cydno clade, we will exclude outgroup (*H. numata*) individuals from our dataset. This can be done using the `--keep` option in plink.

```shell
# To avoid having lots of copies of this vcf file, you can directly specify the full path to the file without copying it to your folder
VCF="/scratchsan/C_computacion/fs20sanger_ac/06_PopulationStructure/wgenome.martin2019.biallelic.mac2.vcf.gz"

# Get the information file:
cp /scratchsan/C_computacion/jm66sanger_ac/Heliconius_info.txt ./

# create a file listing individuals in ingroup (melpomene, timareta, cydno). This file will be used in plink to keep only these individuals
conda activate bcftools
bcftools query -l $VCF | \
    egrep -v Hnum | \
    awk '{print $1}' > mel_tim_cyd.keep

# Subset to ingroup individuals
module load apps/plink/2.0.0-a
plink2 \
    --vcf $VCF \
    --threads 8 \
    --allow-extra-chr \
    --keep mel_tim_cyd.keep \
    --export vcf id-paste=iid \
    --out wgenome.martin2019.ingroup
```

Note that since we excluded outgroup individuals, some of the sites might not be variable in the ingroup dataset. We need to apply again filters to get only variable sites `--min-alleles 2 --max-alleles 2 --mac 2`. We will also remove sites where not all individuals have information `--geno 0`. (Note: you may want to be more or less stringent with the missing data filter --geno, depending on your dataset).
```shell
# Include bi-allelic sites only (excluding singletons)
plink2 \
    --vcf wgenome.martin2019.ingroup.vcf \
    --threads 8 \
    --allow-extra-chr \
    --geno 0 \
    --min-alleles 2 \
    --max-alleles 2 \
    --mac 2 \
    --export vcf id-paste=iid \
    --out wgenome.martin2019.ingroup.mac2
```

Finally, we can prune the dataset based on physical linkage.
```shell
plink2 \
    --vcf $VCF \
    --threads 8 \
    --allow-extra-chr \
    --keep mel_tim_cyd.keep \
    --bp-space 10000 \
    --export vcf id-paste=iid \
    --out wgenome.martin2019.ingroup.mac2.prune10kb
```

Since in plink, filtering commands are processed in a pre-defined [order](https://www.cog-genomics.org/plink/2.0/order) we can run all the previous filters in one command.
```shell
plink2 \
    --vcf $VCF \
    --threads 8 \
    --allow-extra-chr \
    --keep mel_tim_cyd.keep \
    --min-alleles 2 \
    --max-alleles 2 \
    --mac 2 \
    --bp-space 10000 \
    --export vcf id-paste=iid \
    --out wgenome.martin2019.ingroup.mac2.prune10kb
```

So for our plink command, we did the following:
- `--vcf` - specified the location of our VCF file.
- `--threads` - number of compute threads to use.
- `--allow-extra-chr` - allow additional chromosomes beyond the human chromosome set. This is necessary as otherwise plink expects chromosomes 1-22 and the human X chromosome.
- `--keep` - Filter out all samples not named in a file.
- `--geno` - Filter out variants with missing call rates exceeding the provided value.
- `--min-alleles` - Filter out variants with fewer than given # of alleles
- `-max-alleles` - Filter out variants with more than given # of alleles
- `--mac` - Filter out variants with minor allele count lower than #
- `--bp-space` - Remove variants so each pair is no closer than the given distance.

If you wish to prune the dataset based on LD directly estimated from your dataset, you can find a series of commands to do so at the end of this tutorial.


### Perform a PCA

Next we rerun plink with a few additional arguments to get it to conduct a PCA. First, we need to produce a file with allele frequencies per SNP, necessary for the PCA analyses (note this is not necessary in earlier versions of plink)

First, we generate a file with allele frequencies for all sites
```shell
# generate an allele frequency file
plink2 \
    --vcf wgenome.martin2019.ingroup.mac2.prune10kb.vcf \
    --threads 8 \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --freq \
    --out wgenome.martin2019.ingroup.mac2.prune10kb
```

Now we can create our PCA:
```shell
# create pca
plink2 \
    --vcf wgenome.martin2019.ingroup.mac2.prune10kb.vcf \
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
- `--pca` - fairly self explanatory, this tells plink to calculate a principal components analysis.

Once the command is run, we will see a series of new files. We will break these down too:

PCA output:

- wgenome.martin2019.ingroup.mac2.prune10kb.**eigenval** - the eigenvalues from our analysis
- wgenome.martin2019.ingroup.mac2.prune10kb.**eigenvec**- the eigenvectors from our analysis



### Plotting the PCA output

Let's now download the relevant files to our local computers to plot the PCA in R on your own computer. In a new terminal, write (changing the user number to your user number and the IP by the correct IP number)

```shell
scp -r -J <username>@168.176.34.122 <username>@perseus:~/popStructure/wgenome.martin2019.ingroup.mac2.prune10kb.eigenvec ./
scp -r -J <username>@168.176.34.122 <username>@perseus:~/popStructure/wgenome.martin2019.ingroup.mac2.prune10kb.eigenval ./
scp -r -J <username>@168.176.34.122 <username>@perseus:~/popStructure/Heliconius_info.txt ./
```


#### Setting up the R environment
First load the tidyverse package and ensure you have moved the plink output into the working directory you are operating in. You may want to set up an RStudio Project to manage this analysis. See here for a guide on how to do this.


```shell
# load tidyverse package
library(tidyverse)
```

Then we will use a combination of readr and the standard scan function to read in the data.

```shell
# read in data
pca <- read_table2("Heliconius.eigenvec", col_names = FALSE)
eigenval <- scan("Heliconius.eigenval")
info <- read_table2("Heliconius.info")
```

#### Cleaning up the data
Unfortunately, we need to do a bit of legwork to get our data into reasonable shape. First we will remove a nuisance column (plink outputs the individual ID twice). We will also give our pca data.frame proper column names.

```shell
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# add the species information
pca <- as_tibble(merge(pca, info, by="ind"))
```
#### Plotting the data
Now that we have done our housekeeping, we have everything in place to actually visualise the data properly. First we will plot the eigenvalues. It is quite straightforward to translate these into percentage variance explained (although note, you could just plot these raw if you wished).


```shell
# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
```

With that done, it is very simple to create a bar plot showing the percentage of variance each principal component explains.


```shell
# make plot
ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
ylab("Percentage variance explained") + theme_light()
```

Cumulatively, they explain 100% of the variance but PC1, PC2 and possible PC3 together explain about 54% of the variance. We could calculate this with the cumsum function, like so:

```shell
# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)
```

Next we move on to actually plotting our PCA. Given the work we did earlier to get our data into shape, this doesn't take much effort at all.

```shell
# plot pca
ggplot(pca, aes(PC1, PC2, col = species)) + geom_point(size = 3) +
    coord_equal() + theme_light() +
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
