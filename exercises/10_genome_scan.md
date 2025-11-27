# Performing a genome scan

In the previous sessions, we investigated population structure with PCA and phylogenomics. Now we want to see where in the genome these butterfly species are particularly divergent. In comparisons of species that still show gene flow (are hybridising sometimes), highly divergent genomic regions are likely to show reduced gene flow e.g. due to divergent selection or reproductive isolation barriers. However, if the species are not showing gene flow anymore, highly divergent regions are more likely to have higher mutation rates or selection.

There are many different ways to detect regions under divergent selection or that confer barriers to gene flow. In this tutorial, we are going to compute four of them in genomic windows:
- pi, a measure of genetic variation within a population or species
- Fst, a measure of genomic differentiation between populations or species
- dxy, a measure of absolute divergence between populations or species
- fd, a measure of gene flow/introgression from one population/species into another

For a tutorial on long-range haplotype statistics to infer selective sweeps, see [here](https://speciationgenomics.github.io/haplotypes/). You may also want to consider more complex methods such as [SweeD](https://academic.oup.com/mbe/article/30/9/2224/999783#74416771) to infer sweeps, or to detect barrier loci: [gIMble](https://europepmc.org/article/ppr/ppr564457), [diem](https://www.biorxiv.org/content/10.1101/2022.03.24.485605v3) or ancestral recombination graph methods [(review)](https://academic.oup.com/genetics/article/221/1/iyac044/6554197).

Note that dxy and pi require monomorphic sites to be present in the dataset, whereas Fst and fd are only computed on bi-allelic sites. It is thus important to filter out indels and multi-allelic sites and to keep monomorphic sites (no maf filter).

As SNP Fst values are very noisy, it is better to compute Fst estimates for entire regions. Selection is expected to not only affect a single SNP and the power to detect a selective sweep is thus higher for genomic regions. How large the genomic region should be depends on the SNP density, how fast linkage disequilibrium decays, how recent the sweep is and other factors. It is thus advisable to try different window sizes. Here, we will use 20 kb windows. An alternative is to use windows of fixed number of sites instead of fixed size. We will use scripts written by [Simon Martin](https://simonmartinlab.org/) which you can download [here](https://github.com/simonhmartin/genomics_general).
Note, that the scripts by Simon are written in Python2 (not Python3 which may be standard in your working environment). If the scripts do not run, you may have to adjust the first line "#!/usr/bin/env python" to your Python2 path.

First, let's convert the vcf file into Simon Martin's geno file. You can download the script [here](https://github.com/simonhmartin/genomics_general/raw/master/VCF_processing/parseVCF.py).

```shell
# Connect to the server
ssh genomics@toko.uncu.edu.ar
ssh toko05

# make a new folder in your users folder
cd /home/genomics/scratch/users/<yourname>
mkdir genome_scans
cd genome_scans

# We will use a vcf file of just a part of chromosome 18 to speed up this exercise.
# Convert the vcf file to geno.gz which is the format that Simon's scripts require
# Note, we do not filter for bi-allelic sites as we need to include monomorphic sites for pi and Dxy. This file does have a filter on missing data (max 10%).
VCF="/home/genomics/scratch/data/martin2019/Hmel218003o.subset.vcf.gz"

parseVCF.py \
    -i $VCF \
    --skipIndels \
    -o Hmel218003o.geno.gz

# create a file assigning individuals to populations
bcftools query -l $VCF | awk '{print $1"\t"substr($1,1,12)}' > popmap.txt

```
Note: `bcftools query` allows you to manipulate VCF files, extracting specimens names from the header (when using `-l` flag ). The awk command tells it to print the first column "$1", which here is the individual names (e.g. Hcyd.ali.ecu.001), and then a tab "\t" and then the first to 12th character of the first column "substr($1,1,12)", which here means the population name (e.g. Hcyd.ali.ecu).


First, we will calculate pi for each species and Fst and dxy for each pair of species all in one go.
```shell
popgenWindows.py \
    --windType coordinate \
    -g Hmel218003o.geno.gz \
    -o Hmel218003o.popgen.w20s20.csv.gz \
    -w 20000 \
    -s 20000 \
    -m 10000 \
    -f phased \
    -T 2 \
    -p Hmel.mal.col -p Hmel.agl.per -p Hmel.ama.per -p Hmel.mel.gui \
    -p Hnum.bsl.bra -p Htim.flo.per -p Htim.the.per -p Hcyd.chi.pan -p Hcyd.zel.col \
    --popsFile popmap.txt
```

Note that -w 20000 specifies a window size of 20 kb that is sliding by 20 kb (-s 20000) and -m 10000 requests these windows to have a minimum number of 10 kb sites covered (50% of sites present). The way we have encoded the genotypes (e.g. A/T) in our geno.gz file is called "phased" and we specify that with "-f phased" even though our data is actually not phased. Instead of writing all the individual names into the command, we could give only the species names in the code (e.g. -p lysimnia -p polymnia) and with `--popsFile` specify a file that contains a line for each individual with its name and species in a text file.

Next, we calculate fd to test for introgression between H. melpomene amaryllis and H. timareta thelxiopea using H. numata as outgroup. fd is a measure of introgression suitable for small windows.

```shell
ABBABABAwindows.py \
    -g Hmel218003o.geno.gz \
    -o Hmel218003o.dstats.w20s20.csv.gz \
    -f phased \
    -w 20000 \
    -s 20000 \
    -m 100 \
    --minData 0.5 \
    -T 2 \
    -P1 Hmel.mel.gui -P2 Hmel.ama.per -P3 Htim.the.per -O Hnum.bsl.bra \
    --popsFile popmap.txt \
    --writeFailedWindows
```

To speed up the calculation of these statistics, the script can be run on multiple threads by specifying -T \<thread number>.

For this script we need to specify that at least 50% of the individuals of each population need to have data for a site to be considered (-\-minData 0.5) and we reduce m to 100 as it only considers polymorphic sites.

To plot the results, we need will use the files I prepared for the complete chr18 (Hmel218003o) and read it into R. You can download all files found in the Share/genome_scan_results folder. So in a separate terminal, type:


```shell
#Remember change the user name:
scp genomics@toko.uncu.edu.ar:/home/genomics/scratch/data/martin2019/Hmel218003o.hmelv25.mel_tim_cyd_num.dstats.w20s20.csv.gz ./
scp genomics@toko.uncu.edu.ar:/home/genomics/scratch/data/martin2019/Hmel218003o.hmelv25.mel_tim_cyd_num.popgen.w20s20.csv.gz ./

# Unzip the file
gunzip Hmel218003o.hmelv25.mel_tim_cyd_num.dstats.w20s20.csv.gz
gunzip Hmel218003o.hmelv25.mel_tim_cyd_num.popgen.w20s20.csv.gz
```

Then we can start plotting:

```r
# Prepare input files:
# Read in the file with sliding window estimates of FST, pi and dxy
windowStats<-read.csv("Hmel218003o.hmelv25.mel_tim_cyd_num.popgen.w20s20.csv",header=T)

# Read in the fd estimates of 20 kb windows
fd<-read.csv("Hmel218003o.hmelv25.mel_tim_cyd_num.dstats.w20s20.csv",header=T,na.strings = "NaN")

# Let's have a look at the FST and fd datasets
head(windowStats)
head(fd)

# Let's plot FST, dxy and fd between the two younger species
require(ggplot2)
fst<-ggplot(windowStats,aes(mid,Fst_Hmel.mal.col_Hmel.ama.per))+geom_point()
fst
dxy<-ggplot(windowStats,aes(mid,dxy_Hmel.mal.col_Hmel.ama.per))+geom_point()
dxy
fd<-ggplot(fd,aes(mid,fd))+geom_point()
fd
# Let's compare the stats on chr18 next to each other
require(gridExtra)
grid.arrange(fst, dxy, fd, nrow=3)
```

Note, if we had more than one chromosome, we would use a package that allows us to plot the chromosomes next to each other, like the manhattan function of the [qqman R package](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html).
