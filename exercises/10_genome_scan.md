# Performing a genome scan

In the previous sessions, we investigated population structure with PCA and phylogenomics. Now we want to see where in the genome these butterfly species are particularly divergent. In comparisons of species that still show gene flow (are hybridising sometimes), highly divergent genomic regions are likely to show reduced gene flow e.g. due to divergent selection or reproductive isolation barriers. However, if the species are not showing gene flow anymore, highly divergent regions are more likely to have higher mutation rates or selection.

There are many different ways to detect regions under divergent selection or that confer barriers to gene flow. In this tutorial, we are going to compute four of them in genomic windows:
- pi, a measure of genetic variation within a population or species
- Fst, a measure of genomic differentiation between populations or species
- dxy, a measure of absolute divergence between populations or species
- fd/fdM, a measure of gene flow/introgression from one population/species into another

For a tutorial on long-range haplotype statistics to infer selective sweeps, see [here](https://speciationgenomics.github.io/haplotypes/). You may also want to consider more complex methods such as [SweeD](https://academic.oup.com/mbe/article/30/9/2224/999783#74416771) to infer sweeps, or to detect barrier loci: [gIMble](https://europepmc.org/article/ppr/ppr564457), [diem](https://www.biorxiv.org/content/10.1101/2022.03.24.485605v3) or ancestral recombination graph methods [(review)](https://academic.oup.com/genetics/article/221/1/iyac044/6554197).

Note that dxy and pi require monomorphic sites to be present in the dataset, whereas Fst and fd are only computed on bi-allelic sites. It is thus important to filter out indels and multi-allelic sites and to keep monomorphic sites (no maf filter).

As SNP Fst values are very noisy, it is better to compute Fst estimates for entire regions. Selection is expected to not only affect a single SNP and the power to detect a selective sweep is thus higher for genomic regions. How large the genomic region should be depends on the SNP density, how fast linkage disequilibrium decays, how recent the sweep is and other factors. It is thus advisable to try different window sizes. Here, we will use 20 kb windows. An alternative is to use windows of fixed number of sites instead of fixed size. We will use scripts written by [Simon Martin](https://simonmartinlab.org/) which you can download [here](https://github.com/simonhmartin/genomics_general).
Note, that the scripts by Simon are written in Python2 (not Python3 which may be standard in your working environment). If the scripts do not run, you may have to adjust the first line "#!/usr/bin/env python" to your Python2 path.

First, let's convert the vcf file into Simon Martin's geno file. You can download the script [here](https://github.com/simonhmartin/genomics_general/raw/master/VCF_processing/parseVCF.py).

```shell
# make a new folder in your users folder
cd ~
mkdir genome_scans
cd genome_scans

# We will use a vcf file of just a part of chromosome 18 to speed up this exercise.
# Convert the vcf file to geno.gz which is the format that Simon's scripts require
# Note, we do not filter for bi-allelic sites as we need to include monomorphic sites for pi and Dxy. This file does have a filter on missing data (max 10%).
VCF="/scratchsan/C_computacion/fs20sanger_ac/martin2019/Hmel218003o.subset.vcf.gz"

# Convert the vcf file to geno.gz which is the format that Simon Martin's script requires
module load apps/genomic-general/v0.5
parseVCF.py -i $VCF --skipIndels -o Hmel218003o.geno.gz

#/ create a file assigning individuals to populations
conda activate bcftools
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
    -p Hmel.mal.col -p Hmel.agl.per -p Hmel.ama.per -p Hmel.mel.gui -p Hnum.bsl.bra -p Htim.flo.per -p Htim.the.per -p Hcyd.chi.pan -p Hcyd.zel.col \
    --popsFile popmap.txt
```

Note that -w 20000 specifies a window size of 20 kb that is sliding by 20 kb (-s 20000) and -m 10000 requests these windows to have a minimum number of 10 kb sites covered (50% of sites present). 

The way we have encoded the genotypes (e.g. A/T) in our geno.gz file is called "phased" and we specify that with "-f phased" even though our data is actually not phased. Instead of writing all the individual names into the command, we can give only the species names in the code (e.g. -p Hmel.mal.col -p Hmel.agl.per) and with `--popsFile` specify a file that contains a line for each individual with its name and species in a text file.

Next, we calculate fd and fdM to test for introgression between H. melpomene amaryllis and H. timareta thelxiopea using H. numata as outgroup. fd and fdM are measures of introgression suitable for small windows.

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
    -P1 Hmel.mal.col -P2 Hmel.ama.per -P3 Htim.flo.per -O Hnum.bsl.bra \
    --popsFile popmap.txt \
    --writeFailedWindows
```

To speed up the calculation of these statistics, the script can be run on multiple threads by specifying -T \<thread number>.

For this script we need to specify that at least 50% of the individuals of each population need to have data for a site to be considered (-\-minData 0.5) and we reduce m to 100 as it only considers polymorphic sites.

To plot the results, we will use the files I prepared for the complete chr18 (Hmel218003o) and read it into R. You can download all files found in the Share/genome_scan_results folder. So in a separate terminal, type:



```shell
# Copy the popgen and dstat files to your local computer
scp -r -J <user>@168.176.34.122 <user>@perseus:/scratchsan/C_computacion/<user>/genome_scans/popgen.w20s20.csv.gz .
scp -r -J <user>@168.176.34.122 <user>@perseus:/scratchsan/C_computacion/<user>/genome_scans/dstats.w20s20.csv.gz .
# Copy file containing the location of genes of interest in the genome 
scp -r -J fs20sanger_ac@168.176.34.122 fs20sanger_ac@perseus:/scratchsan/C_computacion/fs20sanger_ac/share/colorPatternGenes.csv .

# Unzip the file
gunzip Hmel218003o.popgen.w20s20.csv.gz 
gunzip Hmel218003o.dstats.w20s20.csv.gz 
```

Then we can start plotting:

```r
# load libraries
require(ggplot2)
require(cowplot)
# install.packages('cowplot')

# First set your working directory to that where the files were downloaded to
setwd('/Users/fs20/Documents/2025.BioDivGenomics/exercises/genomic_scans/')

# Read file with information about genes of interest
genes = read.csv('/Users/fs20/Documents/2025.BioDivGenomics/heliconius_martin2019/08_genomeScans/colorPatternGenes.csv')

# Prepare input files:
# Read in the file with sliding window estimates of FST, pi and dxy
windowStats<-read.csv("Hmel218003o.popgen.w20s20.csv.gz",header=T)

# Read in the fd and fdM estimates of 20 kb windows
fstats <- read.csv("Hmel218003o.dstats.w20s20.csv.gz",header=T,na.strings = "NaN")

# Let's have a look at the FST and fdM datasets
head(windowStats)
head(fstats)

# Let's plot FST, dxy and fd between the two younger species
fst_plot = ggplot(windowStats,aes(mid/1000000, Fst_Hmel.mal.col_Hmel.ama.per))+
  geom_point(data = genesSub, aes(x=(start+end)/2000000, y=1), col='red', shape=25) +
  # geom_vline(data = genesSub, aes(xintercept=(start+end)/2000000 ), col='red') +
  geom_point() +
  scale_x_continuous(limits=c(0,NA), expand = c(0.01,0.01), breaks=seq(0,100,1), name="Chromosome Position (Mb)") +
  scale_y_continuous(name=expression(F[ST])) 
fst_plot

dxy_plot = ggplot(windowStats, aes(mid/1000000, dxy_Hmel.mal.col_Hmel.ama.per*100))+
  geom_vline(data = genesSub, aes(xintercept=(start+end)/2000000 ), col='red') +
  geom_point() +
  scale_x_continuous(limits=c(0,NA), expand = c(0.01,0.01), breaks=seq(0,100,1), name="Chromosome Position (Mb)") +
  scale_y_continuous(name=expression(d[XY]))
dxy_plot

fdm_plot = ggplot(fstats, aes(mid/1000000, fdM)) +
  geom_vline(data = genesSub, aes(xintercept=(start+end)/2000000 ), col='red') +
  # geom_point(data = genesSub, aes(x=(start+end)/2000000, y=1), col='red', shape=25) +
  geom_point() +
  scale_x_continuous(limits=c(0,NA), expand = c(0.01,0.01), breaks=seq(0,100,1), name="Chromosome Position (Mb)")


# Let's compare the stats on chr18 next to each other
plot_grid(fst_plot, dxy_plot, fdm_plot, nrow=3, align='hv')

```

Note, if we had more than one chromosome, we can use a package that allows us to plot the chromosomes next to each other, like the manhattan function of the [qqman R package](https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html). This can also be achieved with e.g. ggplot (if your skills in R allow).
