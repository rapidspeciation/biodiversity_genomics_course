# Phylogenomics

For datasets where there are multiple species of different divergence times, the data can be used to infer a phylogenetic tree of these species. Ideally, the monomorphic sites (=invariant sites, sites where all individuals are the same) should be included as iqtree will use a substitution model to infer the best phylogeny. You can specify the substitution model or you can let iqtree infer the best one. Here, we will use just the beginning of chromosome 9 (file called Mechanitis_subset.vcf.gz) but in reality, I would recommend to use the entire genome and subsample the genome randomly, e.g. with the vcfrandomsample tool from [vcflib](https://github.com/vcflib/vcflib). Subsampling is recommended as the file with all sites might be too large to run iqtree on.

```shell
cd
mkdir iqtree
cd iqtree
```
In order to infer a phylogenetic tree from our vcf file, we first need to convert the vcf file into a phylip file (or another alignment format like nexus). For that we can use a tool I wrote:

```shell
python3 vcf2phylip.py -i ~/Share/Mechanitis/Mechanitis.vcf.gz -o Mechanitis.phylip
```

Now let's run iqtree in its simplest form:
```shell
iqtree2 -s Mechanitis.phylip -m GTR -T 2 --prefix Mechanitis
```
we could also have done bootstrapping, by adding `-bb 1000` (1000 rapid bootstraps) but in the interest of time, we will not do that as this takes quite long. I would always recommend bootstrapping though so that you know how trustworthy different parts of the tree are.

Now let's download the tree file in a separate terminal not connected to the Amazon server.

```shell
scp -i c1.pem user1@<IP>:~/iqtree/Mechanitis.treefile ./
```

Now you can plot the phylogeny in R or if you prefer, e.g. in [Figtree](http://tree.bio.ed.ac.uk/software/figtree/).

Here an example R script to plot the phylogeny.

```r
#install.packages("BiocManager") # if you do not have ggtree or phytools, run this line by removing the #
#BiocManager::install("ggtree") # if you do not have ggtree, run this line by removing the #
#BiocManager::install("phytools") # if you do not have phytools, run this line by removing the #

# load the ggtree library
library(ggtree)

#get the tree
tree <- read.tree("Mechanitis.treefile")

# Plot the phylogeny with ggtree (note that as iqtree gives an unrooted phylogeny, we root it first with phytools::reroot)
ggtree(phytools::reroot(tree, node.number = 30, position = 0.005),layout='rectangular') +
       geom_tiplab()
```
