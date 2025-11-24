# Phylogenomics

For datasets where there are multiple species of different divergence times, the data can be used to infer a phylogenetic tree of these species. Ideally, the monomorphic sites (=invariant sites, sites where all individuals are the same) should be included as iqtree will use a substitution model to infer the best phylogeny. You can specify the substitution model or you can let iqtree infer the best one. Here, we will use just the beginning of chromosome 9 (file called Mechanitis_subset.vcf.gz) but in reality, I would recommend to use the entire genome and subsample the genome randomly, e.g. with the vcfrandomsample tool from [vcflib](https://github.com/vcflib/vcflib). Subsampling is recommended as the file with all sites might be too large to run iqtree on. Also, ideally one would use different methods to infer phylogenies, including methods that are better at capturing variance in phylogenies across the genome, such as [ASTRAL](https://github.com/smirarab/ASTRAL) that combines the phylogenies from different regions of the genome. Some tools such as [BEAST](https://beast.community/) also allow us to infer divergence times in years by calibrating the phylogeny.

```shell
cd
mkdir iqtree
cd iqtree
```
In order to infer a phylogenetic tree from our vcf file, we first need to convert the vcf file into a phylip file (or another alignment format like nexus). For that we can use a tool I wrote. You can download it [here](https://github.com/joanam/scripts/blob/master/vcf2phylip.py):

```shell
python3 vcf2phylip.py -i ~/Share/Mechanitis/Mechanitis_subset.vcf.gz -o Mechanitis.phylip
```

Now let's run iqtree in its simplest form:
```shell
iqtree2 -s Mechanitis.phylip -m GTR -T 2 --prefix Mechanitis
```
The parameter -m is used to specify the substitution model, here GTR. By not specifying this, iqtree will infer the best model. However, this takes very long and in the interest of time, we will thus not do it.
We could also have done bootstrapping, by adding `-bb 1000` (1000 rapid bootstraps) but in the interest of time, we will not do that as this takes quite long. I would always recommend bootstrapping though so that you know how trustworthy different parts of the tree are.

Now let's download the tree file in a separate terminal not connected to the Amazon server.

```shell
scp -i c1.pem user1@<IP>:~/iqtree/Mechanitis.treefile ./
```

Now you can plot the phylogeny in R or if you prefer, e.g. in [Figtree](http://tree.bio.ed.ac.uk/software/figtree/).

Here an example R script to plot the phylogeny.

```R
#install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggtree)

#variable
#get the tree
tree <- read.tree("Mechanitis.treefile")

# Root and plot the phylogeny
ggtree(phytools::reroot(tree, node.number = 30, position = 0.005),
       layout='rectangular') +
  geom_tiplab() +
  xlim(c(0,0.03))

```
iqtree gives us an unrooted phylogeny, which means that we do not know which group is an outgroup to the others. The code above thus specifies that the root is to be placed on the node number 30, which roots the tree making Mechanitis messenoides the outgroup. Alternatively, we could also list all Mechanitis messenoides individuals as part of the outgroup.  
