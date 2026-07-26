---
title: IQ-tree
---

# Phylogenomics

## Intro
Whole-genome information can be used to reconstruct evolutionary history and divergence times of genes and species, variable genealogy and branch length variation along the genome, presence of horisontal transfer, infer gene family evolution and gene function, ancestral sequence reconstruction, inference of selection, to name some. Whole genome alignment for tree reconstruction in computationally expensive, so depending on the final target for the analysis specific regions of markers in the genome can be used.

Here we will use BUSCO to find single copy orthologs in our genomes. We infer orthogroups with OrthoFinder to reduce the risk of including paralogous genes. Paralogs have a different divergence time compared to the orthologs, which per definition should have the same divergence time as the speciation event.
OrthoFinder nicely output each single copy orthogroup as a multi fasta file that we can directly use for multiple sequence sequence alignment.

To save us some time and computer power we already have set of multi-fasta single copy orthologues from BUSCO from six Ithomiini butterflies from three genera (*Melinaea*, *Mechanitis* and *Napeogenes*) and we are using the monarch (*Danaus plexippus*) as outgroup. Information on how to run BUSCO is in exercise 10_synteny.

## 1 - Orthologs

```
#set up a new working directory from our home directory
cd
mkdir gene_trees
cd gene_trees

cp ~/Share/gene_trees/input_orthofinder/*.fa ./

ls
#take a look at one of the files
less ilMecMaza1.1.primary.fa

#run orthofinder, this time in we are using nucleotides so we must include option -d
orthofinder -f ./ -d
```
As earlier OrthoFinder will produce a large directory with many interesting files. Are there any species that have less genes assigned to orthogroups?
Why do you think this is? How many shared single copy genes do we have?
We will not use all in the following section, we will use a subset of 10 genes.

## 2 - Alignment
First we need to align the gene sequences to infer homology and character state, and align each base so that the differences between the sequences represent evolutionary changes. There are several benchmarked muliple sequence aligner like[ Clustal Omega](http://www.clustal.org/omega/), [Muscle](https://www.drive5.com/muscle/), [T-coffee](https://tcoffee.crg.eu) and [Mafft](https://mafft.cbrc.jp/alignment/software/). They use different methods and are suitable for different types of datasets.If codon information is wanted in downstream applications, codon aware aligners like [PRANK](https://ariloytynoja.github.io/prank-msa/) or [MACSEv2](https://www.agap-ge2pop.org/macse/) should be used, but they are slower and more computer intenstive.

### Mafft
Mafft is a fast multiple sequence aligner which is also available as webinterface.
```
mafft --help

```
There are several options for specifying which algoritm to use, with accuracy-oriented methods suitable for up to 200 sequences and speed-oriented methods for very large alignments up to 50,000 sequences. It also have the very nice option --auto where it determines the most suitable method for the dataset. Let's align our data.

```
#set up directory
cd ~/gene_trees
mkdir mafft
cd mafft
mkdir input output
```
```shell
#copy ten of the genes from orthofinder single copy sequences to our input (here all that starts with OG000004), before copying check how many:

ls ../OrthoFinder/Results_Jul24/Single_Copy_Orthologue_Sequences/OG000004*

# Note, the exact folder name "Results_Jul24" will be different if you run it on another day
#if ok then copy
cp ../OrthoFinder/Results_Jul24/Single_Copy_Orthologue_Sequences/OG000004* input/

#run mafft on one fasta file to see if it works
mafft --auto input/OG0000041.fa > output/OG0000041.msa.fa
```
Take a look
```
less output/OG0000041.msa.fa
```
We should now have an aligned sequence. There are now inserted dashes (-) in the sequence file as representation of gaps between the alignments. The alignments can be visualised in softwares like Bioedit, MEGA or AliView, and edited if needed.

Run mafft for all ten files in a for-loop
```
#align the rest in a for loop
#first create a list of the input files
ls input/ > list_sco.txt
less list_sco.txt

for GENE in $(cat list_sco.txt);do mafft --auto input/${GENE} > output/${GENE%.*}.msa.fa;done

#if mafft do not run you can get the files from the Share folder
cp ~/Share/gene_trees/mafft/output/* output/
```

Trimming the alignment with [trimAL](https://vicfero.github.io/trimal/), which is a very handy tool for handling, filtering, and quality control of sequence alignments. We can remove gaps or remove species if we want to subset the dataset. It can also be used for format conversion.

```
cd ../
mkdir trimal
mkdir output log
ls ../mafft/output/ > list_trimal_input.txt
for FILE in $(cat list_trimal_input.txt)
 do trimal -in ../mafft/output/${FILE} -fasta -out output/${FILE%.*}_trimmed.fa -sgt -nogaps -htmlout log/${FILE%.*}_trimmed.html > log/${FILE%.*}_trimmed.summary.txt
done

#check removal
grep -A2 "Residues" log/*txt
```


## 3 - A gene tree

### IQ-tree
You have already familiarized yourself with [IQ-tree](http://www.iqtree.org/doc/iqtree-doc.pdf). Now we will use it to reconstruct our gene trees.

```
#start from the gene tree folder
cd ~/gene_tree/
mkdir iqtree
cd iqtree
mkdir output

```

Now we can run iqtree. To keep track of the trees we can use a --prefix setting, that specifies the path and the base name of the output files.

```
iqtree -s ../trimal/output/OG0000041.msa_trimmed.fa --prefix output/OG0000041
```
Take a look in the `output` dir

```
ls output/
#check out the main iqtree result file
less output/OG0000041.treefile
less output/OG0000041.iqtree
```

There is a lot of useful information about the alignment and the sequence composition, it also gives the best substitution model (default: lowest BIC) and the rate parameters and substitution matrix under this model, for more details on the models, look [here](http://www.iqtree.org/doc/Substitution-Models). And you can find a representation of the tree.

### Substitution model
A substitution model infering the actual number of substitutions based on the observed differences between the sequences while accounting for multiple hits and different substitution probabilities.The simplest model is Jukes-Cantor (JC) which assumes equal probabilites of for all types of nucleotide substitutions. Other models acount for AT-biased mutation rate due to increased rate of deamination of C -> U, and that the transition and transversion rates are different. The General Time Reverable (GTR) uses different rates for all substitutions. There are also complex models accounting for rate variation across sites and lineages. To determine the most suitable model IQ-tree computes the log-likelihoods for different substitution models together with the the Akaike information criterion (AIC), corrected Akaike information criterion (AICc), and the Bayesian information criterion (BIC). The information criteria are likelihood-based estimates of model fit that penalises number of additional parameters.The lower AIC or BIC the better model fit. You can also use specific models by specifying the option `-m`.


### How robust is our tree hypothesis?

Assess branch support using bootstrapping. Bootstrapping uses resampling with replacement from the dataset to create rearranged datasets and repeats the generation of the phylogenetic tree. The values represent how many times out of a hundred the same branch is observed, it shows the robustness of the recontructed tree from the given data.
We can specify the best model from the last run. We also have to specify a new prefix otherwise iqtree will raise an error to avoid overwriting files.

```shell
iqtree -s  ../trimal/output/OG0000041.msa_trimmed.fa -m TN+F+I -B 1000 --prefix output/OG0000041_bs1000
```
Check the bootstrap support values in the output/OG0000041_bs1000.iqtree.

## 4 - Multiple gene trees

We will now take a look at the trees from different gene alignments or loci, do they have the same topology and branch length?
First we need to change the fasta header in all alignments so that is contains only the taxa name, same in all alignments. For this we can use `sed` and remove everything after underscore on each line in the file.

```shell

sed 's/_.*//' ../trimal/output/OG0000041.msa_trimmed.fa |head

#change all files with a for loop
#make a list
ls ../trimal/output/ > list_input_msa.txt
mkdir input
for file in $(cat list_input_msa.txt);do sed 's/_.*//' ../trimal/output/$file > input/${file%.*}.renamed.fa;done
#check a file
less OG0000041.msa.fa

```
Now we have the separate alignments for each locus in a folder, so we can perform the following commands:
```shell
#we can specify a variable for out input directory
INPUT=input/

# infer a concatenation-based species tree with 1000 ultrafast bootstrap
iqtree -p $INPUT --prefix output/concat -B 1000

# infer the locus/gene trees
iqtree -S $INPUT --prefix output/loci

# compute gene concordance factors
iqtree -t output/concat.treefile --gcf output/loci.treefile --prefix output/concord

```

## 5 - tree visualisation

We will use ggtree and [ggdensitree](https://rdrr.io/bioc/ggtree/man/ggdensitree.html) for visualisation of the trees.
Copy the script densitree.R to you `iqtree` folder or to your local computer.

Copy the loci.treefile and concord.cf.tree to your local computer.
```
scp -i c1.pem user1@35.92.168.2:~/Share/gene_trees/densitree.R ./
scp -i c1.pem user1@35.92.168.2:~/gene_trees/iqtree/loci.treefile ./
scp -i c1.pem user1@35.92.168.2:~/gene_trees/iqtree/concord.cf.tree ./
```

```r
#load the libraries we need
library(ggtree)
library(phytools)
library(dplyr)
library(ggdensitree)

#read in the species tree, in newick format
tree <- read.tree("concord.cf.tree")

#create the plot
concat_tree <-
ggtree(tree, ladderize = T) +
	geom_tiplab() +
	geom_treescale()

#take a look
concat_tree

#we can add some improvement
concat_tree <-
ggtree(tree, ladderize = T) +
  geom_tiplab(align = T) +
  geom_treescale() +
  xlim(c(0,0.4)) +		#adjust the limits of the x axis so that we can see the tiplabels
  geom_nodelab(hjust = 1.1, vjust = -0.6) + #plots the node information (bootstrap and concordance factor)
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue") #add the branch length

concat_tree

#we can save the tree if we want
ggsave(plot=concat_tree, filename = "concat_tree.png",
       height = 8,
       width = 6)

```

IQ-tree like many ML-tree builders renders an unrooted tree. We have an outgroup so we can root the tree on the outgroup.
```r
#check tiplabel or nodelabel index
tree$tip.label

#reroot with phytools with outgroup
ggtree(phytools::reroot(tree,  node.number = 1), ladderize = T) +
  geom_tiplab(align = T) +
  geom_nodelab() +
scale_x_continuous(limits = c(0,0.4))+
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue")

#specify the position of the root,
#node this position is arbitrary, we do not have the information unless we an another outgroup
concat_tree_rerooted <-
ggtree(phytools::reroot(tree,  node.number = 1, position = 0.1), ladderize = T) +
  geom_tiplab(align = T) +
  geom_nodelab() +
  scale_x_continuous(limits = c(0,0.4))+
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue")

ggsave(plot=concat_tree_rerooted, filename = "concat_tree_rerooted.png",
       height = 8,
       width = 6)
```

Take a look at the genetrees

```r
#read in the file with multiple trees
trees <- read.tree("loci.treefile")

#create the plot
gene_trees <-
ggdensitree(layout = "slanted", trees, color="lightblue") +
  geom_tiplab() +
  geom_treescale()

gene_trees

```
Now we can see as in the output of iqtree (.iqtree), that the trees have very different branch length. Why is that do you think?

We can make a cladogram to visualise the different topologies of the trees more clearly.

```r
#set branch.length="none", this will render a cladogram
ggdensitree(layout = "slanted", trees, branch.length="none", color="blue", alpha = 0.3) +
  geom_tiplab() +
  geom_treescale()


#save the tree if you want
trees_slanted <-
ggdensitree(layout = "slanted", trees, branch.length="none", color="blue", alpha = 0.3) +
  geom_tiplab() +
  geom_treescale()

ggsave(plot = trees_slanted,
       device = "png",
       filename = "trees_rerooted_slanted.png")

```
We can also plot the concatenated tree on top of the multi-tree plot.

```r
##change the concatenated tree to cladogram by setting branch.length "none" and layout "slanted"
concat_tree <-
  ggtree(tree, ladderize = T, layout = "slanted", branch.length = "none") +
  geom_tiplab(align = T) +
  geom_treescale() +
  #xlim(c(0,0.4)) +
  geom_nodelab(hjust = 1.1, vjust = -0.6) +
  geom_text2(aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "blue")

#extract the data information to add to the densitree plot
d1 <- trees_slanted$data
d2 <- concat_tree$data
#move the start point of the root
d2$x <- d2$x - max(d2$x)
#check that they are the same
d1$x
d2$x

#Plot the trees
trees_slanted +
  geom_tree(data = d2, layout = "slanted") +
  geom_nodelab(data = d2, hjust = 1.1, vjust = -0.6) +
  geom_text2(data = d2, aes(label = round(branch.length, 4)), hjust = 0.6, vjust = 1.5, colour = "orange4")

```
Does the concatenated tree give a good representation of the gene trees?
For many applications it is more informative to reconstruct a tree for each loci to represent the evolutionary history of that particular loci. The discrepancy between gene trees both in branch length and topology gives valuable information on the evolutionary processes acting across the genome, while a concatenated tree is at best an approximate summary.
