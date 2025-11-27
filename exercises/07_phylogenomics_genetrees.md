# Phylogenomics: gene trees, species trees and phylogenetic networks
Please note that this tutorial has been taken and modified from [Kevin Sanchez github](https://k-sanchez.github.io/workshop_networks_xcbh/ )

## Gene tree reconstruction with RAxML

For gene tree reconstructions, we will use a dataset of Australasian monitor lizards (genus _Varanus_) from Pavón-Vázquez et al. ([2021](https://doi.org/10.1093/sysbio/syaa102)). It consists of 388 nuclear loci obtained through anchored hybrid enrchment, a technique for capturing orthologous regions of the genome. To estimate trees from these loci, we will rely on <span style="font-variant: small-caps;">RAxML</span>, a program for efficient tree inference based on maximum likelihood (ML). We will also estimate node supports based on bootstrap calculations.

Login into toko server and connect to toko05 node. Then go into your folder and copy all the data in phylogenomics folder provided by this course

```shell
NAME=melisa_o
cd ~/scratch/users/$NAME
cp -r ~/scratch/data/phylogenomics .
```
Check what's in

```shell
cd phylogenomics
ls
```
The folder monitors is the one we are going to use at this session. Get inside this folder and check what's in

```shell
cd monitors
ls
```
There are two loci in PHYLIP format provided as examples, along with a .zip file containing all 388 loci. Here, we will work only with the two example loci, but we also provide the code to automatically generate gene trees for all 388 loci. 
For Toko server, we first need to load raxml program as module and then run it to get a tree.

Now let's run iqtree in its simplest form:
```shell
module load ~/modules/raxml
raxmlHPC-PTHREADS-SSE3 -s locus177.phylip -n 177.boot -m GTRGAMMA -f a -N 100 -p 2334 -x 563454
```

- `-s`: name of the sequence file (include the path to the file if it is located in a different folder)
- `-n`: name of the output files (the files generated during the run will have `.177.stand` appended to the end)
- `-m`: substitution model
- `-f`: Specify one of the different algorithms available in <span style="font-variant: small-caps;">RAxML</span>. If nothing is specified (like in our first run), by default it executes the standard hill climbing algorithm to perform the tree search (which is equivalent to `-f d`). The `a` option tells <span style="font-variant: small-caps;">RAxML</span> to conduct a rapid Bootstrap analysis and search for the best-scoring ML tree in a single run
- `-N`: number of bootstrap pseudoreplicates
- `-p`: random number seed to generate a parsimony starting tree (can be any integer)
- `-x`: specify an integer number (random seed) and turn on rapid bootstrapping

Further command options are detailed in the software manual, or can be explored using:
```sh
raxmlHPC-PTHREADS-SSE3 -help
```
The maximum likelihood tree is printed in the `RAxML_bestTree.1.stand` file. We can visualize the tree in <span style="font-variant: small-caps;">FigTree</span> (download from [here](https://github.com/rambaut/figtree/releases/tag/v1.4.4)) and, optionally, export in any image format. To visualize this tree and the support values open the file in <span style="font-variant: small-caps;">FigTree</span>. On the left-hand side of the screen select: <button>Branch Labels</button> &rarr; <button>Display</button> &rarr; `label`.

Let's estimate a tree for a different locus:
```sh
raxmlHPC-PTHREADS-SSE3 -s locus256.phylip -n 256.boot -m GTRGAMMA -f a -N 100 -p 2334 -x 563454
```
Check out the files using ls

### Automatizing gene tree inference using a loop
It is possible to automatically set a run for all 388 gene trees using the code for a loop. Note that all `.phy` in the dataset folder are named `L_1.phy`, `L_2.phy` ... `L_388.phy`. Thus, we can set a loop with an iterator `i` taking values from 1 to 388 to call all the input `.phy` into <span style="font-variant: small-caps;">RAxML</span>:

First unzip the files
```sh
unzip all_388loci.zip
cd all_loci
ls
```

Then run a loop to iterate across all 388 loci as follows

```sh
for i in {1..388}
do
./raxmlHPC -s L_$i.phy -n $i.boot -m GTRGAMMA -f a -N 100 -p 2334 -x 563454
done
```

## Species tree reconstruction based on gene trees with ASTRAL
Species tree estimation is mainly based on the multispecies coalescent model (MSC; [Liu et al. 2021](https://doi.org/10.1007/978-1-4939-9074-0_7)). This model accomodates gene trees within species trees, while allowing for incomplete lineage sorting (ILS).

<span style="font-variant: small-caps;">astral</span> belongs to a family of species tree methods known as two-step because it uses estimated gene trees from sequence alignments. Here, we will use the maximum likelihood trees inferred from the 388 alignments of monitors.

The software can be downloaded from [GitHub](https://github.com/smirarab/ASTRAL/archive/refs/heads/master.zip)

In Toko server it is already installed, and we can simply run:

```sh
astral -i monitors_trees.tre -o monitor_sptree.tre 
```

- `-i`: file containing input gene trees in newick format (a single file where each gene tree is in a different line)
- `-o`: filename for storing the output species tree

ASTRAL should have created a file monitor_sptree.tre with the estimated species tree, check using 

```sh
ls
```

## Phylogenetic network reconstruction based on gene trees with PhyloNetworks

Phylogenetic networks are an extension of phylogenetic trees used to model gene flow events between species or populations. Specifically, these  events are modelled by reticulation edges that summarizes gene flow that might have occurred over a period of time into a single instantaneous event. These edges have a parameter associated ($\varphi$) that represents the proportion of alleles transferred during the entire period of gene flow.

A popular method to estimate phylogenetic networks is <span style="font-variant: small-caps;">PhyloNetworks</span> ([Solís-Lemus et al. 2017](https://doi.org/10.1093/molbev/msx235)). It is based on a pseudolikelihood function over concordance factors (CFs) of quartets of taxa (four taxa), wich increase computational tractability. CFs are calculated as the proportion of gene trees supporting the three possible splits in a given quartet:

We are going to estimate a phylogenetic network taking raxml gene trees as input data and the starting tree inferred from ASTRAL.

```sh
julia
```

Now we are inside Julia. To load the package type:

```julia
using Distributed
addprocs(1)
@everywhere using PhyloNetworks
using CSV, DataFrames, SNaQ
```

Then set working directory in julia and read gene trees
```julia
cd("/home/genomics/scratch/users/melisa_o/phylogenomics/monitors")
CF = readtrees2CF("monitors_trees.tre")
```

Now read ASTRAL species tree 

```julia
sppTree = readTopology("monitor_sptree.tre");
```

and run a species tree (no yet with gene flow) estimation with SNaQ!
```julia
net_h0 = snaq!(sppTree, CF, hmax = 0, filename = "net0", runs = 1)
```

and finally lets use this tree (net_h0) as starting tree for the phylogenetic network estimation under the multispecies network coalescent model
```julia
net_h1 = snaq!(net_h0, CF, hmax = 1, filename = "net1", runs = 1)
```
Note: it is recommended to run several runs (usually at least runs=20). Each run can be estimated on a different core, so set addprocs() accordingly


As the server does not have a display, to plot results, we are going to need to download the network to our computer. Use scp to download the network to your computer, or simply copy the network you want to plot and paste in julia as shown below.

Once you have the network, in your computer open julia. 

```julia
Pkg.add("PhyloPlots") # if you do not have PhyloPlots then (package to visualize the networks), then install it
Pkg.add("PhyloNetworks") # if you do not have PhyloPlots then (package to visualize the networks), then install it
using PhyloPlots, PhyloNetworks 
net_h1 = readTopology("(Paste,Network)Here))")
plot(net_h1, showgamma = true, style = :majortree, arrowlen = 0.2)

```

The networks returned by the method are not rooted, so it is convinient to include an outgroup species in the datset to properly root the networks after the estimation.

```julia
net_h1.names # explore the names of the terminals
rootatnode!(net_h1, "timorensis") # root at the outgroup node
plot(net_h1, showgamma = true, style = :majortree, arrowlen = 0.2) # plot again
```



