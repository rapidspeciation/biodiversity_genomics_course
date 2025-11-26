# Phylogenomics using SNP data: species trees and phylogenetic networks
Please note that this tutorial has been taken and modified from [Kevin Sanchez github](https://k-sanchez.github.io/workshop_networks_xcbh/ )

## Species tree inference from SNP data using SVDquartets

`svdquartets` is an algorithm that computes species trees directly from SNP data. However, it is not a full-likelihood approach since the data is summarized as pooled site-pattern counts. This algorithm is implemented in <span style="font-variant: small-caps;">PAUP*</span> (it can be downloaded from [http://phylosolutions.com/paup-test/](http://phylosolutions.com/paup-test/)).

We will use a SNPs matrix collected through RADseq from species belonging to the *Liolaemus kingii* group. This group comprises lizards distributed in the Patagonian Steppe and is characterized by a complex diversification history as a result of rapid diversifications and gene flow between species ([Sánchez et al. 2023](https://doi.org/10.1093/sysbio/syad019)).

Now, we want to go to the liolaemus folder

```shell
cd /home/genomics/scratch/users/melisa_o/phylogenomics/liolaemus
ls
```

The file liolaemus_snps.nex is a nexus matriz with one SNP retained per RAD locus, and a list of commands at the end to execute SVDquartets. Explore it using:
```shell
head liolaemus_snps.nex
tail liolaemus_snps.nex
```

At the bottom of the nexus file you should read

outgroup lineomaculatus_0 lineomaculatus_1;

svdq evalQuartets=random nquartets=100000 taxpartition=species bootstrap=standard nreps=100 nthreads=2;

rootTrees rootMethod=outgroup;

savetrees file=SVDquartets.tre format=Newick brlens=yes

- `outgroup`: define outgroup samples in the matrix
- `svdq`: calls the SVDquartets algorithm
- `evalQuartets`: use "x" random quartets (number specified in the next flag)
- `nquartets`: number of quartets to sample
- `taxpartition`: this is the partition that specifies the individual-species associations (already included at the bottom of the `.nex` matrix, you can check this in a text editor, e.g. Notepad)
- `bootstrap`: perform standard bootstrap
- `nreps`: number of pseudoreplicates for bootstrap support
- `nthreads`: number of threads to run in parallel
- `rootTrees`: root tree using the outgroup
- `savetrees`: save trees under `SVDquartets.tre` name

As the nexus file already have all the commands, we can simply run paup

```shell
~/programs/paup4a169_ubuntu64 liolaemus_snps.nex
```

Once it is done, skip paup typing
```shell
quit
```
To check the new file use
```shell
ls
```

## Phylogenetic network reconstruction based on SNP data with PhyloNetworks

Phylogenetic networks are an extension of phylogenetic trees used to model gene flow events between species or populations. Specifically, these  events are modelled by reticulation edges that summarizes gene flow that might have occurred over a period of time into a single instantaneous event. These edges have a parameter associated ($\varphi$) that represents the proportion of alleles transferred during the entire period of gene flow.

A popular method to estimate phylogenetic networks is <span style="font-variant: small-caps;">PhyloNetworks</span> ([Solís-Lemus et al. 2017](https://doi.org/10.1093/molbev/msx235)). It is based on a pseudolikelihood function over concordance factors (CFs) of quartets of taxa (four taxa), wich increase computational tractability. CFs are calculated as the proportion of gene trees supporting the three possible splits in a given quartet:

### Calculation of Concordance Factors from SNP data

It is also possible to get a CF table bypassing gene tree reconstructions, by simply having SNP data as input. Note that each SNP is assumed to be unlinked. Thus, if you have a dataset such as RADseq you want to randomly sample one SNP per locus. For whole genome sequences, sample one SNP separated far enough to reduce linkage, for example 1 SNP every 10 Kb or 50 Kb. 

Different R functions were developed to compute the table of CFs directly from a SNPs matrix ([Olave and Meyer 2020](https://doi.org/10.1093/sysbio/syaa005)), which we will be using here on the *Liolaemus* dataset.


In Toko server we need to load R as module and then start R
```shell
module load ~/modules/R
R
```
Once in R, load the functions with the SNPs2CF that we already uploaded to Toko. This functions are needed to perform calculations of concordance factors from SNPs and generate a CF table as input to PhyloNetworks

```R
source("/home/genomics/programs/functions_v1.7.R")
setwd("/home/genomics/scratch/users/melisa_o/phylogenomics/liolaemus")

```

If we want to convert the vcf to a phylip format, we can use the function vcf2phylip that is on the functions just loaded

```R
vcf2phylip(wd=getwd(), vcf.name="liolaemus_snps.vcf", total.SNPs=8645, random.phase = T, replace.missing = T, output.name=NULL, cores=1)
```

This will create a phylip matrix with the SNPs contained in the vcf

Then run SNPs2CF using this phylip matrix
```R
SNPs2CF(seqMatrix = "liolaemus_snps.phy", ImapName = "Imap.txt", outgroupSp = "lineomaculatus",
         indels.as.fifth.state = FALSE,  
         bootstrap = TRUE, boots.rep = 100, 
         outputName = "SNPs2CF.csv",
         n.quartets = "all", between.sp.only = TRUE,
         save.progress = FALSE,
         cores = 2)
```

# once its done, quit R using
```R
quit()
```

Check what's new on the folder
```sh
ls
```

Now, everything is ready to infer a phylogenetic network in PhyloNetworks using this CF table calculated based on SNP data

```sh
julia
```

Now we are inside Julia. To load the package type:

```julia
using Distributed
addprocs(1)
@everywhere using PhyloNetworks
using CSV, DataFrames, RCall, SNaQ
```

Then set working directory in julia and read the CF table
```julia
CF = readtableCF("SNPs2CF.csv")
```

Now read SVDquartets species tree 

```julia
sppTree = readTopology("SVDquartets.tre");
```

and run a species tree estimation (no yet with gene flow) with SNaQ!
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
using PhyloPlots
net_h1 = readTopology("(Paste,Network)Here))")
plot(net_h1, showgamma = true, style = :majortree, arrowlen = 0.2)

```

The networks returned by the method are not rooted, so it is convinient to include an outgroup species in the datset to properly root the networks after the estimation. In our dataset, the outgroup species is *Liolaemus lineomaculatus*:

```julia
net_h1.names # explore the names of the terminals
rootatnode!(net_h1, "lineomaculatus") # root at the outgroup node
plot(net_h1, showgamma = true, style = :majortree, arrowlen = 0.2) # plot again
```
