# Selection

## Introduction

To test if a gene, or a part of a gene, is affected by selection, the ratio between the non-synonymous (dN) and the synonymous (dS) substitution rates can be estimated. If there is no selection acting on the sites then only random fixation should determine the rate and that would be on average the same for both dN and dS so the ratio should be around 1. According to the nearly neutral model, which states that most mutations are neutral or slightly deleterious, non-synonymous mutations should be removed by selection at a higher rate than synonymous substitutions, who evolve mostly neutrally, leading to a dN/dS << 1, negative purifying selection. In rare cases, for example after a gene duplication leading to new functions or after change in the environment for a population, mutations changing the function of a protein could be beneficial. Then the rate of fixation of non-synonymous mutations will be higher than under neutrality, dN/dS > 1.
 
dN/dS < 1	Purifying selection

dN/dS = 1	Neutral evolution, relaxed selection

dN/dS > 1	Positive selection


There are multiple programs developed test if there are signs of positive selection, here we will use *codeml* in the program suite [PaML](https://github.com/abacus-gene/paml). For more information on how to run different models [Álvarez-Carretero et al. 2023](https://doi.org/10.1093/molbev/msad041), [PaML-manual](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf), [PaML_FAQ](https://ocw.mit.edu/courses/6-877j-computational-evolutionary-biology-fall-2005/9a6d5e515fb1e7608eb3919855b01880_pamlfaqs.pdf).
Codeml uses a codon model of evolution, where the codon triplet is the unit of evolution. It is not only inferring the dN/dS ratio across the gene and across the phylogeny but also includes branch models for detecting elevated rates in specific branches compared to the background branches, and site models that allows the dN/dS to vary across the gene, specifying the likelihood of a specific codon being under positive selection. 

dN/dS is an estimator of omega (w), and so you will see all three notations in the tutorial.

This tutorial consists of two parts: 

First retriving and aligning orthologous genes, and run a basic model to infer the average evolutionary rate in a gene family. For this we will use a dataset of Ithomiini butterflies. 

For part two we will use the branch-site model A to test for signs of selection in a vertebrate gene family.


## Part one - model M0

We will first specify a simple model (model_M0) with the same average dN/dS across all branches to infer the evolutionary rate and the general statistics of the gene family in general. We do not have gene annotation for all the species yet, so we first need to find orthologous genes to be able to infer the evolutionary rates. One tool to detect and locate single copy genes is [BUSCO](https://busco.ezlab.org/busco_userguide.html) (Benchmarking Universal Single-Copy Orthologs). BUSCO uses lineage specific databases containing genes present in > 90 % of the taxa and occur as single copy in >90% of the taxa in each lineage. It is commonly used to assess the completeness of genome assemblies, and we can use our genomes as input to retrive single copy orthologs. Since these are conserved genes our hypothesis is that they will be under strong purifying selection in general.


Start by organising the directory

```bash
mkdir selection_part_1
cd selection_part_1


```

### Input
The input for codeml is a phylogenetic tree, multiple sequence alignments, and a control file with the ending .ctl. The control file to tell codeml which models and parameters to use.

### Alignment

#### Step 1: get the single copy orthologs in our genomes
Here we will use single copy orthologs detected with BUSCO in our genomes. We infer orthogroups with OrthoFinder to reduce the risk of including paralogous genes. Paralogs have a different divergence time compared to the orthologs, which per definition should have the same divergence time as the speciation event. OrthoFinder uses different modules for detecting sequence similarities and cluster genes together in orthogroups or gene families, reconstructing species trees and use phylogenetic information to distinguish between orthologs and paralogs. OrthoFinder nicely output each single copy orthogroup as a multi fasta file that we can directly use for multiple sequence alignment.

Input for orthofinder are multi-fasta files, one for each taxa that we will concatenate from the single copy sequences from the busco output.
This is the code I used run BUSCO on all species at once. Now we want nucleotide sequences so I need to add the flag --metaeuk (using a different database), and will also take much longer time.

```bash
##Do not run during the course! It takes too long time...
busco -i ../renamed_fasta \
    -l lepidoptera_odb12 \
    -m geno \
    --metaeuk \
    -o busco_output \
    -c 12

```
To save us some time and computer power we already have set of multi-fasta single copy orthologues from BUSCO from six Ithomiini butterflies from three genera (Melinaea, Mechanitis and Napeogenes) and we are using the monarch (Danaus plexippus) as outgroup.


#### Step 2: run OrthoFinder to get orthologs genes

```bash
#set up a new working directory from your working directory
pwd # check that you are in the right place 

mkdir orthofinder
cd orthofinder

# copy all the fasta files 
cp ../orthofinder/OrthoFinder/Results_Nov28/Single_Copy_Orthologue_Sequences/${MY_GENE_FAMILY}.fa ./

ls
#take a look at one of the files
less ilMecMaza1.sco.fa

# How many sequences is there in the file?

```

```bash
#run orthofinder, this time in we are using nucleotides so we must include option -d

orthofinder -f ./ -d

```

This takes a couple of minutes.

As earlier OrthoFinder will produce a large directory with many interesting files. Are there any species that have less genes assigned to orthogroups? Why do you think this is? How many shared single copy genes do we have? 



#### Step 3: Align the sequences

We will use a multi sequence aligner that specifically accounts for codons, [PRANK](https://github.com/ariloytynoja/prank-msa/tree/master), which have a lot of other useful applications. PRANK uses evolutionary information for the placement of gaps and modelling of the substitution process. It infers a guide tree from genetic distances estimated from pairwise alignments using the neighbour-joining (NJ) algorithm and then iterates the alignment using an improved guide tree estimated from the first multiple alignment. 


Prepare input sequences

```bash
# go back to the selection_part_1 folder
cd ..
mkdir prank
cd prank

# select an orthogroup, this variable will be reused
MY_GENE_FAMILY=OG0000036

#copy the orthogroup, ops change the name of the result folder!!
cp ../orthofinder/Results_Nov28/Single_Copy_Orthologue_Sequences/${MY_GENE_FAMILY}.fa ./

#check how many sequences there are
grep ">" ${MY_GENE_FAMILY}.fa

#check total length of all sequences (to divide by 7)
wc -c ${MY_GENE_FAMILY}.fa

#change name of the sequences, remove everything after species name
cut -f1 -d"_" ${MY_GENE_FAMILY}.fa > ${MY_GENE_FAMILY}_mod.fa

```

Now we can align our sequences. We will only use one iteration in this exercise, set with the flag -iterate=1, otherwise it could take too long time. 

```bash

# run the aligner
prank -d=${MY_GENE_FAMILY}_mod.fa -o=${MY_GENE_FAMILY}_codon -codon -iterate=1 -F -f='paml'

```

This will take a while, maybe time for a quick break?

```bash
# Take a look at the alignment
less ${MY_GENE_FAMILY}_codon.phy
```


### Phylogenetic tree

The program codeml also needs a guide tree to infer the rate of substitutions (it could infer one but it takes longer time and it is not made to be a tree inference program, so not recommended). We will use iqtree to infer a tree from our alignment


```bash

mkdir iqtree
cd iqtree

#convert to interleaved phy, iqtree do not like the paml format
prank -convert -d=../prank/${MY_GENE_FAMILY}_codon.best.phy -o=${MY_GENE_FAMILY}_codon.phylipi -f=phylipi

iqtree -s ${MY_GENE_FAMILY}_codon.phylipi.phy --prefix ${MY_GENE_FAMILY}

```


### Format control file
Now we have the all the data needed to run codeml, so it is time to edit the control file.


```bash

mkdir codeml
cd codeml

```


You will have a template form for the control file (.ctl) in the common repository /home/genomics/scratch/. 
This file has to specify: 

input tree

alignment file

output file

This is important to change when running the different models or genes otherwise it will overwrite the previous result.

You also have to specify your models and other parameters of interest. See below how a general control file looks like:


      seqfile = ALN            * Path to the alignment file
     treefile = TREE           * Path to the tree file
      outfile = OUT            * Path to the output file
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = NDAT           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = CLEAN              * Remove sites with ambiguity data?
		
        model = CODMOD         * Models for ω varying across lineages
	  NSsites = NSSIT          * Models for ω varying across sites
    CodonFreq = CODFREQ        * Codon frequencies
	  estFreq = ESTFREQ        * Use observed freqs or estimate freqs by ML
        clock = CLOCK          * Clock model
    fix_omega = FIXOME         * Estimate or fix omega
        omega = INITOME        * Initial or fixed omega



```bash

# Copy the .ctl template to the codeml directory

cp /home/genomics/scratch/data/comparative_genomics/selection/template.ctl ./

# we will copy the tree file and the alignment file here to for simplicity

cp ../prank/${MY_GENE_FAMILY}_codon.best.phy ./

cp ../iqtree/${MY_GENE_FAMILY}.treefile ./


```

### Run codeml

#### Model M0
First we will run with model M0 to estimate the general level of evolutionary rate for the genes of interest. This is the simplest model with fewest parameters and produces the average w across the tree for this group of genes.
Here we get estimates of branch lengths and ω under the model and basic statistics of the data, such as the divergence levels, the base composition and codon usage bias.


```bash
mkdir model_M0
cd model_M0

cp ../template.ctl codeml-M0.ctl

# take a look at the ctl file 

less codeml-M0.ctl

# Replace variable names with the 
# values needed to run the analysis 

# Set path to input and output files (../XXXX.phy, ../XXX.treefile, XXX_M0.out)
nano codeml-M0.ctl
# exit and save

# we use sed to change the file in place with the flag -i

# Set data to 1 (only one loci)
sed -i 's/NDAT/1/' codeml-M0.ctl 
sed -i 's/CLEAN/1/' codeml-M0.ctl # Remove sites with ambiguity data? Yes=1, No=0

# Specify all the sites models 
sed -i 's/CODMOD/0/' codeml-M0.ctl  # Models for ω varying across lineages, 0 is one ratio
sed -i 's/NSSIT/0/' codeml-M0.ctl   # Models for ω varying across sites, 0 is one ratio
sed -i 's/CODFREQ/7/' codeml-M0.ctl  # Codon frequencies, use mutation-selection model
sed -i 's/ESTFREQ/0/' codeml-M0.ctl  # Use observed freqs or estimate freqs by ML
sed -i 's/CLOCK/0/' codeml-M0.ctl  # Assume no clock

# Starting values to be used when the model parameters are estimated
sed -i 's/FIXOME/0/' codeml-M0.ctl  # Enables option to estimate omega yes=0, fixed=1 
sed -i 's/INITOME/0\.5/' codeml-M0.ctl # Initial or fixed omega


```

Now we are ready to run codeml!

```bash
codeml codeml-M0.ctl > codeml-M0.log

```
This redirects (>) the information from the standard out (screen) to a log file.
 
This may take a few minutes.

Take a look at the output file.

It shows the alignment before and after removing ambiguous sites, and the site pattern count. Other useful pieces of information are the length of the alignment, nucleotide composition, codon usage. Check the codon usage table, if you have a very biased codon usage (only 20), you might have used the wrong type of aligner. Some programs do protein alignment and then use a codon table to back-translate instead of using the actual nucleotide alignment.

Pairwise comparison can be used check if the divergence in the tree is appropriate for using dN/dS as test of selection. With high divergence there will be risk of substitution saturation, if too closely related they will not have accumulated enough substitutions to be informative. This can also be inspected the treelength, branch length dN and dS for the tree. 

Omega (w) is what we are interested in here.
What does the average evolutionary rate in this gene suggest?



## Part 2 test for positive selection using branch-site model A

The Model-M0 is quite unrealistic, is is rather unlikely that all branches and all sites would have the same evolutionary rate. A more realistic model is a model that allow the w to vary among branches and account for different selection pressures on codons in the gene. Branch-site model A is a test of positive selection in a proportion of sites in the foreground branch relative to the background branches. An increased rate could for example suggest local adaption to novel environment in our species of interest compare to its relatives. 

The test is indicative of positive selection if of a proportion of codons at which w > 1, and if the likelihood of this model is higher than the null model. The null model uses the same parameters except we fix omega to w = 1 instead of allowing for w > 1. The Likelihood Ratio Test (LRT) statistic is used against a chi-square distribution with 1 degrees of freedom for significance testing. The LRT is constructed to compare nested models, so a null model that does not allow for any codons with w > 1, against a more general model that does.

The alternative branch-site model A allows omega to vary in the branch of interest by setting fix_omega=0, estimate from the data.

Model A: model = 2, NSsites = 2, fix_omega = 0

The null model is also the branch-site model A but with w = 1 fixed to 1, specified by

Model A1: model = 2, NSsites = 2, fix_omega = 1, omega = 1

#### Branch-site model: Alternative model

Prepare for the analysis.

```bash
# Go back to your working directory and make a new directory
mkdir selection_part_2
cd selection_part_2
```

Copy alignment and tree file, we will use these for both models

```bash
cp /home/genomics/scratch/data/comparative_genomics/selection/input_model_A/vertebrate* ./

```

In the tree file we need to add a label to the branch or branches we want as foreground branch.

```bash

sed 's/Chicken_Mx/Chicken_Mx \#1/' vertebrate.tree > vertebrate_Chicken_Mx.tree

```

We will first run the alternative model:

```bash

mkdir model_A_est
cd model_A_est

#copy and rename the template file
cp ../../selection/codeml/template.ctl model_A_est.ctl

# change input and output paths and filenames

nano model_A_est.ctl

```
Now we need to set the parameters in the control file to

Nsites = 2
Model = 2
Fix(Omega) = 0 (estimate from data)
Omega = .4 (initial value)


```bash

# Set data to 1 (only one loci)
sed -i 's/NDAT/1/' model_A_est.ctl 
sed -i 's/CLEAN/1/' model_A_est.ctl # Remove sites with ambiguity data? Yes=1, No=0

# Specify all the sites models 
sed -i 's/CODMOD/2/' model_A_est.ctl  # Models for ω varying across lineages
sed -i 's/NSSIT/2/' model_A_est.ctl   # Models for ω varying across sites
sed -i 's/CODFREQ/7/' model_A_est.ctl  # Codon frequencies, use mutation-selection model
sed -i 's/ESTFREQ/0/' model_A_est.ctl  # Use observed freqs or estimate freqs by ML
sed -i 's/CLOCK/0/' model_A_est.ctl  # Assume no clock

# Starting values to be used when the model parameters are estimated
sed -i 's/FIXOME/0/' model_A_est.ctl  # Enables option to estimate omega yes=0, fixed=1 
sed -i 's/INITOME/.4/' model_A_est.ctl # Initial or fixed omega, start value for the model

```

Time to run the model

```bash
codeml model_A_est.ctl > model_A_est.log

```

This can take some minutes depending the size of the alignment.
Some times the program have trouble with convergence of the model, and then it could help to use the tree estimated in the simplest model M0 as guide tree. 

Check the output in the .out file.

```bash
less model_A_est_chicken.out

grep -A5 "MLE" model_A_est_chicken.out
```

Here we an see the the proportion of sites in the different classes and the estimated omega values for each of the classes.
We can also see the lnL log-likelihood for the model.


#### Branch-site model: Null model
Now we can run the alternative model for the branch-site test of positive selection.

Input variables are the same, but the output file name needs to be changed, and we need to change to fixed omega. So we force omega to be 1 in the class that we previously allowed omega to be above 1, all other model paramenters are the same.

```bash

cd ../
mkdir model_A_fixed
cd  model_A_fixed

# copy the control file and change the name of the file, we will change the output and omega settings
cp ../model_A_est/model_A_est.ctl model_A_fixed.ctl

# change the name of the output file NAME is your orthogroup (NAME_modA_fixed.out)
# Fix(Omega) = 1 (fixed omega)
# Omega = 1 (initial value, fixed to 1)

nano model_A_fixed.ctl
```
Run the null model

```bash

codeml model_A_fixed.ctl > model_A_fixed.log
```
This can take a while again.

```bash
less model_A_est_chicken.out

grep -A5 "MLE" model_A_est_chicken.out
```




#### LRT

Use the LRT to test if the different in likelihood is statistically significant. We are testing a restricted model to a more permissive with more free parameters so it should by design have a higher likelihood, but is this more than expected?
The LRT statistic is calculated 2x(lnL_est-lnL_fix). If lnL_est=-1130 and lnL_fix=-1146 is χ2 = 2×(-1130 + 1146) = 16. The degrees of freedom are k = 4-3 = 1. We can use the 'chi2' program from the PAML package to assess the significance of the LRT statistic. The LRT test statistic does not strictly follow a chi-square distribution but using a χ2 with one degree of freedom makes the test conservative. 

```bash
# get the likelihood for both model_A:s
grep "lnL" model_A*/*out
```
Calculate LRT = 2x(lnL_est-lnL_fix)

```bash
# check out chi2
chi2 

```

And you should see the following output:

```bash
Chi-square critical values

	Significance level

DF 0.9950 0.9750 0.9000 0.5000 0.1000 0.0500 0.0100 0.0010

1 0.0000 0.0010 0.0158 0.4549 2.7055 3.8415 6.6349 10.8276
2 0.0100 0.0506 0.2107 1.3863 4.6052 5.9915 9.2103 13.8155
...
```

The critical value for one-degree of freedom and significance level α=0.05 is 3.8415, so if you LRT statistic is larger than that we can reject our null model.

```bash
chi2 --help
# d.f. & Chi^2 value (Ctrl-c to break)?
# Type in d.f. (degree of freedom), in our case 1, and the LRT-statistics. Here is an example where the value is 4
1 4
# output
# df =  1  prob = 0.045500265 = 4.550e-02
# to run directly
chi2 1 4

#df =  1  prob = 0.045500265 = 4.550e-02

```

Is the LRT significant?

Check the output file:
If the LRT suggests presence of codons under positive selection in the foreground branch then we could check for the result of the Bayes empirical Bayes (BEB) method that calculates the posterior probabilities that each codon is from the site class of positive selection. 

In each line, the first column shows the site position followed by the amino acid at this site in the first sequence (this is for identification of the site in the sequence). The third column (Pr (w > 1)) shows the posterior probability for the site to be from the positive-selection class (i.e., with ω > 1).

Are there sites in this gene under positive selection? Which gene is it?






For more reading:
This tutorial was using part of the material from https://github.com/abacus-gene/paml-tutorial/tree/main 

Ref
Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118

Yang, Z. 2007. PAML 4: Phylogenetic analysis by maximum likelihood. Mol. Biol. Evol. 24:1586-1591.

For more reading on the data set in part 2 [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/).
