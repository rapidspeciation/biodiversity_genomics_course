## Dsuite

We will use [Dsuite](https://github.com/millanek/Dsuite) to infer if introgression has occured in the past between non-sister species, here *H. melpomene* and *H. timareta*. In this example, we will test for introgression between sympatric *H. timareta thelxionoe* and *H. melpomene amaryllis* versus allopatric melpomene (*H. m. melpomene*) 
```shell
# Connect to the server
ssh genomics@toko.uncu.edu.ar
ssh toko05

# go to your user directory and make a folder for Dsuite
cd /home/genomics/scratch/users/<yourname>
mkdir Dsuite
cd Dsuite
```

Dsuite requires a SET file containing the names of all individuals in the first column and the names of the species or populations in the second column. The outgroup needs to be called "Outgroup". We will thus rename *H. numata* individuals to Outgroup. To do that, we will use a command called `sed` that allows replacing one word by another word in each line.

Note that the $ sign after `Hnum.bsl.bra` specifies that we only want to replace the word `Hnum.bsl.bra` if it is the last word of the line. We write the output into a new file called `Melpomene.sets.txt`. We also need to remove the header as Dsuite does not like if the file has a header (here `ind species`). For this, we will use `grep -v "ind"` which will print all lines that do not contain the word `ind`.

```shell
# prepare the SETS file
cp /home/genomics/scratch/data/martin2019/martin2019_info.txt ./
sed 's/Hnum.bsl.bra$/Outgroup/' martin2019_info.txt | grep -v "ind"  > melpomene.sets.txt
tail melpomene.sets.txt
```

*D* statistics only work with **bi-allelic** SNPs. So if you have not done it already, make sure to remove multi-allelic SNPs (i.e. positions with more than one alternative allele), indels and monomorphic sites, e.g. with plink. 

We compute D statistics with [Dsuite](https://github.com/millanek/Dsuite) with the tool `Dtrios` using the vcf file and our newly created file with the species names and outgroup.

```shell
VCF="/home/genomics/scratch/data/martin2019/wgenome.martin2019.biallelic.mac2.vcf.gz"
Dsuite Dtrios $VCF melpomene.sets.txt
```

Let's have a look at the output files. By piping it into column -t, we can align the columns so that the header and the results are nicely aligned.

```shell
cat melpomene.sets_BBAA.txt | column -t
cat melpomene.sets_Dmin.txt | column -t
```

The file with the `*_BAAA.txt` suffix orders each trio assuming that the correct tree is the one where the BBAA pattern is more common than the discordant ABBA and BABA patterns.
The file with the `*_Dmin.txt` suffix outputs the minimum D for each trio regardless of any assumptions about the tree topology.

If we had more trios, we could parallelise the analysis using `DtriosParallel` (see the [Dsuite website](https://github.com/millanek/Dsuite) for instructions on how to do this).

If you want to explore ***D* statistics** more, I would recommend using the `admixr` R-package. Here a [tutorial](https://speciationgenomics.github.io/ADMIXTOOLS_admixr/). To infer the direction of gene flow, I recommend [**Dfoil**](https://github.com/jbpease/dfoil). If you have many species that might have hybridised, check out **Fbranch**, which is part of Dsuite and allows the visualisation of Dstatistics across many different species comparisons. It could also be useful to run an **ADMIXTURE** or **STRUCTURE** plot in order to figure out if gene flow is still ongoing. Here a [tutorial](https://speciationgenomics.github.io/ADMIXTURE/). If gene flow is ongoing, ADMIXTURE will show that some individuals are introgressed. However, if gene flow is ancestral, only *D* statistics will show it.
