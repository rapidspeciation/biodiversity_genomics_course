# Dsuite
We will use [Dsuite](https://github.com/millanek/Dsuite) to infer if introgression has occured in the past between non-sister species, here Mechanitis nesaea and M. lysimnia.

```shell
cd
mkdir Dsuite
cd Dsuite
```
Dsuite requires that the outgroup is called "Outgroup". We will thus rename messenoides to Outgroup. To do that, we will use a command called `sed` that allows replacing one word by another word in each line. Note that the $ sign after messenoides specifies that we only want to replace the word messenoides if it is the last word of the line. We write the output into a new file called Mechanitis.outgroup.txt. We also need to remove the header as Dsuite does not like if the file has a header (here `ind species`). For this, we will use `grep -v "ind"` which will print all lines that do not contain the word `ind`.

```shell
cp
sed 's/messenoides$/Outgroup/' Mechanitis.info | grep -v "ind"  > Mechanitis.outgroup.txt
```

Dstatistics only work with bi-allelic SNPs, we will thus remove multi-allelic SNPs (=positions with more than one alternative allele), indels and monomorphic sites with vcftools.

```shell
vcftools --gzvcf ~/Share/Mechanitis/Mechanitis.vcf.gz --mac 3 --remove-indels --max-alleles 2 --recode --stdout | gzip > Mechanitis.snps.vcf.gz
```

Now we can run Dsuite with the tool Dtrios using the vcf file and our newly created file with the species names and outgroup.
```shell
Dsuite Dtrios Mechanitis.snps.vcf.gz Mechanitis.outgroup.txt
```


Let's have a look at the output files. By piping it into `column -t`, we can align the columns so that the header and the results are nicely aligned.

```shell
cat Mechanitis.outgroup_BBAA.txt | column -t
```

If you want to explore D statistics more, I would recommend using the admixr R-package. Here a [tutorial](https://speciationgenomics.github.io/ADMIXTOOLS_admixr/). To infer the direction of gene flow, I recommend [Dfoil](https://github.com/jbpease/dfoil). It could also be useful to run an ADMIXTURE or STRUCTURE plot in order to figure out if gene flow is still ongoing. Here a [tutorial](https://speciationgenomics.github.io/ADMIXTURE/). If gene flow is ongoing, ADMIXTURE will show that some individuals are introgressed. However, if gene flow is ancestral, only D statistics will show it.
