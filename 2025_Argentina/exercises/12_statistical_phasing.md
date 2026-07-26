## Statistical phasing with shapeit2
(Note, this tutorial is an adapted version of the speciation genomics tutorial by Mark Ravinet and Joana Meier: https://speciationgenomics.github.io/phasing)

Many tools require haplotype information, e.g. extended haplotype statistics we are running in this tutorial, but also tools that do local ancestry inference like Relate, tskit or fineSTRUCTURE.

Phasing basically means figuring out for heterozygous positions which of the alleles are part of the same haplotype or chromosome (e.g. which of the alleles were inherited together on a chromosome from the mother). For instance if an individual is heterozygous at two SNPs with genotypes AG and TC, phasing would tell us if the allele A at SNP1 one was inherited on the same chromosome like T or like C at the second SNP. Phased genotypes would be given as A\|G and C\|T meaning that A and C are on the same chromosome (e.g. maternal) and G and T are on the same chromosome (e.g. paternal).

In the absence of long reads that span these SNPs, we can use statistical phasing using all sequenced individuals of the same population (the more the better). There are lots of different tools for phasing and most of them also impute missing genotypes. This means that they infer missing genotypes statistically resulting in a dataset without missing data.

The most accurate way to phase is to sequence [**parent-offspring trios**](https://genome.cshlp.org/content/23/1/142.full.html). However this is simply not possible for a large number of study systems. As an alternative, we can instead [statistically phase](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3217888/) using computational methods. If you have a linkage map, it is recommended to use it for making phasing more accurate as it accounts for recombination rate variation.

There are several different packages available for phasing such as [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) and [fastPHASE](http://scheet.org/software.html). For this tutorial, we will use [shapeit2](http://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html). Note that there are more [recent versions of shapeit](https://jmarchini.org/shapeit3/) available but these typically only offer significant improvements if you have very large sample sizes (i.e. 100s or 1000s of individuals).

Shapeit2 does not require recombination rate information.

Our aim in this tutorial is to phase whole genome resequencing data from house sparrows (*Passer domesticus domesticus*) and a wild subspecies, the Bactrianius sparrow (*Passer domesticus bactrianus*). The dataset is taken from [Ravinet *et al* (2018)](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2018.1246). We will ultimately use this phased data to perform a selection scan using extended haplotype statistics.

### Getting the vcf file to the directory where you want to work

```shell
# Let's make a working directory and move into it
mkdir phasing
cd phasing
# Get the file to work with. It is a small vcf of sparrow genotypes
cp ~/workshop_materials/28_ehh/data/sparrows_chr8.vcf.gz ./
```

Like for most tools, bcftools requires an indexed file, so let's index our vcf file first:
```shell
bcftools index sparrows_chr8.vcf.gz
```

Now, let's get to know this file. Let's see how many individuals and SNPs we have here.

```shell
bcftools query -l sparrows_chr8.vcf.gz | wc -l  # number of individuals
bcftools query -l sparrows_chr8.vcf.gz   # list of individuals
bcftools index -n sparrows_chr8.vcf.gz   # number of sites in this file
```

### Phasing the data

As we mentioned previously, `shapeit` is a tool for [statistically estimating haplotype phase from genotypes](https://en.wikipedia.org/wiki/Haplotype_estimation). Phasing is actually quite easy to do - what is much harder is knowing whether you did it correctly...


With this taken care of, we are ready to run `shapeit`. Note that here we are running the program only for a single chromosome - this is the most practical (and sensible) way to phase - it makes it much easier to parallelise.

We do this like so:

```shell
shapeit --input-vcf sparrows_chr8.vcf.gz \
 -O sparrows_chr8_phased \
--window 0.5 -T 2
```

Shapeit complains that there are individuals with too high missing data proportion. Shapeit will impute missing genotypes, i.e. replace the missing genotypes by the best guess given the nearby genotypes of that individual. Too much missing data is an issue as there will be more and more genotypes that need to be inferred and not much information to use. Figure out which individual has >10% missing data and remove it. If you don't remember how to indentify the missing data per individual, you can look it up in Q3 of the genomics [primer tutorial](https://evomics.org/learning/population-and-speciation-genomics/2025-population-and-speciation-genomics/population-genomics-primer-i/#ex1.1). To remove it, there are many ways. Here one suggestion with vcftools:

```shell
vcftools --gzvcf sparrows_chr8.vcf.gz --remove-indv <add individual name> \
  --recode --stdout | bgzip > sparrows_chr8_v2.vcf.gz

# Perhaps not best practice, but let's now overwrite the original file.
mv sparrows_chr8_v2.vcf.gz sparrows_chr8.vcf.gz
```

Now let's run Shapeit again. As it complains again, we will add --force because we know that the missing data is low enough.

```shell
shapeit --input-vcf sparrows_chr8.vcf.gz \
 -O sparrows_chr8_phased \
--window 0.5 -T 2 --force
```

What did we do with this command?

* `--input-vcf` - this flag allows us to read in the vcf we want to phase.
* `--window` - sets the window size within which the phasing is carried out. The default is 1 Mb, here we set it to 0.5 Mb.
* `-T` - set the number of threads to use for parallelisation - we used 2 here


### Examining the phased data and converting it to a phased vcf

Once `shapeit` has finished running we can see it has produced two files, `sparrows_chr8_phased.haps` and `sparrows_chr8_phased.sample`. We can look at these in more detail:

```shell
head sparrows_chr8_phased.sample
```

The `sparrows_chr8_phased.sample` file is simply a list of samples with the ID for each sample repeated in two columns and a third column showing the proportion of missing data. Since there is no missing data in this dataset, this column is just full of zeros.

We can also look at the `.haps` output.

```shell
less -S sparrows_chr8_phased.haps
```

This is basically a matrix with the first five columns identical to those in a vcf - i.e. chromosome, ID, position, reference allele, alternative allele. After this, each entry is the phased allele for each individual, where `0` is the reference allele and `1` is the alternative.

If we want to use this data downstream, particularly in the `R` package `rehh`, we need to convert it to a vcf. Luckily this is easy using `shapeit`:

```shell
shapeit -convert \
--input-haps sparrows_chr8_phased \
--output-vcf sparrows_chr8_phased.vcf
```

We can then compress and index the vcf, like so:

```shell
bgzip sparrows_chr8_phased.vcf
bcftools index sparrows_chr8_phased.vcf.gz
```

Before moving on to the next step, let's have a quick look at our phased vcf to see how differs from an unphased vcf file.

```shell
bcftools view -H sparrows_chr8.vcf.gz | head | cut -f 1-12
bcftools view -H sparrows_chr8_phased.vcf.gz | head | cut -f 1-12
```

Comparing the two, we can see that the phased vcf only contains the genotypes - all the other information has been stripped out. Furthermore, in the phased vcf, genotypes are encoded as `0|0`, `0|1` or `1|1` instead of `0/0`, `0/1` or `1/1`. This is because in a vcf, `|` is typically used to denote that the phase of these loci are known. Thus, it is possible to read the haplotype of an individual by reading downwards across loci on either side of the `|`.

### Subsetting the vcf for a selection scans

We will be using our phased vcf for long-range haplotype statistic estimation in `rehh`. While it is possible to split the populations in the vcf apart in R, it is a bit clumbsy to do so. Instead, it is easier to split the vcf using `bcftools`. To do this, we first need the sample names

```shell
# look at the sample names
bcftools query -l sparrows_chr8_phased.vcf.gz
# extract sample names
bcftools query -l sparrows_chr8_phased.vcf.gz > samples
```

This vcf contains data from two house sparrow subspecies, the European house sparrow (from France and Norway) and the wild Bactrianius sparrow (from Kazahkstan and Iran). We split the population data like so:

```shell
grep "house" samples > house # samples starting with 8 are from Norway, samples with F are from France
grep "bac" samples > bac # samples starting with P are from Iran, samples with K from Kazahkstan
```
Let's check that it worked by counting the number of lines in each set. There should be around 20 in both of them.
```shell
wc -l house
wc -l bac
```

Next we split the vcf:

```shell
bcftools view -S house -O z -o house_chr8.vcf.gz sparrows_chr8_phased.vcf.gz
bcftools view -S bac -O z -o bac_chr8.vcf.gz sparrows_chr8_phased.vcf.gz
```

What did we do here? We used the `bcftools view` command to extract the samples for each population. The `-S` flag extracts the samples listed in each file. `-O z` specifies that we want a compressed vcf. Finally `-o` tells `bcftools` where to write the output. Now all we need to do is index the vcfs.

```shell
bcftools index house_chr8.vcf.gz
bcftools index bac_chr8.vcf.gz
```

Now we're ready to read this into `R` for a selection scan analysis. However, in order to make the phasing faster, we have so far worked with a very small dataset. For the extended haplotype statistics to work, we require all SNPs on chr8.
