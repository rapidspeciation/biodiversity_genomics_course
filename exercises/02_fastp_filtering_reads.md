## Filtering Illumina reads

Sometimes Illumina adapter sequences are still present in some reads because adapters can form adapter dimers and then one of them gets sequenced or if a DNA fragment is shorter than the read length, the sequencer continues to "read-through" into the adapter at the end of the DNA fragment. In the latter case the forward and the reverse read will contain adapter sequences, which is called a "palindrome". While a full adapter sequence can be identified relatively easily, reliably identifying a short partial adapter sequence is inherently difficult. However, if there is a short partial adapter present at the end of the forward read and the beginning of the reverse read, that is a good sign for a "palindrome" sequence.
To anyone using Novaseq, or any NextSeq Illumina technology. Watch out for overrepresented polyG sequences (weirdly long sequences of GGGGGGG), particularly in the reverse reads. This is a problem of the latest Illumina instruments that use a [two-colour system](https://sequencing.qcfail.com/articles/illumina-2-colour-chemistry-can-overcall-high-confidence-g-bases/) to infer the bases. A lack of signal is called as G with high confidence. These polyG tails need to be removed or the read will not map well to the reference genome.
Reads that start or end with very low quality can be aligned better if the bad quality parts are trimmed off. We will use [fastp](https://github.com/OpenGene/fastp) to fix all of these issues. fastp can remove low quality reads, adapters and polyG tails. It even automatically detects what adapters were used.
There are also other excellent read filtering and trimming tools such as [Trimmomatic](https://speciationgenomics.github.io/Trimmomatic) or the fast tool [Ktrim](https://academic.oup.com/bioinformatics/article/36/11/3561/5803071).
fastp also generates a html file that shows the read quality before and after filtering.

```shell

# Check the options of fastp
fastp -h

# Now let's again make a folder to work in
cd ~
mkdir filteredReads
cd filteredReads

# Let's get the wgs read files:
cp /home/data/wgs_raw/wgs1.R*.fastq.gz ./

# Run fastp
fastp --in1 wgs1.R1.fastq.gz --in2 wgs1.R2.fastq.gz --out1 wgs1.R1.trimmed.fastq.gz --out2 wgs1.R2.trimmed.fastq.gz -l 50 -h wgs.html &> wgs1.log

# Note &> redirects the information on what it did into the file wgs.log (both stderror and stdout are written into this file)

# Let's have a look at what files fastp produced:
ls
# You can see that it produced the two output files we specified and also an html and a json file which shows some quality information

# To take a look at the HTML file, we first need to download it to our computer using the following command line. Make sure to use your user number and replace IP by the actual IP address.

scp -i c1.pem user1@IP:xx~/wgs.html ./
```
### Parameters specified here:
* \-\-in1 and \-\-in2: specify your files of forward (1) reads and of the reverse (2) reads.
* \-\-out1 and \-\-out2: specify the output files for forward and reverse reads that are still Paired.
* -l 50: this specifies that if a read is shorter than 50 basepairs after all filters, it should be removed.
* -h: specifies name for the html file with plots showing the read quality before and after filtering


### PolyG tail trimming
This feature removes the polyG tails that arise from lack of signal in NextSeq/NovaSeq technologies. It is enabled for Nextseq/Novaseq data by default, and you can specify -g to enable it for any data, or specify -G to disable it.

### Removal of adapter sequences
Adapter trimming is enabled by default, but you can disable it with -A. Adapter sequences can be automatically detected for both PE/SE data.

### Length filter
Reads below the length threshold (e.g. due to adapter removal) are removed. Length filtering is enabled by default. The minimum length requirement is specified with -l.

### Quality filtering
Quality filtering is enabled by default, but you can disable it with -Q. Currently fastp supports filtering by limiting the number of uncalled (N) bases (-n, Default 5) and the percentage of unqualified bases.
To filter reads by its percentage of unqualified bases, two options should be provided:
* -q : Quality threshold per base required. Default: 15, which means that a Phred quality score of at least 15 is required
* -u : Percent of bases allowed to be below the quality threshold to keep the read (0~100). Default 40 means 40% bases can fail the quality threshold. If more bases fail, the read is removed.

### Hands-on: Run fastp on the following data and interpret the results.

```shell
RAD1.fastq.gz and RAD2.fastq.gz

```
