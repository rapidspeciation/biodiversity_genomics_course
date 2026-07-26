# Aligning reads to the reference genome

Now that we have our filtered reads, we need to map (or align) them to the reference genome.

Using alignment software, we essentially find where in the genome our reads originate from and then once these reads are aligned, we are able to either call variants or construct a consensus sequence for our set of aligned reads.

### Getting access to the reference genome

We will be aligning our sequence data to the *Heliconius sara* reference genome, first published by [Rueda-M *et al.* (2024)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011318).

### Hands-on: Downloading the genome reference

On the National Center for Biotechnology Information (NCBI) website, search for the reference genome of the species *Heliconius sara* v1.2 and download it directly to the home directory of the cluster using the wget command:

```shell
mkdir reference
cd reference/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/917/862/395/GCA_917862395.2_iHelSar1.2/GCA_917862395.2_iHelSar1.2_genomic.fna.gz
```

The file is compressed with `gzip`, so before we can do anything with it, we need to decompress it. Usually we avoid decompressing files but it was compressed here because of it's large size and we need it to be uncompressed for it to be used properly in our analysis. So to do this we simply use `gunzip`:

```shell
gunzip GCA_917862395.2_iHelSar1.2_genomic.fna.gz
```

In order to align reads to the genome, we are going to use `bwa` which is a very fast and straightforward aligner. See [here](http://bio-bwa.sourceforge.net/) for more details on `bwa`.

#### Indexing the reference genome

Before we can actually perform an alignment, we need to index the reference genome we just copied to our home directories. This essentially produces an index for rapid searching and aligning. We use the `bwa-mem2 index` tool to achieve this.

```shell
bwa-mem2 index GCA_917862395.2_iHelSar1.2_genomic.fna
```
The `bwa-mem2 index` tool simply requires the reference fasta file from which to build our genome index. So it is a very simple command.

If it takes too long, you can copy the file as shown below.

```shell
cp ~/Share/reference/GCA_* ./
```

Use `ls` to take a look, but this will have copied in about 5 files all with the `GCA_917862395.2_iHelSar1.2_genomic.fna.` prefix that we will use for a reference alignment.

When `bwa-mem2` aligns reads, it needs access to these files, so they should be in the same directory as the reference genome. Then when we actually run the alignment, we tell `bwa-mem2` where the reference is and it does the rest. To make this easier, we will make a variable pointing to the reference.

```shell
REF=~/reference/GCA_917862395.2_iHelSar1.2_genomic.fna
```

#### Performing a paired end alignment

Now we are ready to align our sequences! To simplify matters, we will first try this on a single individual. First, we will create a directory to hold our aligned data:

```shell
cd ~
mkdir align
cd align
```
As a side note, it is good practice to keep your data well organised like this, otherwise things can become confusing and difficult in the future.

To align our individual we will use `bwa`. You might want to first have a look at the options available for it simply by calling `bwa`. We are actually going to use `bwa mem` which is the best option for short reads.

We will use the individual - `wgs` which we have already trimmed. There are two files for this individual - `R1` and `R2` which are forward and reverse reads respectively.

Let's go ahead and align our data, we will break down what we did shortly after. Note that we run this command from the home directory.

```shell
bwa-mem2 mem -t 4 $REF \
~/Share/filteredReads/wgs.R1.trimmed.fastq.gz \
~/Share/filteredReads/ wgs.R2.trimmed.fastq.gz > wgs1.sam
```
Since we are only using a shortened fastq file, with 100K reads in it, this should just take a couple of minutes. In the meantime, we can breakdown what we actually did here.

`-t` tells `bwa` how many threads (cores) on a cluster to use - this effectively determines its speed.

Following these options, we then specify the reference genome, the forward reads and the reverse reads. Finally we write the output of the alignment to a **SAM file**.

Once your alignment has ended, you will see some alignment statistics written to the screen. We will come back to these later - first, we will learn about what a SAM file actually is.

Since these analyses take quite a bit of time, we will stop the analysis and copy the output file from the Share folder.

    cp ~/Share/align/wgs1.sam ./

#### SAM files

Lets take a closer look at the output. To do this we will use `samtools`. More details on `samtools` can be found [here](http://www.htslib.org/doc/samtools.html)

```shell
cd align
samtools view -h wgs1.sam | head
samtools view wgs1.sam | head
```

This is a SAM file - or sequence alignment/map format. It is basically a text format for storing an alignment. The format is made up of a header where each line begins with `@`, which we saw when we used the `-h` flag and an alignment section.

The header gives us details on what has been done to the SAM, where it has been aligned and so on. We will not focus on it too much here but there is a very detailed SAM format specification [here](https://samtools.github.io/hts-specs/SAMv1.pdf).

The alignment section is more informative at this stage and it has at least 11 fields. The most important for now are the first four. Take a closer look at them.

```shell
samtools view wgs1.sam | head | cut -f 1-4
```

Here we have:

* The sequence ID
* The flag - these have various meanings, 0 = mapped, 4 = unmapped
* Reference name - reference scaffold the read maps to
* Position - the mapping position/coordinate for the read

**Note!!** The other fields are important, but for our purposes we will not examine them in detail here.

Now, lets look at the mapping statistics again:

```shell
samtools flagstat wgs1.sam
```

This shows us that a total of 600k reads were read in (forward and reverse), that around 95% mapped successfully, 84% mapped with their mate pair, 1.15% were singletons and the rest did not map.

#### BAM files

One problem with SAM files is that for large amounts of sequence data, they can rapidly become HUGE in size. As such, you will often see them stored as BAM files (Binary Aligment/Mapping format). A lot of downstream analyses require the BAM format so our next task is to convert our SAM to a BAM. Note, there is an even more space-efficient format is CRAM (typically 30-60% smaller than BAM, but can only be read if the reference genome is also provided).

```shell
samtools view wgs1.sam -b -o wgs1.bam
```

Here the `-b` flag tells `samtools` to output a BAM and `-o` identifies the output path. Take a look at the files with `ls -lah` - you will see a substantial difference in their size. This would be even more striking if we were to map a full dataset.

You can view bamfiles in the same way as before.

```shell
samtools view wgs1.bam | head
```
Note that the following will not work (although it does for SAM) because it is a binary format

```shell
head wgs1.bam
```

Before we can proceed to more exciting analyses, there is one more step we need to perform - sorting the bamfile. Most downstream analyses require this in order to function efficiently

```shell
samtools sort wgs1.bam -o wgs1_sort.bam
```

Once this is run, we will have a sorted bam. One point to note here, could we have done this is a more efficient manner? The answer is yes, actually we could have run all of these commands in a single line using pipes like so:

```shell
bwa-mem2 mem -t 4 $REF ~/Share/wgs_raw/wgs.R1.fastq.gz ~/Share/wgs_raw/wgs.R2.fastq.gz | samtools view -b | samtools sort -T wgs_sort > ./align/wgs_sort.bam
```

However as you may have noticed, we have only performed this on a single individual so far... what if we want to do it on multiple individuals? Do we need to type all this everytime? The answer is no - we could do this much more efficiently.

### Hands-on: Map the individual 2 (wgs2.R1.fastq.gz wgs2.R2.fastq.gz)

### Working with multiple individuals

There are many ways to work with multiple individuals. One way is to use a bash script to loop through all the files and align them one by one. We will examine this way in detail together. An alternative way is to use some form of parallelisation and actually map them all in parallel (=at the same time). We will also demonstrate an example of this, but it is quite advanced and is really only to give you some familiarity with the approach. An even more advanced way is to use workflow managers like [Snakemake](https://snakemake.readthedocs.io/en/stable/).

#### Using a bash script to loop through and align individuals

We are going to work together to write a bash script which you can use to map multiple individuals one by one. It is typically a lot more straightforward to write a script locally on your machine and then upload it to the cluster to make sure it does what you want.

The first thing we will do is initiate the script with the line telling the interpreter that the file must be treated as a bash script:

```shell
#!/bin/sh
```
Next, we will declare an array to ensure that we have all our individuals

```shell
INDS=($(for i in ~/Share/wgs_raw/*R1.fastq.gz; do echo $(basename ${i%.R*}); done))
```
This will create a list of individuals which we can then loop through in order to map each individual. Here we used bash substitution to take each forward read name, remove the directory and leave only the individual name.

If you want to, decleare the array in your command line and then test it (i.e. type `echo ${IND[@]}`). You will see that we have only individual names, which gives us some flexibility to take our individual name and edit it inside our `for` loop (i.e. it makes defining input and output files much easier.

Next we will add the actual `for` loop to our script. We will use the following:

```shell
for IND in ${INDS[@]};
do
	# declare variables
	REF=~/reference/GCA_917862395.2_iHelSar1.2_genomic.fna
	FORWARD=~/Share/wgs_raw/${IND}.R1.fastq.gz
	REVERSE=~/Share/wgs_raw/${IND}.R2.fastq.gz
	OUTPUT=~/align/${IND}_sort.bam

done
```
In this first version of our loop, we are making the `$REF`, `$FORWARD`, `$REVERSE` and `$OUTPUT` variables, making this much easier for us to edit later. It is a good idea to test the variables we are making in each loop, using the `echo` command. Feel free to test this yourself now - you can just copy and paste from your script into the command line to make it easier.

After we have tested the loop to make sure it is working properly, all we have to do is add the `bwa mem` command we made earlier but with our declared variables in place.

```shell
for IND in ${INDS[@]};
do
	# declare variables
	REF=~/reference/GCA_917862395.2_iHelSar1.2_genomic.fna
	FORWARD=~/Share/wgs_raw/${IND}.R1.fastq.gz
	REVERSE=~/Share/wgs_raw/${IND}.R2.fastq.gz
	OUTPUT=~/align/${IND}_sort.bam

	# then align and sort
	echo "Aligning $IND with bwa"
	bwa-mem2 mem -t 4 $REF $FORWARD \
	$REVERSE | samtools view -b | \
	samtools sort -T ${IND} > $OUTPUT

done
```

With this completed `for` loop, we have a command that will create variables, align and sort our reads and also echo to the screen which individual it is working on, so that we know how well it is progessing.

Although it takes a bit more work to make a script like this, it is worth it. This script is very general - you would only need to edit the variables in order to make it work on almost any dataset on any cluster. Reusability (and clarity) of scripts is something to strive for.

Now that we have built our script, it is time to use it. Save it and name it `align_sort.sh` We will move it on to the cluster, either using `scp` or `filezilla`, a process explained in [this tutorial](). Be sure to move it to your `home` folder.

Once the script is on the cluster, open a screen and call it `align`; (i.e. `screen -S align`). Then run the script like so:

```shell
bash align_sort.sh
```

You will now see the script running as the sequences align. Press `Ctrl + A + D` in order to leave the screen. Now is a good time to take a break as you wait for the job to complete.

Since these analyses take quite a bit of time, we will stop the analysis and copy the output files from the Share folder to run the next step:

```shell
cp ~/Share/align/*.bam ./
```

#### Remove duplicates reads

Finally, we need to remove duplicate reads from the dataset to avoid PCR duplicates and technical duplicates which inflate our sequencing depth and give us false certainty in the genotype calls. We can use [Picard Tools](https://broadinstitute.github.io/picard/) to do that.

```shell
java -Xmx1g -jar /home/scripts/picard.jar \
 MarkDuplicates REMOVE_DUPLICATES=true \
 ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
 INPUT={}.bam \
 OUTPUT={}.rmd.bam \
 METRICS_FILE={}.rmd.bam.metrics

# Now we need to index all bam files again and that's it!
samtools index *.rmd.bam
 ```
