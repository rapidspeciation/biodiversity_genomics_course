### fastq format

Let's have a look at the first sequence from our raw read files which are stored in the [fastq format]([https://en.wikipedia.org/wiki/FASTQ_format). As we saw in the (https://github.com/rapidspeciation/biodiversity_genomics_course/blob/main/slide_presentations/03_Raw_sequences_and_quality_control.pdf), each DNA sequence is composed of four lines. Therefore, we need to visualize the first four lines to have a look at the information stored for the first sequence.

First we set the name of the fastq file that we will work with as the variable `FILE`. Then, we copy that file to our directory. Finally, we will examine the first 4 lines. However, we cannot just directly write `head -4 $FILE` like we might with a normal text file because the fastq file is actually compressed. It is thus a binary file which cannot just be read. Luckily, there are many commands that can directly read binary files. Instead of `cat`, we use `zcat`, instead of `grep`, we use `zgrep`. If we want to use any other command, we need to read the file with `zcat` and then pipe the output into our command of choice such as `head` or `tail`.

```shell

# First, we will make a folder to work in
mkdir fastqc
cd fastqc

# Now let's specify FILE as the name of the file containing the forward reads
FILE="wgs1_R1.fastq.gz"
cp /scratchsan/C_computacion/nr10sanger_ac/biodiversity_genomics_course/data/Heliconius/$FILE ./

# Let's have a look at the first read:
zcat $FILE | head -4
```

These are the four lines which make up the information for each read. You can learn more about what they each mean [here](https://en.wikipedia.org/wiki/FASTQ_format). Now the file is in our directory and readable, let's count the number of lines:

```shell
zcat $FILE | wc -l
```

The number of sequences is thus this number divided by 4, or we can count the number of lines starting with the header

```shell
zgrep -c "@E00" $FILE
```

We might think that we could have just counted the number of `@` - i.e. the first symbols for each header.

```shell
zgrep -c "@" $FILE
```

However, we see that this does not give us the same number. The reason is that `@` is also used as a symbol for encoding [quality scores](https://en.wikipedia.org/wiki/Phred_quality_score) and thus some quality score lines were also counted.

### Assessing read quality with fastqc

To assess the read quality, we use `fastqc` which is extremely easy to run and can be run with the name of the fastq file as the only argument. It can aslo handle gzipped files.

To get help on `fastqc`:

```shell
fastqc -h | less
```

Let's run `fastqc` on our read subsets:

```shell
fastqc $FILE
```

We should now also run fastqc on the file of reverse reads. As we do not need copies of these files in all of your personal directories, we will just write the file names with the paths.

`fastqc` allows an output directory with the `-o` flag. We will thus just work in our home directories and run `fastqc` giving the file name with its path and specifying the output folder as the current directory (i.e. `-o ./`).

```shell
# Reverse reads
fastqc -o ./ /scratchsan/C_computacion/nr10sanger_ac/biodiversity_genomics_course/data/Heliconius/wgs1_R2.fastq.gz 

```

When working with multiple individuals, it is more efficient to run FastQC on all of them using a single command with a wildcard such as `*.fastq.gz`. This tells FastQC to process all files in the directory that end with `*.fastq.gz`.

To save the results in a specific folder, you can first create the folder in your directory and then specify its path in the FastQC command.

```shell
# fastqc for all the .fastq.gz files in the folder
mkdir fastqc_results
fastqc -o ./fastqc_results/ /home/genomics/scratch/data/Heliconius/*.fastq.gz

```
If you have many individuals, checking each HTML file manually can become very tedious. A better option is to run MultiQC. This program summarizes the results from FastQC (and other tools) into a single, interactive HTML report. To run it, navigate to the folder where the .html files are located and use the following command:

```shell
cd fastqc_results
multiqc .

```

To visualize the HTML files, you need to download them from the cluster using the `scp` command on your local machine. For example:

```shell
scp genomics@toko.uncu.edu.ar:/home/genomics/scratch/users/<yourname>/fastqc/fastqc_results/multiqc_report.html ./

```

Here some [slides](https://github.com/speciationgenomics/presentations/blob/master/fastqc_interpretation.pdf) on interpreting fastqc html output.

### Hands-on: Run FastQC on the following data (RAD1.fastq.gz and RAD2.fastq.gz), interpret the results and answer these questions:

How does the mean quality score change along the sequence?
is there a warning for the per-base sequence content graphs?
is there a fail for the per sequence GC content graphs?

```shell
# Remember that you can run `FastQC` on all files with a single command by specifying the path to the directory containing the reads:
mkdir fastqc_results_RADs
fastqc -o ./fastqc_results_RADs/ /home/genomics/scratch/data/Heliconius/RADs/*.fastq.gz

```

### Challenging exercises for the bash wizards and those with extra time left

In the `RAD2.fastq.gz `there are some reads with very low GC content which likely represent reads of contaminants. Find the 10 reads with the lowest GC content and check what they are by blasting them.

Here one very condensed solution: Try to find your own solution first!
```shell
FILE=RAD2

cp /home/genomics/scratch/data/Heliconius/RADs/$FILE.fastq.gz ./

#Add GC content to each read in fastq file to check reads with highest or lowest GC contents:

zcat ${FILE}.fastq.gz | awk 'NR%4==2' | awk '{split($1,seq,""); gc=0; at=0; n=0; for(base in seq){if(seq[base]=="A"||seq[base]=="T") at++; else if(seq[base]=="G"||seq[base]=="C") gc++; else n++}; print $0,gc/(at+gc+n)*100,n}' > ${FILE}.gc

#Lowest GC content:
sort -k 2 -t " " $FILE.gc | head

#Highest GC content:
sort -k 2 -t " " $FILE.gc | tail

#Get the worst 10 sequences with all information:
zcat $FILE.fastq.gz | grep -f <(sort -k 2 -t " " $FILE.gc | tail | cut -d" " -f 1) -A 2 -B 1 > $FILE.lowGC

# Make a new fastq file with these reads:
grep -v "^--" $FILE.lowGC | gzip > $FILE.lowGC.fastq.gz
```

As a second exercise, try to generate a new file from the fastqz file containing every 1000th read. This is useful as subsampling is often needed to test software. Fastqc will take very long and a lot of memory if it needs to read in a giant file. It is thus better to subsample if you have large fastq files.

```shell
# Forward (R1) reads
zcat /home/genomics/scratch/data/Heliconius/wgs1_R1.fastq.gz | awk '{printf("%s",$0); n++; if(n%4==0){
printf("\n")}else{printf("\t")} }' | awk 'NR == 1 || NR % 1000 == 0' | tr "\t" "\n" | gzip > wgs.R1.subsampled.fastq.gz &

# Reverse (R2) reads
zcat /home/genomics/scratch/data/Heliconius/wgs1_R2.fastq.gz | awk '{printf("%s",$0); n++; if(n%4==0){
printf("\n")}else{printf("\t")} }' | awk 'NR == 1 || NR % 1000 == 0' | tr "\t" "\n" | gzip > wgs.R2.subsampled.fastq.gz &
```
