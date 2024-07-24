---
title: "Synteny"
layout: archive


---


# Synteny


## Intro
The aim of this tutorial is to characterise and visualise large scale rearrangements in genome structure using chromosome-level genomes in fasta format. This tutorial will show one way to visualise the macro synteny in a comparison between two genomes. There are many different ways to arrive at the same result.

We have two new *Mechanitis* genomes, and we want to know if they have konserved karyotype or if there have been any chromosomal rearrangments since the divergence. Are there any chromosomes with rearrangments? Are there any chromosomes with conserved macrosynteny?

## Whole genome alignment with Minimap2
A whole genome alignment provides detection of large scale  rearrangements between relatively closely related genomes with sufficient sequence similarity. [Minimap2](https://github.com/lh3/minimap2) is a well documented long-read aligner, which works well for comparing taxa with divergence below ~15%. Minimap2 is memory hungry (5 -25GB) but is very fast.

**Input**: Two chromosome-level genomes in fasta format

**Output**: The output can be either [sam](https://samtools.github.io/hts-specs/SAMv1.pdf) (Sequence Alignment/Map format) format or [paf](https://github.com/lh3/miniasm/blob/master/PAF.md) (Pairwise mApping Format) format depending on what you need for downstream analysis.


```
minimap2 ref.fa query.fa > approx-mapping.paf
#or for sam format
minimap2 -a ref.fa query.fa > approx-mapping.sam
```
Here we will choose the paf-format so that we easliy can visualise it in R.


### 1 - Get assemblies

```
#log in and begin in you home directory on the cloud
cd

#first organise the directory
mkdir synteny
cd synteny
mkdir genomes
cd genomes

```
The genomes can be accessed on ncbi, Mechanitis mazaeus and Mechanitis messenoides. A little repetition on how to download genomes: Go to ncbi (https://www.ncbi.nlm.nih.gov/datasets/genome/) and search for *Mechanitis*. We are looking for the primary assemblies of *Mechanitis mazaeus* and *M. messenoides*. Click on one of them, then chose FTP, right-click on the link to the file that ends in genomic.fna.gz and use `wget <link_to_file>`in the command line to retrieve the file to the Amazon cloud.

```
#download the files from ncbi (go to https://www.ncbi.nlm.nih.gov/datasets/genome/)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/959/347/395/GCA_959347395.1_ilMecMaza1.1/GCA_959347395.1_ilMecMaza1.1_genomic.fna.gz

#Check that the file looks ok
zcat GCA_959347395.1_ilMecMaza1.1_genomic.fna.gz | head

#check how many chromosomes and scaffolds there are by grepping the fasta header, which always starts with a >
#very important to use quotes around the >
zcat GCA_959347395.1_ilMecMaza1.1_genomic.fna.gz | grep ">"

```
Now do the same for *M. messenoides*. Do you think there are any rearrangements between these genomes?

To make it easy for ourselves in the downstream analyses we will rename the fasta sequences in the genome files. There are several software available to manipulate fasta files [seqkit](https://bioinf.shenwei.me/seqkit/) and [seqtk](https://github.com/lh3/seqtk) to name two. Here we will use seqkit and rename the chromosomes and create a copy of the genome file ending with \_renamed.fa. Here I am using the abbreviation in the genome file as species identifier (ilMecMaza1).

```
#make variables, genome without gz
GENOME=GCA_959347395.1_ilMecMaza1.1_genomic.fna
TAXA_NAME=ilMecMaza1

#make a tab separated key-value file with the old seq name and the new seq name
zcat ${GENOME}.gz | grep ">" | awk -v taxa_name=$TAXA_NAME '{print $0"\t"taxa_name"_"$7}' | tr -d ">" > ../list_chr_names_$TAXA_NAME.txt

#run seqkit
zcat ${GENOME}.gz | seqkit replace -p "(.+)" -r '{kv}' -k ../list_chr_names_$TAXA_NAME.txt - > ${GENOME%.*}_renamed.fa

#check the fasta headers, remember the qoutes ">"!
grep ">" ${GENOME%.*}_renamed.fa
```
Options in `seqkit`
`-p` pattern
`-r` replacement
`'{kv}'` key - value


Do the same for *M. messenoides*.

```
#go back to the synteny folder
cd ../
```

Now we should be ready to run Minimap2.

### 2 - Run Minimap2

```
#check that you are in your synteny folder
pwd
#if not go there
cd ~/synteny/
#make a directory for minimap
mkdir minimap
cd minimap
#make directories for the output and for log-files
mkdir output log
```

Check that minimap2 is installed and working:
```
minimap2 --help
```

As you can see there are many options to refine the mapping and use different optimisations depending on the data. Running Minimap2 with default settings and no options results in an approximate alignment coordinates without base level detail, but usually enough for detecting and visualising large scale rearrangements.


```
minimap2 -t 2 ../genomes/GCA_959347395.1_ilMecMaza1.1_genomic_renamed.fa ../genomes/GCA_959347415.1_ilMecMess1.1_genomic_renamed.fa > output/MecMaza_MecMess.paf

```

This takes 10-20 minutes on a large cluster, but here it might take much longer and Minimap2 is memory demanding. With the settings -t 2 it took 15 min on our cluster with maximum memory usage of approximately 10 GB. During the wait you can take a look at Step 3 visualisation, or get a coffee.

If the script is not finished or memory demands are to high we have prepared results in the `Share` folder. Copy the result file `MecMaza_MecMess.paf` to you output directory and look at the output:
```shell
cp ~/Share/synteny/MecMaza_MecMess.paf output/
head output/MecMaza_MecMess.paf

```
This file have a lot of columns. The important ones for us right now are the query and ref sequence names and length for each alignment, start and end of alignment, strand, gaps, total aligned bases, mapping quality. To interpret the output this site is informative: [paf](https://github.com/lh3/miniasm/blob/master/PAF.md).
It is hard to directly interpret the output from the table so we need additional tools for exploring the results.

### 3 - Visualise in R SyntenyPlotteR
There are many different ways to visualise genome alignments, such as dotplots, chromosome painting, ribbon plots, circos plots etc. This is reflected in the number of different packages and scripts available for alignment exploration and plotting, like the R-packages  [pafr](https://github.com/dwinter/pafr), [gggenomes](https://thackl.github.io/gggenomes/), [Rideogram](https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html), [circlize](https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#customize-chromosome-track) and [plotsr](https://github.com/schneebergerlab/plotsr) and [pycircilise](https://moshi4.github.io/pyCirclize/) for Python just to name a few.

Here we will use [SyntenyPlotteR](https://github.com/Farre-lab/syntenyPlotteR) to create a ribbon plot.

All synteny visualisation require information about the name and size of the chromosomes and the links between the genomes (start and end of alignments or positions of common markers). SyntenyPlotter uses two intermediate files, a chromosome size file and a chain file containing the connections betwen the two genomes. We will use a wrapper script that  create these files and makes the ribbon plot.

```shell
#organise our directory
#go back to synteny folder
cd ~/synteny/

mkdir syntenyplotter
cd syntenyplotter
#create two files, one for intermediate files and one for the resultant plots
mkdir intermediate plots

```
Copy the Syntenyplotter_paf_wrapper.R to your `syntenyplotter` directory.
```shell
#copy the script from the Share folder
cp  ~/Share/synteny/Syntenyplotter_paf_wrapper.R ./
```

This script formates the output from minimap2 to fit the input of SyntenyplotteR, prints the two intermediate files in the folder `intermediate`, then plots and saves the figure in the `plots` directory. You can run this script on the commandline.

It requies three arguments:
1. the paf-file (including the path to it)
2. a name for the reference taxa
3. a name for the query taxa.

The order of the arguments is important in this case.

```
#copy the whole syntenyplotter folder to your local computer
scp -r -i user1.pem  user1@35.161.175.22:~/synteny/syntenyplotter/ ./
#go into the syntenyplotter folder
cd syntenyplotter/
#check
ls
#run the script
Rscript Syntenyplotter_paf_wrapper.R ../minimap/output/MecMaza_MecMess.paf ilMecMaza1 ilMecMess1
```
Your alignment plot should now be in your `plots` directory. Describe what you see. Can you answer some of the questions we asked in the beginning?

If you want you can add the variables directly in the script (hardcode). The hardcoding could be useful if running your own genomes with a different formatting of the sequence (chromosome) names, the script is taylored for chromosome names in the format string_number or just a number.

**Extra**: Sometimes you need to refine to plot to increase the visibilty of the rearrangements to facilitate interpretation. For example, I want the Z-chromosome to be displayed last in both taxa. Try and change the order of the chromosomes in the plot.


## Synteny plot with BUSCO genes
For more divergent taxa whole genome alignment can be difficult due to low sequence similarity and low mappability. However, genes with important and basal functions are usually conserved due to strong puryfing selection. We can use this to find anchor points or common markers in the genomes we want to compare. Are the genes in the same order and on the same chromsome in the genomes? Or have rearrangements changed their positions? If we have an annotation there are several tools for synteny analyses like McscanX and GeneScape. If we do not have an annotation we can still find conserved genes. One tool to detect and locate conserved genes is [BUSCO](https://busco.ezlab.org/busco_userguide.html) (Benchmarking Universal Single-Copy Orthologs). BUSCO uses lineage specific databases containing genes present in > 90 % of the taxa and occur as single copy in >90% of the taxa in each lineage. It is commonly used to assess the completeness of genome assemblies, but we can use our genomes as input to retrive the position of busco genes in the genomes of interest.

### 4 - BUSCO

```
#go back to the synteny folder
cd ../   #you can always check with pwd where you are
pwd
mkdir busco
cd busco

#check that busco is installed and works
busco --help

#check version
busco -v
```
BUSCO is widely used and frequently updated, if you want to rerun with new genomes it is important to use the same version of BUSCO and the same version of the lineage database, or rerun all genomes with the new version.

We will use the same *Mechanitis* genomes as in the previous excercise so they should already be present in the `genomes` folder.

You can run BUSCO either on each genome separately or for a batch of genomes by specifying the input directory.

Example:
```
#usage: busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]
busco -i $INPUT_DIR \          #the input file or directory with the input files
    -l lepidoptera_odb10 \   #the lineage database closest to the taxa of interest
    -o $OUTPUT_DIR \            #output directory for the result files
    -m geno \                     #mode of BUSCO, genomes, proteins or transcriptomes
    -c 12                            #number of threads
```

Run it on one of our genomes:
If we do not specify an output directory BUSCO will create an output directory in the directory where you are running the command.
```
busco -i ../genomes/GCA_959347395.1_ilMecMaza1.1_genomic_renamed.fa -l lepidoptera_odb10 -m geno
```

BUSCO takes some time and is memory consuming so you can start one genome to see that it works but if it takes to long we have prepared output files that you can copy to your output directory.

Take a look at the output.
Check the file short_summary.txt.
```
less ~/Share/synteny/busco_out_summary/short_summary.txt
```
Is the quality of the assembly good enough to use for synteny analysis?

One result of interest for our purposes are the full_table.tsv with the genomic position of the best hits of the genes in the database.

```
head ~/Share/synteny/busco_out_summary/full_table.tsv
```
BUSCO also output the nucleotide and protein sequences of the potential single copy genes. We will use the single_copy_busco_sequences in the next step.

### 5 - OrthoFinder
We will use the single copy sequences from BUSCO to get the orthogroups from [orthofinder](https://github.com/davidemms/OrthoFinder ). OrthoFinder uses different modules for detecting sequence similarities and cluster genes together in orthogroups or gene families, reconstructing species trees and use phylogenetic information to distinguish between orthologs and paralogs.

Input for orthofinder are multi-fasta files, one for each taxa that we will concatenate from the single copy sequences from the busco output.
```
#start from the synteny folder
cd ~/synteny/
mkdir orthofinder
cd orthofinder

#check that orthofinder works
orthofinder -h
```

Do not run the next section, we already have the files prepared in the `Share/synteny`folder. But I am adding the commands here that we used to create these files, so you can see how it was done.

```shell
#make a list of the busco runs that we want to include
ls ../busco/ |grep "fa" > list_fasta.txt
#did we get the right files
less list_fasta.txt

for file in $(cat list_fasta.txt)
do
cat ../busco/${file}/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences/*.faa > ${file%.*}_sco_cat.fa
wait
done

#check that we have the files
ls
```
This is what we will do:
```
#we will copy the sequences from the Share folder
cp ~/Share/synteny/busco/busco_sco_cat/*fa ./

#what do they look like
less ilMecMaza1.1_sco_cat.fa

#run orthofinder
orthofinder -f ./ -t 2
```

This will run for a couple of minutes.

Explore the output, there is a lot of information here. Take a look at the comparative statistics. Here you can see number of genes in orthogroups, number of unassigned genes etc. This is not so interesting now when we are using BUSCOs but if you have a genome annotation this can tell you of total number of genes, genes clustering in orthogroups, number of species specific genes, etc.
Note: exchange the path to your actual result folder.
```
ls OrthoFinder/

#exchange the path to your actual result folder.
less OrthoFinder/Results_Jul23/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv
```

What we are interested in now are the single copy orthologues and their position in each genome.

This file shows the name of orthogroups with single copy orthologues. We can count how many they are.
```
#count the number of seq
wc -l OrthoFinder/Results_Jul23/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
```
Are these markers enough to detect the  rearrangements we are interested in?
You can do a rough estimate of marker density per MB by dividing the number of markers with the genome size.


### 6 -  Visualise in R circlize
This time we will use a circular graph to show the chromosomes and the links between them. Circos was originally developed in a perl script, but now there are both R and Python versions. we will use [circlize](https://jokergoo.github.io/circlize_book/book/introduction.html), which is an R implementation.
### Prepare input
Here we will use the output from OrthoFinder, we want to use single copy orthologues so we get the position from /Orthogroups.tsv by selecting only those that are present as singel copy in the file Orthogroups_SingleCopyOrthologues.txt. we can use grep with a file of patters with the -f option.

```
grep -f OrthoFinder/Results_Jul23/Orthogroups/Orthogroups_SingleCopyOrthologues.txt OrthoFinder/Results_Jul23/Orthogroups/Orthogroups.tsv > single_copy_orthogroups.tsv
```
Most synteny visualisations requires a chromosome size file and we do not get this information in the single_copy_markers.tsv, so we have to create one. Here we will create a file with chromosome name and length from samtools [faidx]( http://www.htslib.org/doc/samtools-faidx.html).
```shell
#go to the folder with the genomes
cd ../genomes
#check that samtools works
samtools
#run samtools faidx to creat and index file
samtools faidx GCA_959347415.1_ilMecMess1.1_genomic_renamed.fa
```
This command will output a tab-separated index file genome.fa.fai, with the name and length of each sequence (chromosome) in the fasta file. In addition, information about the offset of the postion in the file along with line length in bases and bytes are also given (column 3-5). Here we only care about the first and second column.

```shell
head GCA_959347415.1_ilMecMess1.1_genomic_renamed.fa.fai
```
Do the same for *M. messenoides.*

Now we have prepared the input files so we are ready to visualise our genomes. We will do this in Rstudio on our local computer. We will prepara a folder with all the input files that we can copy to our local computer.

#create a new folder
```shell
#go back to synteny
cd ../
mkdir circlize
cd circlize
#copy the script
cp ~/Share/synteny/circlize_orthofinderR.Rmd ./
#copy the chromosome length files
cp ../genomes/*.fai ./
#copy the link file
cp ../orthofinder/single_copy_orthogroups.tsv ./
#check that we have our files
ls

#make a directory for the plot
mkdir plots
```

**Local computer**

We will copy the `circlize` directory to our local computer.
```
#make sure you are in your workshop dir
#copy the folder from the Amazon cluster
scp -r -i  user1.pem  user1@35.161.175.22:~/synteny/circlize/ ./
```
Open the script in Rstudio.
We will go through the script step by step.

```r
#load libraries
library(circlize)
library(dplyr)
library(tidyr)
library(gtools)
library(forcats)
```

Add variables:
```r
#ref taxa
TAXA_1 <- "ilMecMaza1"
#for plotting
NAME_REF <- "Mechanitis mazaeus"
#chromosome length file from samtools faidx
REF_FAI <- "GCA_959347395.1_ilMecMaza1.1_genomic_renamed.fa.fai"

#query TAXA
TAXA_2 <- "ilMecMess1"
NAME_QUERY <- "Mechanitis messenoides"
QUERY_FAI <- "GCA_959347415.1_ilMecMess1.1_genomic_renamed.fa.fai"


CHAIN_FILE <- "single_copy_orthogroups.tsv"

```
Read in the data.
```r
#read in the data
sco.df <- read.csv(CHAIN_FILE, sep = "\t", header = F)

#check
sco.df

#convert to long format and split columns
sco_long.df <-
    tidyr::pivot_longer(sco.df, cols = c(V2:V5)) %>%
    separate(value, into = c("marker", "position", "strand"), sep = "\\|") %>%
    separate(position, into = c("taxa_id", "chr", "start_end"), sep = "_") %>%
    separate(start_end, into = c("start", "end"), sep = "-", convert = T) %>%
    select(-c(V1, name))
```
Now we need to prepare the input for circos.

```r
sco_long.df$seq_id <- paste(sco_long.df$taxa_id, sco_long.df$chr, sep = "_" )


chr_length.df <- rbind(read.table(REF_FAI)[,1:2],
                       read.table(QUERY_FAI)[,1:2])
colnames(chr_length.df) <- c("seq_id", "chr_length")


#select taxa and filter
map_seq <-
  sco_long.df %>%
  filter(taxa_id %in% c(TAXA_1, TAXA_2))
map_seq <- as.data.frame(map_seq)

#remove na
map_seq <-
  map_seq %>%
  na.omit()

#keep only genes present in both
map_seq <- subset(map_seq, marker %in% Reduce(intersect, split(map_seq$marker, map_seq$taxa_id)))

#order file after ref so markers are in the same order in the chain files
map_seq <-
  map_seq %>%
  mutate(taxa_id, taxa_id=as.factor(taxa_id)) %>%
  arrange(taxa_id, forcats::fct_inorder(marker))

#create chain files, keep the order so that the lines will end up in right place
chain_circ_ref <- map_seq[map_seq$taxa_id==TAXA_1, c("seq_id", "start", "end", "marker") ]
chain_circ_query <- map_seq[map_seq$taxa_id==TAXA_2, c("seq_id", "start", "end", "marker") ]

#order chromosomes and taxa, map_seq table determines order, optional
#reorder after chr number
sub_ref <- map_seq[map_seq$taxa_id==TAXA_1,]
sub_ref <- sub_ref[mixedorder(sub_ref$seq_id),]

#reorder after reverse chr number in query
sub_query <- map_seq[map_seq$taxa_id==TAXA_2,]
sub_query <- sub_query[mixedorder(sub_query$seq_id, decreasing = T),]

#merge again with query first to ordered map_seq
map_seq <- rbind(sub_query, sub_ref)

#add chr length
map_seq <- left_join(map_seq, chr_length.df)

```
Set colours

```r
#set the colour vector after ref chromosomes
synt_col <- rev(as.vector(c(pals::kelly(),pals::polychrome()))[1:c(length(unique(sub_ref$seq_id)))])

#change grey
synt_col[13] <- "#2B4E00"
synt_col[14] <- "#222222"


#make df
synt_col.df <- data.frame(seq_id=c(unique(sub_ref$seq_id)), synt_col=synt_col)

#merge with the ref df
chain_circ_ref <- left_join(chain_circ_ref, synt_col.df)

#text colour
col_text <- "black"

#colour of the chr blocks representing the chromosomes
block_col <- alpha(c(rep("goldenrod1", length(unique(map_seq[map_seq$taxa_id==TAXA_2, "seq_id"]))),
                   rep("royalblue1", length(unique(map_seq[map_seq$taxa_id==TAXA_1, "seq_id"])))),
                   alpha = 0.8)
#link colours
anc_col <- alpha(chain_circ_ref$synt_col, alpha = 0.2)
```

Finally time for the actual plotting
```r
#run circos
#make the image
pdf(paste("plots/circos", taxa_1, taxa_2, Sys.Date(),".pdf", sep = ""))
#png(paste("plots/circos", taxa_1, taxa_2, Sys.Date(),".png", sep = ""))

circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0))
circos.initialize(factors=unique(map_seq$seq_id),
xlim=matrix(c(rep(0, length(unique(map_seq$seq_id))), unique(map_seq$chr_length)),
            ncol=2))

#The xlim matrix defines the start and stop for each genome/chr. Essentially the genome or chr sizes.

circos.track(ylim=c(0, 1), panel.fun=function(x, y) {
chr=gsub(".*_", "", CELL_META$sector.index)
xlim=CELL_META$xlim
ylim=CELL_META$ylim
circos.text(x=mean(xlim), y=mean(ylim), labels = chr,
            cex=0.6, col=col_text, facing="bending.inside", niceFacing=TRUE)
}, bg.col=block_col, bg.border=F, track.height=0.06)

circos.text(sector.index = paste(taxa_1,(round(length(unique(chain_circ_ref$seq_id))/2,0) -1), sep = "_"),
            x=0,y=0,adj=c(0.05,-2.4),labels=name_ref,facing="bending.inside", font = 3)
circos.text(sector.index = paste(taxa_2,(round(length(unique(chain_circ_query$seq_id))/2,0) -2), sep = "_"),
            x=0,y=0,adj=c(0.3, 3),labels=name_query,facing="bending.outside", font = 3)

# rearrangements
circos.genomicLink(chain_circ_ref[,1:4], chain_circ_query, col=anc_col)

dev.off()
```
Your alignment plot should now be in your `plots` directory. Take a look, and describe what you see. Can you answer some of the questions we asked in the beginning? Are there any differences compared to the whole genome alignment we did with Minimap2/SyntenyPlotteR?
