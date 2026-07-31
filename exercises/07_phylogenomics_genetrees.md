# Phylogenomics: from raw data to species trees


## **Objective**

In this practice, you will learn the main steps to perform phylogenetic inference with genomic data. These steps include obtaining genomic data (for this practice it will be from public databases), performing quality control, assembling sequences from short-read data (Illumina), extracting the loci of interest, aligning sequences, inferring gene trees, and building a species tree with support based on quartet frequency.

---



## **1. Downloading sequence data**

To obtain short-read sequence data, we will use **fasterq-dump**, a tool from the **SRA Toolkit** package.

`SRA Toolkit` is a set of tools provided by the National Center for Biotechnology Information (NCBI) that allows downloading and manipulating data from the Sequence Read Archive (SRA). This program is essential for obtaining the raw sequences needed for the analyses proposed in this practice.

NCBI is an institution that provides access to biological databases, including the genomic and transcriptomic sequence repository. We will use their portal to access public data through the SRA Toolkit. SRA is a format defined by NCBI for NGS data.

### **Steps:**

1.1. **Create an environment with SRA Toolkit** (including pigz and parallel) and activate it:

   ```bash
   conda create -n sra-tools -c bioconda sra-tools pigz parallel

   conda activate sra-tools
   ```

1.2. **Download the FASTQ files**

- List of selected accessions:
     
Check the article [Gardner et al. (2023)](https://doi.org/10.1073/pnas.2222035120) and look for where the accession code information is found. Select 10 species and make a list with the accession codes of the selected samples and their respective IDs separated by tabs. To make it easier to explore and select the species to analyze, use the [ENA accession search tool](https://www.ebi.ac.uk/ena/browser/home)

Save this information in a text file called `accesiones.txt`.

   - Expected structure of the text file `accesiones.txt`

It must contain two columns separated by a tab and end with an empty line at the end:
```
SRR24706287 Ficus_apollinaris
SRR24706125 Ficus_austrocaledonica
SRR24706157 Ficus_callosa
SRR24706179 Ficus_assimilis
SRR24706382 Ficus_ingens
SRR24706212 Ficus_platypoda
SRR24706402 Ficus_globosa
SRR24706401 Ficus_gommelleira
SRR24706366 Ficus_lutea
SRR24706331 Ficus_antandronarum
```
    Column 1: Sequence accession number.
    Column 2: Desired name for the downloaded files.


   - Download the sequences using this [script](https://github.com/gsilvaarias/curso-sistematica-biologica/blob/main/download_fastq.sh). Download the script file and save it in your working directory.

Once you have downloaded the script, run it with the following command:

```bash
bash download_fastq.sh
```

This process will generate `.fq.gz` files for each downloaded accession with the appropriate names for `Captus` (`_R1` and `_R2`) and the desired IDs for each sample (species name).

Once the process of downloading the sequences to your machine is finished and you have verified that the `.fq.gz` files are complete in the folder, deactivate the `sra-tools` environment:

```bash
conda deactivate
```


---

### **Install software for data process and analysis**

For this part we can create a conda environment that includes all needed software with this command:

```
conda create -n phylogenetic -c conda-forge -c bioconda captus bbmap=39.52 salmon=1.10.3 mafft=7.526 iqtree aster -y

conda activate phylogenetic
```

---

## **2. Quality control, assembly, extraction, and gene alignment with Captus**

`CAPTUS` is an automated pipeline designed to de novo assemble target gene sequences from raw FASTQ data, especially useful in phylogenomics studies. It integrates tools for quality control, assembly, ortholog/paralog extraction, and multiple alignment, facilitating a reproducible workflow from obtaining/downloading raw data to obtaining phylogenetic matrices ready for phylogenetic analysis. It is compatible with data from hybridization-based capture libraries (target enrichment), transcriptomes, or whole genome sequencing, adapting to studies ranging from population scale to deep phylogenetics (Ortiz et al. 2023).

- Full documentation: [https://edgardomortiz.github.io/captus.docs/](https://www.google.com/url?q=https%3A%2F%2Fedgardomortiz.github.io%2Fcaptus.docs%2F)


### **2.1 Quality control of reads**

Next, we proceed with the quality control of the raw reads. To do this, use the `clean` module of `CAPTUS`, which automates filtering, trimming, and quality evaluation of sequences in `.fastq.gz` format.

To run the process, make sure to correctly indicate the path of the folder containing the downloaded and compressed raw read files (`00_raw_reads`). The module will automatically identify the corresponding file pairs for each sample (suffixes `_R1.fastq.gz` and `_R2.fastq.gz`) and apply the default or user-defined quality filters.

```bash
captus clean -r 00_raw_reads --threads 4
```

This will generate clean `.fq.gz` files in the `01_clean_reads` folder. Examine the report of the cleaning and quality control process for **all samples** by opening the `captus-clean_report.html` file in your browser (you will find it inside the `01_clean_reads` folder).

Additionally, select a specific sample and analyze in detail the individual `FastQC` reports generated before and after filtering. These are found in the `01_clean_reads/01_qc_stats_before/` and `01_clean_reads/02_qc_stats_after/` folders. Within these, locate the `fastqc_report.html` files corresponding to the selected species, before and after cleaning, for both the forward (`R1`) and reverse (`R2`) reads. Evaluate the changes in read quality, the presence of adapters, and the length distribution to determine the effectiveness of the cleaning process applied by `CAPTUS`.

* To guide your interpretation of the results, consult the [documentation of the clean module](https://edgardomortiz.github.io/captus.docs/assembly/clean/report/).


### **2.2 De novo assembly**

We use **Captus assemble** to assemble the sequences into contigs:

```bash
captus assemble -r 01_clean_reads --sample_reads_target 1_000_000 --threads 4
```

This process will perform the de novo assembly based on a maximum of 1 million reads (this subsampling is done in this exercise to reduce computational time). This process will generate assembled files in the `02_assemblies` folder. Examine the report of the assembly process by opening the `captus-extract_report.html` file in your browser (you will find it inside the `02_assemblies` folder).


### **2.3 Extraction of genes of interest**

We use **Captus extract** to recover the genes of interest using reference sequences that you can find in the file `artocarpus_333genes.fasta`, which were obtained in the study by [Gardner et al. (2016)](https://doi.org/10.3732/apps.1600017). Download the file and save it in your working directory.

```bash
captus extract -a 02_assemblies -d artocarpus_333genes.fasta --threads 4
```

This will generate FASTA files in the `03_extractions` folder with the sequences extracted from the set of assembled contigs for each of the selected samples and each of the reference loci. Examine the report of the extraction process by opening the `captus-extract_report.html` file in your browser (you will find it inside the `03_extractions` folder).


### **2.4 Gene alignment**

We use **Captus align** to align the extracted sequences:

```bash
captus align -e 03_extractions --threads 4
```

This will generate alignments in the `04_alignments` folder. Examine the report of the alignment process by opening the `captus-align_report.html` file in your browser (you will find it inside the `04_alignments` folder).


---

## **3. Inference of gene trees with IQ-TREE**

For each of the aligned loci, we will use **IQ-TREE** to infer gene phylogenetic trees.


### **3.1 Obtaining gene trees**

With the alignments obtained for the loci, perform a phylogenetic inference analysis using the Maximum Likelihood approach implemented in `iqtree3`. Since this analysis involves repeatedly performing the same process on the fasta (`.fna`) files (alignments), you can automate the process using a `loop` implemented in the script `run_iqtree.sh` available [here](https://github.com/gsilvaarias/curso-sistematica-biologica/blob/main/run_iqtree.sh):

Once you have the script saved in your working directory, run it in the terminal with the command:

```bash
bash run_iqtree.sh
```

This will generate a new folder called `05_trees` with the output files of the phylogenetic analysis for each locus. Explore the content of the different output files.


---

## **4. Species tree inference with ASTRAL-Pro 3**

We will use `ASTRAL-Pro 3` to build a **species tree** with support based on quartet frequency.


### **4.1. Obtaining the species tree with `ASTRAL-Pro 3`**

Once you have inferred the gene trees with `IQ-TREE 3`, you can estimate the **species tree** using `ASTRAL-Pro 3`. This program takes into account **gene tree discordance** and **gene duplication events (paralogs)**, making it well suited for phylogenomic datasets.

To run `ASTRAL-Pro 3`, you need **two input files**:

1. **A single file containing all gene trees** in Newick format, generated with `IQ-TREE 3`.

   Create a directory called `06_species_tree` and concatenate all gene trees into a single file:

```bash
mkdir 06_species_tree
cat 05_trees/*.treefile > 06_species_tree/gene_trees.tre
```

2. **A gene-to-species mapping file** named `mapping.txt`.

   This file tells `ASTRAL-Pro 3` **which gene copies belong to each species**. It is required because phylogenomic datasets often contain **paralogs**, meaning that a species may be represented by multiple sequences in the same gene tree. For example, datasets generated with `Captus` typically label duplicated copies with the species names (or any sample names you use in your input fastq files) and suffixes such as `_00`, `_01`, `_02`, etc.

   To generate this file, download the script `get_tree_tips.py` and place it in your working directory. The script reads all tip names from the gene trees and automatically creates the required species-to-gene correspondence table (`mapping.txt`).

> **Important:** `ASTRAL-Pro 3` cannot correctly interpret multiple gene copies without this mapping file.

Before running the analysis, make sure the following files are present in your **working directory**:

* `run_astral-pro3.sh`. Available [here](https://github.com/gsilvaarias/curso-sistematica-biologica/blob/main/run_astral-pro3.sh)
* `get_tree_tips.py`. Available [here](https://github.com/gsilvaarias/curso-sistematica-biologica/blob/main/get_tree_tips.py)

Once both scripts are in place, run:

```bash
bash run_astral-pro3.sh
```

> **Expected output:** The script will generate the species tree inferred by `ASTRAL-Pro 3` in the file `species_tree.tre` located in the `06_species_tree` directory.


### **4.2. Tree visualization and analysis of node support based on quartet frequency**

Finally, draw the obtained species tree and include on each node a way to represent the support values obtained. You can do this using this R code:

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("treeio", "ggtree", "ggimage"))


library(treeio)
library(ggtree)
library(ggimage)

tree <- read.astral("species_tree.tre")

# root the phylo part directly
phy <- as.phylo(tree)
phy_rooted <- root(phy, outgroup = "Ficus_apollinaris", resolve.root = TRUE)

# reattach the ASTRAL data (q1/q2/q3 etc.) by node number
tree_rooted <- treeio::full_join(as.treedata(phy_rooted), as_tibble(tree), by = "node")

# plot the tree with quartet frequencies on nodes
pies <- nodepie(as_tibble(tree_rooted), cols = c("q1","q2","q3"))
ggtree(tree_rooted) %<+% as_tibble(tree_rooted) +
  geom_tiplab() +
  geom_inset(pies, width = 0.05, height = 0.05)

```



---

## **5. References**

- Ortiz EM, Höwener A, Shigita G, Raza M, Maurin O, Zuntini A, Forest F, Baker WJ and Schaefer H (2024). A novel phylogenomics pipeline reveals complex pattern of reticulate evolution in Cucurbitales. bioRxiv. DOI: [https://doi.org/10.1101/2023.10.27.564367](https://www.google.com/url?q=https%3A%2F%2Fdoi.org%2F10.1101%2F2023.10.27.564367)

- Nguyen L, Schmidt HA, von Haeseler A, Minh BQ (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution. DOI: [https://doi.org/10.1093/molbev/msu300](https://www.google.com/url?q=https%3A%2F%2Fdoi.org%2F10.1093%2Fmolbev%2Fmsu300)

- Zhang C, Rabiee M, Sayyari E and Mirarab S (2020). ASTRAL-Pro: Quartet-based species-tree inference despite paralogy. Molecular Biology and Evolution. DOI: [https://doi.org/10.1093/molbev/msaa139](https://www.google.com/url?q=https%3A%2F%2Fdoi.org%2F10.1093%2Fmolbev%2Fmsaa139)



