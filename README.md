# Welcome to the Biodiversity Genomics course, held for the first time at IKIAM, in Tena, Ecuador!
## held for the first time at IKIAM, in Tena, Ecuador in July 2024
This course taught by Karin Näsvall, Nicol Rueda and Joana Meier is an adapted version of the [speciation genomics course](https://speciationgenomics.github.io/) by Mark Ravinet and Joana Meier. Here you will find all relevant course materials for the biodiversity genomics course, including all files we used during the course in case you want to redo anything at home.

Day 1:
- [Slides](slide_presentations/01_Welcome_BiodiversityGenomics_introduction.pdf) introducing biodiversity genomics and sequencing technologies
- [Slides](slide_presentations/02_Summary_reads-vcf.pdf) summarising the steps from getting raw Illumina reads, filtering, aligning them to a reference and calling variants and genotypes
- [Slides](slide_presentations/03_Raw_sequences_and_quality_control.pdf) introducing the structure of Illumina reads and fastq files
- [Slides](exercises/Connecting_to_the_Amazon_server.pdf) on logging on to the Amazon server by Carlo Pecoraro from [Physalia Courses](https://www.physalia-courses.org)
- [Exercise](exercises/01_RawReadsExploration_fastqc.md) on exploring Illumina reads and visualising the quality with fastqc which uses these [input files](input_files/raw_reads)
- [Slides](slide_presentations/04_fastqc_interpretation.pdf) on interpreting fastqc output of RAD data

Day 2:
- [Exercise](02_fastp_filtering_reads.md) on filtering Illumina reads with fastp which uses these [input files](input_files/raw_reads)
- [Slides](slide_presentations/06_Aligning_reads_to_reference.pdf) on aligning reads to a reference genome
- [Exercise](03_Mapping_to_a_reference_genome.md) to align Illumina paired-end reads to a reference genome with bwa-mem2
- [Slides](slide_presentations/07_Variant_and_genotype_calling.pdf) about variant and genotype calling

Day 3:
- [Exercise](exercises/04_variant_calling.md) on variant and genotype calling with bcftools.
- [Exercise](exercises/05_filtering_variants.md) on filtering vcf files.
- [Slides](slide_presentations/08_Ithomiini_introduction_PCA.pdf) introducing ithomiini butterflies and PCA
- [Exercise](exercises/06_pca.md) on PCA.
- [Exercise](exercises/07_iqtree.md) on making phylogenies with iqtree.

Day 4:
- [Slides](slide_presentations/09_Detecting_hybridisation_Dstats.pdf) on identifying introgression with D statistics
- [Exercise](exercises/08_Dstatistics.md) on Dstatistics to infer introgression using Dsuite
- [Slides](slide_presentations/10_Genome_scans.pdf) on genome scans
- [Exercise](exercises/09_genome_scan.md) on genome scans
