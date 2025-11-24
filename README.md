# Welcome to the Biodiversity Genomics course
## held in [Mendoza, Argentina](https://biodiversitygenomicslatam.weebly.com/) in November 2025
For the course taught in July 2024 at IKIAM University in Tena, Ecuador, see the folder [2024](2024).

This course is taught by Karin Näsvall, Nicol Rueda, Fernando Seixas, and Joana Meier from the [Wellcome Sanger Institute](https://www.sanger.ac.uk/group/meier-group/) and by Melisa Olave (CONICET, Argentina).
Some of the course material is based on the [speciation genomics course](https://speciationgenomics.github.io/) by Joana Meier and Mark Ravinet.

The course website showing logistics and useful information about Menoza is [here](https://biodiversitygenomicslatam.weebly.com/).


## Course structure
Day 0:
- [Slides](slide_presentations/00_Intro_unix.pdf) basic introduction to unix and the command line
- [Exercise](exercises/00_intro_unix_basic_practical.md) on how to use the command line
- [Exercise](exercises/00_more_unix.md) on more use of the command line
- [Exercise](exercises/00_advanced_unix_awk.md) for more advanced users (awk, variables, arrays, writing a bash script)

Day 1:
- [Slides](slide_presentations/01_Welcome_BiodiversityGenomics_introduction.pdf) introducing biodiversity genomics and sequencing technologies
- [Slides](slide_presentations/02_Summary_reads-vcf.pdf) summarising the steps from getting raw Illumina reads, filtering, aligning them to a reference and calling variants and genotypes
- [Slides](slide_presentations/03_Raw_sequences_and_quality_control.pdf) introducing the structure of Illumina reads and fastq files
- [Exercise](exercises/01_RawReadsExploration_fastqc.md) on exploring Illumina reads and visualising the quality with fastqc which uses these [input files](input_files/raw_reads)
- [Slides](slide_presentations/04_fastqc_interpretation.pdf) on interpreting fastqc output of RAD data

Day 2: TBA


Day 3: TBA


Day 4: TBA


Day 5: TBA


Advanced materials:
- Statistical phasing [Exercise](exercises/12_statistical_phasing.md)
- Extended haplotype statistics [Exercise](exercises/13_identifying_selection_with_haplotype_statistics.md)
- Identifying genes at selection peaks [Exercise](exercises/14_finding_candidate_genes.md)


## Learning more

### Physalia courses
[Physalia courses](https://www.physalia-courses.org/) offers lots of great courses on genomics, bioinformatics, and related fields. The [Speciation Genomics](https://speciationgenomics.github.io/) course was a course taught via Physalia.

### Speciation genomics course website
https://speciationgenomics.github.io/ contains many more tutorials, including topics not covered here, such as demographic modeling, haplotype-based tests for selection, advanced unix and R tutorials, simulating data with SliM, etc.

### Speciation and population genomics workshop
[Workshop website](https://evomics.org/workshops/workshop-on-population-and-speciation-genomics/2025-workshop-on-population-and-speciation-genomics-cesky-krumlov/) contains many useful tutorials that are quite advanced, including machine learning and pangenomics.

### tidypopgen (R package for population genomics):
[website](https://evolecolgroup.github.io/tidypopgen/index.html)

### Biodiversity Genomics Conference
28 October - 1 November 2024: Online, free conference on biodiversity genomics across all time zones:
https://www.biodiversitygenomicsconference.org

### For RAD or UCE data
If you have short-read data for only a subset of the genome because you used a reduced-representation technique (e.g. RAD or UCE), most of the tutorial will still be relevant. For RAD or UCE data you do not necessarily need a reference genome, unless you want to run the genome scans for finding regions with high differentiation or introgression.
If you have RAD (restriction-enzyme associated DNA) data, you can either follow the steps in our tutorial with mapping reads to a reference genome or if you do not have a reference genome, you can do a de novo assembly, i.e. make your own reference for just the RAD loci. The most widely used tool for RAD data analysis is [STACKS](https://catchenlab.life.illinois.edu/stacks). If you are working with polyploids, check out [polyRAD](https://academic.oup.com/g3journal/article/9/3/663/6026786).
If you have UCE (ultra-conserved elements) data, have a look at this [website](https://www.ultraconserved.org/) for guidance.

##
## Publications we recommend:

### Introduction to Unix
- Introduction to the Unix Command Line [Dowling, et al., 2019](Papers/Introduction_unix_command-may2019.pdf)
- Unix and Perl Primer for Biologists [Bradnam, et al 2016](Papers/Unix_Perl.pdf)

### Reviews on biodiversity genomics:
- How genomics can help biodiversity conservation [Theissinger et al. 2023](Papers/Theissinger_et_al-2023_Genomics_Conservation.pdf)
- Genomics and the origin of species [Seehausen et al. 2015](Papers/Seehausen_et_al-15-NatRevGenet.pdf)

### Publications related to the examples in the course:
- Genomics of Neotropical biodiversity indicators: two butterfly radiations with rampant chromosomal rearrangements and hybridisation [van der Heijden et al. 2024](Papers/vanDerHeijden_et_al-2025-Biorxiv.pdf)
- Genomic evidence reveals three W-autosome fusions in Heliconius butterflies [Rueda et al. 2024](Papers/Rueda_et_al-24-PlosGenetics.pdf)

### Papers or manuals for some of the tools mentioned in this course:
- Arima-HiC Mapping Pipeline [Arima Genomics](Papers/Arima_Mapping_UserGuide.pdf)
- Tutorial on ASTRAL [Mirarab, et al](Papers/Species_tree_Astral-tutorial.pdf)
- DensiTree Manual: Making sense of sets of trees [Bouckaert, et_al 2014](Papers/DensiTree_Manual.v2.2.pdf)
- Local PCA Shows How the Effect of Population Structure Differs Along the Genome [Li, et al., 2019](Papers/Li_2019_Local_PCA_Shows_How_the_Effect_of_Population_lostruct.pdf)
- The Sequence Alignment/Map format and SAMtools [Li, et al., 2019](Papers/Li_H_2009_The_Sequence_Alignment_Map_format_and_SAMtools.pdf)
- Dsuite - Fast D-statistics and related admixture evidence from VCF files [Malinsky, et al., 2020](Papers/Malinsky_2020_Dsuite.pdf)
- Sequence Alignment Map Format Specification [SAM BAM group](Papers/Sequence_Alignment_Map_Format_Specification.pdf)
- Quality Scores for Next-Generation Sequencing [Illumina](Papers/technote_Q-Scores.pdf)
- Accurate, scalable and integrative haplotype estimation (SHAPEIT4) [Delaneau, et al., 2019](Papers/Accurate_scalable_and_integrative_haplotype_estimation.pdf)
- WhatsHap: Weighted Haplotype Assembly for Future-Generation Sequencing Reads [Patterson, et al., 2015](Papers/WhatsHap.pdf)
- Modern technologies and algorithms for scaffolding assembled genomes [Ghurye, et al., 2019](Papers/Modern_technologies_and_algorithms_for_scaffolding_assembled_genomes.pdf)
- BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes [Manni, et al., 2021](Papers/BUSCO_Update_Novel_and_Streamlined_Workflows.pdf)
- OrthoFinder and recap on orthologs and paralogs [github.com/davidemms/OrthoFinder](https://github.com/davidemms/OrthoFinder)
- IQtree manual [iqtree.org](http://www.iqtree.org/doc/Home)

### Reviews on phylogenomics
- Phylogenomics and the reconstruction of the tree of life [Delsuc, et al., 2005](Papers/Delsuc_F.etal.2005_PHYLOGENOMICS.pdf)
- Phylogenetic tree building in the genomic age [Kapli, et al., 2020](Papers/Kapli_phylogenomics_in_genomic_era.pdf)
- Genetic Terminology [Elston, et al 2012](Papers/Elston_2012_Genetic_Terminology.pdf)
- Phylogenomics reveals the evolutionary timing and pattern of butterflies and moths [Kawahara, et al 2019](Papers/kawahara_etal-2019_phylogenomic_reveals_evolutionary_timing_pattern_butterflies_and_moths.pdf)
- A method for genome-wide genealogy estimation for thousands of samples [Speidel, et al., 2019](Papers/A_method_for_genome_wide_genealogy.pdf)

### Reviews on speciation
- What does Drosophila genetics tell us about speciation? [Mallet, J., 2006](Papers/Mallet,J_2006_speciation.pdf)
- Defining the speciation continuum. [Stankowski, et al., 2021](Papers/Stankowski_etal,2021_Defining_the_speactiation_continuum.pdf)

### Reviews or examples about introgression
- Phylogenomic approaches to detecting and characterizing introgression [Hibbins, et al., 2022](Papers/Hibbins_M._2021._Phylogenomic_approaches_to_detecting_and_characterizing.pdf)
- Evaluating the Use of ABBA–BABA Statistics to Locate Introgressed Loci [Martin, et al., 2015](Papers/Martin_etal_2015_Evaluating_Use_ABBA–BABA_Statistics.pdf)
- The genetic consequences of hybridization [Moran, et al., 2021](Papers/The_genomic_consequences_of_hybridization.pdf)
- How reticulated are species? [Mallet, et al., 2015](Mallet_How_reticulated_are_species)

### Reviews or examples on chromosome evolution
- Chromosomal rearrangements and speciation. [Rieseberg, et al., 2021](Papers/Rieseberg,etal.,_2001_Chromosomal_rearrangement_and_speciation.pdf)
- Chromosome evolution [Ingo Schubert, 2019](Papers/Chromosome_evolution.pdf)
- The Impact of Chromosomal Rearrangements in Speciation: From Micro- to Macroevolution [Lucek, et al., 2023](Papers/The_Impact_of_Chromosomal_Rearrangements_Speciation.pdf)

### Review on next-generation sequencies technologies
- Next-Generation Sequencing Technology: Current Trends and Advancements [Satam, et al., 2023](Papers/Satam_et_al_2024_Next_generation_sequencing.pdf)

### Population genomic diversity and divergence statistics
- Patterns of Z chromosome divergence among Heliconius species highlight the importance of historical demography [Van Belleghem, et al., 2017](Papers/Van_Belleghem_2018_Patterns_of_z_chromosome.pdf)
- Complex modular architecture around a simple toolkit of wing pattern genes [Van Belleghem, 2018](Papers/VanBelleghem_2017_Complex_modular_architecture_around_simple_toolkit_of_wing_pattern_genes.pdf)
