#script by Karin Näsvall using syntenyPlotter
#https://github.com/Farre-lab/syntenyPlotteR
#This script uses the output from minimap2 whole genome alignments.
#Creates chromosome size and chain files for syntenyPlotter and plots ribbon plot
#Chromosome names must be in the format of a number or name underscore number (name_number)
#dir organisation: intermediate/, plots/
#only keeps alignments > 100kb, if short alignments are needed change in the script

#install.packages("devtools")
#devtools::install_github("marta-fb/syntenyPlotteR")
library(syntenyPlotteR)
library(dplyr)
library(ggplot2)

rm(list = ls())


######
####If using command line and argument uncomment these (and comment away the variable lines below):
cmd_args <- commandArgs(trailingOnly = TRUE)
table_file <- cmd_args[1]
taxa1=cmd_args[2]
taxa2=cmd_args[3]

#################
#Variables if using hardcoding, change these variables:
#output from minimap2 in paf-format
# table_file <- "file.paf"
#reference (the taxa in the sixth column of the output in minimap2)
# taxa1="ref_taxa"
#query (the taxa in the first column of the output in minimap2)
# taxa2="query_taxa"

#################

#read in table
paf_table <- read.csv2(table_file, sep = "\t", header = F)[,1:12]
#add header
colnames(paf_table) <- c("query", "Qseq_length", "Qstart", "Qend", "strand", "reference", "Rseq_length", "Rstart", "Rend",  "matches", "total","MQ")

#head(paf_table)

#check distribution of alignments before filtering 
# hist(paf_table$total, breaks = 3000)
# hist(paf_table[paf_table$total > 500000,"total"])
# length(paf_table[paf_table$total > 500000,c("total")])

#filter , used 100-500 kb alignments since closely related taxa, needs to be adjusted
paf_table <- paf_table[paf_table$total > 100000 & paf_table$MQ==60,]
#paf_table <- paf_table[paf_table$total > 100000 & paf_table$MQ==60,]


paf_table$refID <- taxa1
paf_table$queryID <- taxa2

#correct column order for the chain file, switching the position of query and reference
chain_table <- paf_table[, c("reference", "Rstart", "Rend", "query", "Qstart", "Qend", "strand", "refID", "queryID", "Rseq_length", "Qseq_length")]


#change names of chr OBS customise!!!
#get default colours if chr names are just numbers
#otherwise have to assign colours manually for the plot

#removes everything before underscore
chain_table$reference <- sub(".*_", "", chain_table$reference)
chain_table$query <- sub(".*_", "", chain_table$query)


#order df after chr length in ref taxa, will also order the other taxa after first appearance in the chain file
chain_table <- 
  chain_table %>% arrange(desc(Rseq_length))


#change orientation if majority is on negative strand 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


chain_1 <- chain_table
chain_1$Rstart_init <- chain_1$Rstart
chain_1$Rend_init <- chain_1$Rend

chain_1$Qstart_init <- chain_1$Qstart
chain_1$Qend_init <- chain_1$Qend


for (i in unique(chain_table$query)) {
  if (Mode(chain_table[chain_table$query==i, "strand"])=="-") {
    print(i)
    
    chain_1[chain_1$query==i, "Qend"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qstart"]

    chain_1[chain_1$query==i, "Qstart"] <- chain_1[chain_1$query==i, "Qseq_length"] - chain_1[chain_1$query==i, "Qend_init"]
    
    chain_1[chain_1$query==i & chain_1$strand=="-", "strand"] <- "temp"
    chain_1[chain_1$query==i & chain_1$strand=="+", "strand"] <- "-"
    chain_1[chain_1$query==i & chain_1$strand=="temp", "strand"] <- "+"
  }
}

chain_table <- chain_1


#make size_file for input (and order of chr in the plot)

chr_size_q <- unique(chain_table[,c("query", "Qseq_length", "queryID")])
colnames(chr_size_q) <- c("chr", "seq_length", "ID")
chr_size_r <- unique(chain_table[,c("reference", "Rseq_length", "refID")])
colnames(chr_size_r) <- c("chr", "seq_length", "ID")
#the query chr are first and the ref are last
chr_size <- rbind(chr_size_q, chr_size_r)


#writing input files
names(chr_size) <- NULL
write.table(chr_size, file = paste("intermediate/chr_length_", taxa1, taxa2, ".txt", sep = ""), sep = "\t", row.names = F)

names(chain_table) <- NULL
write.table(chain_table, file = paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""), sep = "\t", row.names = F)

#produce the synteny figures
draw.linear(directory = "plots/", 
            output = paste("synt_", taxa1, taxa2, Sys.Date(), sep = ""), 
            paste("intermediate/chr_length_", taxa1, taxa2, ".txt", sep = ""), 
            paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""),  
            fileformat = "png", w=13, h=5)

# ggsave(filename = paste("plots/synt_", taxa1, taxa2, Sys.Date(), ".png", sep = ""), 
#        device = "png", width = 13, height = 5)
# 
