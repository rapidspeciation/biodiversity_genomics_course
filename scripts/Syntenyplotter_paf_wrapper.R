#script by Karin Näsvall using syntenyPlotter
#https://github.com/Farre-lab/syntenyPlotteR
#This script uses the output from minimap2 whole genome alignments.
#Creates chromosome size and chain files for syntenyPlotter and plots ribbon plot
#Chromosome names must be in the format of a number or name underscore number (name_number)
#dir organisation: intermediate/, plots/
#only keeps alignments > 100kb, if short alignments are needed change in the script

#install.packages("devtools")
#devtools::install_github("marta-fb/syntenyPlotteR")
library(dplyr)
#library not needed, added the function
#library(syntenyPlotteR)
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

# function
draw.linear <- function(output, sizefile, ..., directory = NULL, fileformat = "png", colours = colours.default, w = 13, h = 5, opacity = .5) {

  if (is.null(directory)) {
    directory <- tempdir()
  }

  synteny.data.reframing <- function(data, tar.y, ref.y, compiled.size) {
    synteny <- data.frame()
    for (i in c(1:nrow(data))) {
      reference <- data[i, "ref.species"]
      target <- data[i, "tar.species"]
      tar_chr <- data[i, "tarchr"]
      ref_chr <- data[i, "refchr"]
      dir <- data[i, "dir"]
      tar_sizes <- compiled.size[compiled.size$species == target, ]
      names(tar_sizes) <- c("tarchr", "size", "species", "xstart", "xend")
      ref_sizes <- compiled.size[compiled.size$species == reference, ]
      names(ref_sizes) <- c("refchr", "size", "species", "xstart", "xend")
      tar_add <- tar_sizes[as.character(tar_sizes$tarchr) == as.character(tar_chr), ]$xstart
      ref_add <- ref_sizes[as.character(ref_sizes$refchr) == as.character(ref_chr), ]$xstart
      tar_y <- tar.y
      ref_y <- ref.y
      tar_xstart <- data[i, "tarstart"] + tar_add
      tar_xend <- data[i, "tarend"] + tar_add
      ref_xstart <- data[i, "refstart"] + ref_add
      ref_xend <- data[i, "refend"] + ref_add

      inverted <- grepl("-", dir, fixed = TRUE)
      if (inverted == TRUE) {
        df <- data.frame(
          x = c(tar_xstart, tar_xend, ref_xstart, ref_xend), y = c(tar_y, tar_y, ref_y, ref_y),
          fill = ref_chr, group = paste0("s", i), ref = reference, tar = target
        )
      } else {
        df <- data.frame(
          x = c(tar_xstart, ref_xstart, ref_xend, tar_xend), y = c(tar_y, ref_y, ref_y, tar_y),
          fill = ref_chr, group = paste0("s", i), ref = reference, tar = target
        )
      }
      synteny <- rbind(synteny, df)
    }
    return(synteny)
  }

  colours.default <- c(
    "1" = "#BFD73B", "2" = "#39ACE2", "3" = "#F16E8A",
    "4" = "#2DB995", "5" = "#855823", "6" = "#A085BD",
    "7" = "#2EB560", "8" = "#D79128", "9" = "#FDBB63",
    "10" = "#AFDFE5", "11" = "#BF1E2D", "12" = "purple4",
    "13" = "#B59F31", "14" = "#F68B1F", "15" = "#EF374B",
    "16" = "#D376FF", "17" = "#009445", "18" = "#CE4699",
    "19" = "#7C9ACD", "20" = "#84C441", "21" = "#404F23",
    "22" = "#607F4B", "23" = "#EBB4A9", "24" = "#F6EB83",
    "25" = "#915F6D", "26" = "#602F92", "27" = "#81CEC6",
    "28" = "#F8DA04", "29" = "peachpuff2", "30" = "gray85", "33" = "peachpuff3",
    "W" = "#9590FF", "Z" = "#666666", "Y" = "#9590FF", "X" = "#666666",
    "LGE22" = "grey", "LGE64" = "gray64",
    "1A" = "pink", "1B" = "dark blue", "4A" = "light green",
    "Gap" = "white", "LG2" = "black", "LG5" = "#CC99CC"
  )

  xstart <- xend <- refchr <- tarchr <- x <- y <- group <- fill <- chromosome <- species <- NULL
  sizes <- utils::read.delim(sizefile, header = FALSE) # to be consistent with naming in EH
  names(sizes) <- c("chromosome", "size", "species")
  sizes$size <- as.numeric(gsub(",", "", sizes$size))

  count <- 0
  compiled.size <- data.frame()
  for (i in unique(sizes$species)) {
    size.intermediate <- sizes[sizes$species == i, ]
    for (x in c(1:nrow(size.intermediate))) {
      if (x == 1) {
        total_start <- 1
        total_end <- size.intermediate[x, "size"]
      } else {
        total_start <- total_end + 6000000
        total_end <- total_start + size.intermediate[x, "size"]
      }
      size.intermediate[x, "xstart"] <- total_start
      size.intermediate[x, "xend"] <- total_end
    }
    compiled.size <- rbind(compiled.size, size.intermediate)
  }

  for (z in unique(compiled.size$species)) {
    compiled.size$y[compiled.size$species == z] <- count
    count <- count + 2
  }

  list.of.files <- list()
  for (i in list(...)) {
    list.of.files[[i]] <- i
  }

  listsynt <- list()
  for (i in 1:length(list.of.files)) {
    num <- i
    file <- list.of.files[[num]]
    dataTMP <- utils::read.delim(file, header = FALSE)
    data2 <- dataTMP[, c(4, 5, 6, 1, 2, 3, 7, 8, 9)]
    colnames(data2) <- c("tarchr", "tarstart", "tarend", "refchr", "refstart", "refend", "dir", "ref.species", "tar.species")
    data2$tarstart <- as.numeric(gsub(",", "", data2$tarstart))
    data2$tarend <- as.numeric(gsub(",", "", data2$tarend))
    data2$refstart <- as.numeric(gsub(",", "", data2$refstart))
    data2$refend <- as.numeric(gsub(",", "", data2$refend))
    reference <- data2[1, "ref.species"]
    target <- data2[1, "tar.species"]
    ref_y <- compiled.size[compiled.size$species == reference, "y"]
    tar_y <- compiled.size[compiled.size$species == target, "y"]
    if (tar_y[1] > ref_y[1]){
      ref_y <- ref_y[1] + 0.1
      tar_y <- tar_y[1]
    } else{
      ref_y <- ref_y[1]
      tar_y <- tar_y[1] + 0.1
    }
    x <- synteny.data.reframing(data2, tar_y, ref_y, compiled.size)
    x$fill <- as.factor(x$fill)
    listsynt[[i]] <- x
  }

  compiled.size$chromosome <- as.factor(compiled.size$chromosome)

  p <- ggplot2::ggplot()

  for (i in 1:length(listsynt)) {
    data <- listsynt[[i]]
    reference <- data[1, "ref"]
    target <- data[1, "tar"]
    ref_sizes <- compiled.size[compiled.size$species == reference, ]
    tar_sizes <- compiled.size[compiled.size$species == target, ]
    p <- p + ggplot2::geom_rect(
      data = ref_sizes, mapping = ggplot2::aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10, fill = chromosome),
      color = "black", alpha = 0.85, linewidth = 0.2
    ) +
      ggplot2::geom_text(data = ref_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.2, label = chromosome), size = 2, angle = 45) +
      ggplot2::geom_text(data = ref_sizes, mapping = ggplot2::aes(x = 2, y = y, label = species), size = 3, hjust = 1) +
      ggplot2::geom_rect(
        data = tar_sizes, mapping = ggplot2::aes(xmin = xstart, xmax = xend, ymin = y, ymax = y + 0.10), fill = "grey85",
        color = "black", alpha = 0.85, linewidth = 0.2
      ) +
      ggplot2::geom_text(data = tar_sizes, ggplot2::aes(x = (xstart + xend) / 2, y = y + 0.2, label = chromosome), size = 2, angle = 45) +
      ggplot2::geom_text(data = tar_sizes, mapping = ggplot2::aes(x = 2, y = y, label = species), size = 3, hjust = 1) +
      ggplot2::geom_polygon(data = data, alpha = opacity, ggplot2::aes(x = x, y = y, group = group, fill = fill))
  }

  p <- p + ggplot2::scale_fill_manual(values = colours) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  message(paste0("Saving linear image to ", directory))
  print(p)
  ggplot2::ggsave(paste0(directory,"/",output, ".", fileformat), p, device = fileformat, width = w, height = h)
}


#produce the synteny figures
draw.linear(directory = "plots/", 
            output = paste("synt_", taxa1, taxa2, Sys.Date(), sep = ""), 
            paste("intermediate/chr_length_", taxa1, taxa2, ".txt", sep = ""), 
            paste("intermediate/chain_", taxa1, taxa2, ".txt", sep = ""),  
            fileformat = "png", w=13, h=5)

# ggsave(filename = paste("plots/synt_", taxa1, taxa2, Sys.Date(), ".png", sep = ""), 
#        device = "png", width = 13, height = 5)
# 

