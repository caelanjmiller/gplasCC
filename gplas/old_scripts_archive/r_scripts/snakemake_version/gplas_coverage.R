#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output

suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(Biostrings))

path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["links"]]
classifier <- snakemake@params[["classifier"]]
threshold <- snakemake@params[["threshold"]]
path_prediction <- snakemake@input[["prediction"]]


raw_nodes <- readDNAStringSet(filepath = path_nodes, format="fasta")
raw_contig_names <- names(raw_nodes)

# Checking if the assembly was generated using Unicycler or SPAdes

kc_check <- grep(pattern = 'KC', x = raw_contig_names)

if(length(kc_check) == length(raw_contig_names))
{
  length <- as.numeric(nchar(as.character(paste(raw_nodes))))
  kc_count <- str_split_fixed(string = raw_contig_names, pattern = ':', n = 4)[,3]
  kc_count <- as.numeric(gsub(pattern = '_', replacement = '', x = kc_count))
  kc_coverage <- kc_count/length
  coverage <- kc_coverage/median(kc_coverage)
}


raw_number <- str_split_fixed(string = raw_contig_names, pattern = '_', n = 2)[,1]
number <- gsub(pattern = 'S', replacement = '', x = raw_number)

# Checking if the assembly was generated using Unicycler

if(length(kc_check) != length(raw_contig_names))
{
  raw_length <- str_split_fixed(string = raw_contig_names, pattern = ':', n = 4)[,3]
  length <- gsub(pattern = '_dp', replacement = '', x = raw_length)
  coverage <- str_split_fixed(string = raw_contig_names, pattern = ':', n = 5)[,5]

}

contig_info <- data.frame(number = number,
                          length = length,
                          coverage = coverage,
                          Contig_name = raw_contig_names)


contig_info$length <- as.numeric(as.character(contig_info$length)) # Converting the column length into a numeric column
contig_info$coverage <- as.numeric(as.character(contig_info$coverage)) # Converting the column coverage into a coverage column

graph_pos_contigs <- contig_info
graph_pos_contigs$number <- paste(graph_pos_contigs$number, '+', sep = '')

graph_neg_contigs <- contig_info
graph_neg_contigs$number <- paste(graph_neg_contigs$number, '-', sep = '')

graph_contigs <- rbind(graph_pos_contigs, graph_neg_contigs)

write.table(x = graph_contigs, file = snakemake@output[["graph_contigs"]], row.names = FALSE)

small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

raw_links <- read.table(path_links)
colnames(raw_links) <- c('L','first_node','first_sign','second_node','second_sign','OM')

links <- NULL
for(row in 1:nrow(raw_links))
{
  connection <- raw_links[row,]

  reverse_first_sign <- ifelse(connection$first_sign == '+','-','+')
  reverse_second_sign <- ifelse(connection$second_sign == '+','-','+')

  info_reverse <- c('L',connection$second_node,reverse_second_sign,connection$first_node,reverse_first_sign,'OM')

  clean_first <- paste(connection$first_node, connection$first_sign, sep = '')
  clean_second <- paste(connection$second_node, connection$second_sign, sep = '')

  info_forward <- c(clean_first,'to',clean_second)

  clean_rev_first <- paste(connection$second_node, reverse_second_sign, sep = '')
  clean_rev_second <- paste(connection$first_node, reverse_first_sign, sep = '')

  info_reverse <- c(clean_rev_first,'to',clean_rev_second)

  full_connection <- rbind(info_forward, info_reverse)
  links <- rbind(links, full_connection)

}

links <- as.data.frame(links)

write.table(x = links, file = snakemake@output[["clean_links"]], row.names = FALSE)

unique_nodes <- unique(links$V1)

# Loop to create a dataframe with the number of links from each node

outdegree_info <- NULL
for(node in unique_nodes)
{
  repeat_links <- subset(links, links$V1 == node)
  repeat_links <- repeat_links[! duplicated(repeat_links),]

  node_info <- data.frame(number = node,
                          outdegree = length(repeat_links$V3))

  outdegree_info <- rbind(outdegree_info, node_info)
}


indegree_info <- NULL
for(node in unique_nodes)
{
  repeat_links <- subset(links, links$V3 == node)
  repeat_links <- repeat_links[! duplicated(repeat_links),]

  node_info <- data.frame(number = node,
                          indegree = length(repeat_links$V1))

  indegree_info <- rbind(indegree_info, node_info)
}

repeat_info <- merge(indegree_info, outdegree_info , by = 'number')

repeats <- subset(repeat_info, repeat_info$indegree > 1 | repeat_info$outdegree > 1) # Transposases may be identified as hubs in the graph. They should have more than one link

repeats_graph <- repeats

write.table(x = repeats_graph, file = snakemake@output[["graph_repeats"]])

repeats$number <- gsub(pattern = '\\+',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)
repeats$number <- gsub(pattern = '\\-',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)

pred <- read.table(path_prediction, header = TRUE)

clean_pred <- NULL

if(classifier == 'mlplasmids')
{
  clean_pred <- read.table(file = path_prediction, header = TRUE)
}

if(classifier == 'predict')
{
  clean_pred <- read.table(file = path_prediction, header = TRUE)
}

if(classifier == 'plasflow')
{
  plasflow_prediction <- read.table(file = path_prediction, header = TRUE)
  plasflow_prediction <- subset(plasflow_prediction, plasflow_prediction$contig_length > 1e3)
  labels <- as.character(plasflow_prediction$label)
  for(contig in 1:nrow(plasflow_prediction))
  {
    particular_contig <- plasflow_prediction[contig,]

    val <- max(as.numeric(particular_contig[,c(6:33)]))
    column <- which(particular_contig[,c(1:33)] == val)
    particular_contig$label <- colnames(particular_contig)[column]
    particular_contig$prob <- as.numeric(as.character(particular_contig[,column]))
    particular_contig$Prediction <- str_split_fixed(string = particular_contig$label, pattern = '\\.', n = 2)[,1]

    if(particular_contig$Prediction == 'chromosome')
    {
      particular_contig$Prediction <- 'Chromosome'
      particular_contig$Prob_Chromosome <- particular_contig$prob
      particular_contig$Prob_Plasmid <- 1-particular_contig$Prob_Chromosome

      clean_contig <- data.frame(Prob_Chromosome = particular_contig$Prob_Chromosome,
                                 Prob_Plasmid = particular_contig$Prob_Plasmid,
                                 Prediction = particular_contig$Prediction,
                                 Contig_name = particular_contig$contig_name,
                                 Contig_length = particular_contig$contig_length)

      clean_pred <- rbind(clean_pred, clean_contig)
    }

    if(particular_contig$Prediction == 'plasmid')
    {
      particular_contig$Prediction <- 'Plasmid'
      particular_contig$Prob_Plasmid <- particular_contig$prob
      particular_contig$Prob_Chromosome <- 1-particular_contig$Prob_Plasmid

      clean_contig <- data.frame(Prob_Chromosome = particular_contig$Prob_Chromosome,
                                 Prob_Plasmid = particular_contig$Prob_Plasmid,
                                 Prediction = particular_contig$Prediction,
                                 Contig_name = particular_contig$contig_name,
                                 Contig_length = particular_contig$contig_length)

      clean_pred <- rbind(clean_pred, clean_contig)
    }

  }
}


raw_number <- str_split_fixed(string = clean_pred$Contig_name, pattern = '_', n = 2)[,1]
clean_pred$number <- gsub(pattern = 'S', replacement = '', x = raw_number)

clean_pred$length <- NULL
clean_pred$coverage <- NULL
clean_pred$number <- NULL

clean_pred <- merge(clean_pred, contig_info, by = 'Contig_name')

#get the plasmid nodes without repeats
final_prediction <- clean_pred[! clean_pred$number %in% repeats$number,]
write.table(x = final_prediction, file = snakemake@output[["clean_prediction"]], row.names = FALSE)

#Get the isolated nodes in the graph. -This nodes should appear in the contig list, but not in the graph
unique_nodes_singless<-as.data.frame(links$V1)
names(unique_nodes_singless)<-'node_number'
unique_nodes_singless$node_number<-gsub(unique_nodes_singless$node_number,pattern="[+]",replacement = '')
unique_nodes_singless$node_number<-gsub(unique_nodes_singless$node_number,pattern="[-]",replacement = '')
isolated_nodes<-contig_info[! contig_info$number %in% unique_nodes_singless$node_number ,]
isolated_nodes <- merge(clean_pred, isolated_nodes[4], by = 'Contig_name')
isolated_nodes <- subset(isolated_nodes, isolated_nodes$Prob_Plasmid >= as.numeric(as.character(threshold))) # Selecting only contigs predicted as plasmid-derived
write.table(x = isolated_nodes, file = snakemake@output[["isolated_nodes"]], row.names = FALSE)

#Get the repeat information.
repeats_final <- clean_pred[clean_pred$number %in% repeats$number,]
write.table(x = repeats_final, file = snakemake@output[["clean_repeats"]], row.names = FALSE)

###########################################################################################################################################

# Selecting the plasmid seeds in our graph
pl_nodes <- subset(final_prediction, final_prediction$Prob_Plasmid >= as.numeric(as.character(threshold))) # Selecting only contigs predicted as plasmid-derived
pl_nodes <- subset(pl_nodes, pl_nodes$Contig_length > 5e2) # Filtering the plasmid seeds by length 
pl_nodes <- pl_nodes[! pl_nodes$number %in% repeats$number,] # From these contigs we remove contigs that could correspond to transposases
pl_nodes <- pl_nodes[order(pl_nodes$length, decreasing = TRUE),] # Sorting the contigs based on length

initialize_nodes <- unique(pl_nodes$number) # These contigs are going to be used as seeds to start finding the solutions

write.table(x = initialize_nodes, file = snakemake@output[["initialize_nodes"]], row.names = FALSE)

###########################################################################################################################################

# Select the repeat nodes that will be used to explore their surroundings
repeats_nodes<-unique(repeats_final$number)
write.table(x = repeats_nodes, file = snakemake@output[["repeat_nodes"]], row.names = FALSE)

###########################################################################################################################################

# First we need to remove the contigs with more than one link as well (possible transposases to get a better estimation about what is the coverage variation present in the graph)

chr_contigs <- subset(clean_pred, clean_pred$Prob_Chromosome > 0.7 & clean_pred$Contig_length > 1e3) # The length is tricky since shorter contigs have higher k-mer coverage fluctuations

cov_estimation <- chr_contigs[! chr_contigs$number %in% repeats$number, ]

# Different measures could be consider, we can use the median absolute deviation (mad) which can be more robust in case our population does not follow a normal distribution

# This is the coverage spread that we find in our data

cov <- mad(cov_estimation$coverage) # We assign the mad to the maximum variation k-mer variation that we allow in our data

sd_estimation <- sd(cov_estimation$coverage)

max_variation <- cov*1.0
max_variation_small <- cov*5.0

write.table(x = sd_estimation, file = snakemake@output[["coverage"]], row.names = FALSE)
