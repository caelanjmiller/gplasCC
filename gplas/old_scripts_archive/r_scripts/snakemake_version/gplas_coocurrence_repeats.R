#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(cooccur))
suppressMessages(library(ggrepel))
suppressMessages(library(Biostrings))
suppressMessages(library(seqinr))

# Inputs
path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["clean_links"]]
path_prediction <- snakemake@input[["clean_prediction"]]
path_graph_contigs <- snakemake@input[["graph_contigs"]]
path_graph_repeats <- snakemake@input[["graph_repeats"]]
path_init_nodes <- snakemake@input[["repeat_nodes"]]
path_cov_variation <- snakemake@input[["coverage"]]
input_solutions <- snakemake@input[["solutions_repeat"]] #NEW
classifier <- snakemake@params[["classifier"]]
threshold <- snakemake@params[["threshold"]]
iterations <- snakemake@params[["iterations"]]
modularity_threshold <- snakemake@params[["modularity_threshold"]]
path_bins <- snakemake@input[["bins"]] #NEW
clean_repeats_path<- snakemake@input[["clean_repeats"]] #NEW
bold_sd_coverage <- snakemake@params[["bold_sd_coverage"]]

#import links between contigs
links <- read.table(file = path_links, header = TRUE)
#import contigs information from the graph
graph_contigs <- read.table(file = path_graph_contigs, header = TRUE)
#classify as small contigs shorter than 500 bp.
small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

repeats_graph <- read.table(file = path_graph_repeats, header=TRUE)
repeats <- repeats_graph
repeats$number <- gsub(pattern = '\\+',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)
repeats$number <- gsub(pattern = '\\-',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)

initialize_nodes <- read.table(file = path_init_nodes, header = TRUE)
initialize_nodes <- initialize_nodes[,1]

clean_pred <- read.table(file = path_prediction, header = TRUE)

#import the coverage variance of the chromosomal contigs.
max_variation <- read.table(file = path_cov_variation, header = TRUE)
max_variation <- as.numeric(max_variation[1,1])*as.numeric(bold_sd_coverage)
#allow more variance for small contigs.
max_variation_small <- max_variation*5.0

#import plasmid walks
solutions <- read.table(input_solutions ,sep=" ",fill=TRUE, col.names=c('walks','initial_classification','unitig_classification','path_coverage'))
#separate the nodes in the walk
max_nodes <- max(count.fields(input_solutions, sep = ","))
suppressWarnings(solutions<-separate(solutions, 'walks', into=as.character(1:max_nodes), sep=','))
#get the last node of solutions
solutions$last_node <- solutions[cbind(1:nrow(solutions), max.col(!is.na(solutions[,c(1:max_nodes)]), ties.method = 'last'))]
solutions$last_node_signless<-gsub(solutions$last_node, pattern='\\+',replacement = '')
solutions$last_node_signless<-gsub(solutions$last_node_signless, pattern='\\-',replacement = '')

#====Merge with bin data=============
#Import information from the bins and match it with total pairs
bins_data <- read.table(file = path_bins, header = TRUE)
#merge it
solutions <- merge(solutions, bins_data[,c(1,8)], by.x="last_node_signless",by.y="number", all=FALSE, all.x=TRUE)
#assign 0 to the chromsoome and -1 to repeats
solutions$Bin<-as.character(solutions$Bin)
solutions$Bin[is.na(solutions$Bin)]<-0
#remove connections to unbinned unitigs
solutions<-solutions[solutions$Bin!='Unbinned',]
#remove connections to repeats only
solutions<-solutions[solutions$unitig_classification!='Repeat',]

#get the initial node with out a symbol
solutions$initial_node<-gsub(solutions$`1`, pattern='[+]',replacement = '')
solutions$initial_node<-gsub(solutions$initial_node, pattern='[-]',replacement = '')
#combine inital nodes with bin in a single variable
solutions$repeat_bin<-paste(as.character(solutions$initial_node),as.character(solutions$Bin), sep='-')

#keep cases in which the last node and second node are the same (repeat directly connected to a unitig)
solutions$keep<-ifelse(solutions$`2`==solutions$last_node & solutions$unitig_classification!='Repeat','yes','no')
final_solutions<-solutions[solutions$keep=='yes',]

#----Analyze the remainig cases to check if we have the same bin upstream and downstream from the repeat--
#separate into positive and negative solutions
positive_solutions<-solutions[grepl("[+]",solutions$`1`),]
negative_solutions<-solutions[grepl("[-]",solutions$`1`),]
#Keep only the walks if they lead to same bin, or to the chromosome in both directions
positive_solutions_keep<-positive_solutions[positive_solutions$repeat_bin %in% negative_solutions$repeat_bin,  ]
positive_solutions_keep<-positive_solutions_keep[positive_solutions_keep$keep=='no', ]
negative_solutions_keep<-negative_solutions[negative_solutions$repeat_bin %in% positive_solutions$repeat_bin,  ]
negative_solutions_keep<-negative_solutions_keep[negative_solutions_keep$keep=='no', ]

## Filter only the valid walks
final_walks<-final_solutions[,c(1:max_nodes+1)]
final_positive_walks<-positive_solutions_keep[,c(1:max_nodes+1)]
final_negative_walks<-negative_solutions_keep[,c(1:max_nodes+1)]

solutions<-rbind(final_walks,final_positive_walks,final_negative_walks)

#replace NA values with empty values
solutions[is.na(solutions)]<-''

#create a null vector to contain all the nodes included in the plasmid walks
all_nodes <- NULL

#create a list with all the nodes that appear in the plasmid walks.
  for(solution in 1:nrow(solutions))
  {
    iteration <- solutions[solution,]
    iteration <- t(iteration)
    nodes <- as.character(iteration[,1])
    nodes <- nodes[nodes != '']
    all_nodes <- append(x = all_nodes, values = nodes, after = length(all_nodes))
  }

#filter out nodes that appear multiple times in the nodes
unique_nodes <- unique(all_nodes)
unique_nodes <- unique_nodes[unique_nodes != ""]
unique_nodes <- as.character(na.omit(unique_nodes))

#CREATE A CO-OCURRENCE MATRIX 
##Column names are the nodes included in plasmid walks.
##Each row is a new walk
##Assign a 1 if the node is present in walk and a 0 if node is not present
co_ocurrence <- data.frame(matrix(0, nrow = nrow(solutions), ncol = length(unique_nodes)))
colnames(co_ocurrence) <- unique_nodes
for(row in 1:nrow(co_ocurrence))
{
  iteration <- co_ocurrence[row,]
  particular_solution <- solutions[row,]
  particular_solution <- t(particular_solution)
  particular_solution <- particular_solution[,1]
  particular_solution <- particular_solution[particular_solution != ""]
  presence_absence <- ifelse(colnames(co_ocurrence) %in% particular_solution == TRUE, 1, 0)
  co_ocurrence[row,] <- presence_absence
}

suppressMessages(co_ocurrence <- apply(co_ocurrence, 2, as.integer))

#define the seed nodes
starting_nodes <- subset(unique_nodes, unique_nodes %in% unique(solutions[,1]))

#Number of iterations
number_iterations <- as.numeric(iterations)

#create a blank dataframe for co-ocurrence frequency (in network format)
#Start_node, Connecting_node,nr_occurences
total_pairs <- NULL

scalar1 <- function(x) {x / sqrt(sum(x^2))}

search_solutions <- as.data.frame(co_ocurrence)

#Get the number of times that two nodes co-ocuur in every walk
for(node in starting_nodes)
{
  index_col_node <-  which(colnames(search_solutions) == node)
  presence_node <- subset(search_solutions, search_solutions[index_col_node] == 1)
  index_presence_node <-  which(colnames(presence_node) == node)
  
  presence_node[index_presence_node] <- NULL
  sumatory <- colSums(presence_node)
  df_presence_absence <- data.frame(Starting_node = node, 
                                    Connecting_node = colnames(presence_node),
                                    weight = sumatory)
  
  total_pairs <- rbind(total_pairs, df_presence_absence)
}


total_pairs$Starting_node <- gsub(pattern = '[+]', replacement = '', x = total_pairs$Starting_node)
total_pairs$Starting_node <- gsub(pattern = '[-]', replacement = '', x = total_pairs$Starting_node)

total_pairs$Connecting_node <- gsub(pattern = '[+]', replacement = '', x = total_pairs$Connecting_node)
total_pairs$Connecting_node <- gsub(pattern = '[-]', replacement = '', x = total_pairs$Connecting_node)

#Filter-out cases of no-coocurrence.
total_pairs <- subset(total_pairs,total_pairs$weight > 1)

complete_node_info <- NULL

#Scale weigth
for(node in unique(total_pairs$Starting_node))
{
  first_node <- subset(total_pairs, total_pairs$Starting_node %in% node)
  particular_node <- NULL
  for(connecting_node in unique(first_node$Connecting_node))
  {
    first_second_nodes <- subset(first_node, first_node$Connecting_node %in% connecting_node)
    total_weight <- sum(first_second_nodes$weight)
    info_node <- data.frame(Starting_node = unique(first_second_nodes$Starting_node),
                            Connecting_node = unique(first_second_nodes$Connecting_node),
                            weight = total_weight)
    particular_node <- rbind(particular_node, info_node)
  }
  particular_node$scaled_weight <- scalar1(particular_node$weight)
  complete_node_info <- rbind(complete_node_info, particular_node)
}

total_pairs <- complete_node_info
initial_nodes <- gsub(pattern = '[+]', replacement = '', x = starting_nodes)
initial_nodes <- gsub(pattern = '[-]', replacement = '', x = initial_nodes)

total_pairs$Starting_node <- as.character(total_pairs$Starting_node)
total_pairs$Connecting_node <- as.character(total_pairs$Connecting_node)

#get a list of clean untigs
clean_unitigs<-clean_pred$number

#Filter out connected repeated elements. Keep only connections from starting nodes (repeats) unitigs.
total_pairs <- subset(total_pairs, total_pairs$Starting_node %in% initial_nodes)
total_pairs <- subset(total_pairs, total_pairs$Connecting_node %in% clean_unitigs)

#Import information from the bins and match it with total pairs
bins_data <- read.table(file = path_bins, header = TRUE)
bins_data$number<-as.character(bins_data$number)
bins_coverage<-bins_data %>% group_by(Bin) %>% summarise(bin_coverage=round(mean(coverage),2))
total_pairs <- merge(total_pairs, bins_data[,c(1,8)], by.x="Connecting_node",by.y="number", all=FALSE, all.x=TRUE)
#ASSIGN 0 to the chromosome
total_pairs$Bin<-as.character(total_pairs$Bin)
total_pairs$Bin[is.na(total_pairs$Bin)]<-0
total_pairs<-total_pairs[,c(2,5,3,4)]
colnames(total_pairs)[2]<-'Connecting_node'

#===Obtain the totality of co-courences 

#First check if we actually have co-ocurrence of unitigs.
if (length(total_pairs) > 0 && nrow(total_pairs) != 0) {
for(row in 1:nrow(total_pairs))
{
  total_pairs$Pair<-paste(total_pairs$Connecting_node,total_pairs$Starting_node, sep = '-')
}
} else {
#if gplas did not find any unitig-repeat connection
  print ("gplas couldn't find any walks connecting repeats to plasmid-nodes.")
  quit(status=1)
}
  
single_edge_counting <- total_pairs %>%
  group_by(Pair) %>%
  summarize(Weight = sum(weight))

weight_graph <- data.frame(From_to = as.character(str_split_fixed(string = single_edge_counting$Pair, pattern = '-', n = 2)[,2]),
                           To_from = as.character(str_split_fixed(string = single_edge_counting$Pair, pattern = '-', n = 2)[,1]),
                           weight = single_edge_counting$Weight)

total_scaled_weight <- NULL

full_graph_info <- NULL

#Get the data from coverages. Repeats and bins
clean_repeats <- read.table(file = clean_repeats_path, header = TRUE)
weight_graph <- merge(weight_graph, clean_repeats[,c(6,8)], by.x="From_to",by.y="number", all=FALSE, all.x=TRUE)
weight_graph <- merge(weight_graph, bins_coverage, by.x="To_from",by.y="Bin", all=FALSE, all.x=TRUE)
#assign a coverage of 1 to chromosomal unitigs
weight_graph$bin_coverage[is.na(weight_graph$bin_coverage)]<-1

#Explore if the combination of bins proposed by the algorithm is plausible based on coverage
max_variation <- read.table(file = path_cov_variation, header = TRUE)
max_variation <- as.numeric(max_variation[1,1])*as.numeric(bold_sd_coverage)
#get an empty dataframe
repeat_assignments<-NULL

#loop thru each of the repeats
for(node in unique(weight_graph$From_to))
{
  df_node <- subset(weight_graph, weight_graph$From_to == node) #subset data to include only one repeat at a time
  df_node$weight<-as.numeric(df_node$weight)
  df_node$rank<-base::rank(-df_node$weight, ties.method = 'random') #create a rank of the most likely connections, based on the co-ocurence count
  rank<-1 #start from the first ranking
  accumulated_cov<-0 #set up a variable for the accumulated coverage (this will be used to substract coverage from the repeat)
  while (rank<=max(df_node$rank)) { #loop thru the different ranks
    repeat_bin<-df_node[df_node$rank==rank,] #isolate the data from the highest rank bin
    
    if (as.numeric(repeat_bin$coverage+max_variation-accumulated_cov) >= repeat_bin$bin_coverage) { #if the repeat coverage is higher than the bin coverage then 
      accumulated_cov<-accumulated_cov+repeat_bin$bin_coverage #increse the accumulated coverage
      repeat_assignments<-rbind (repeat_assignments,repeat_bin) #add the information to the latest dataframe.
      rank<-rank+1 #move to the next rank
    } else {
      rank<-rank+1
    }
  }
}
repeat_assignments<-repeat_assignments[,c(1,2)]
names(repeat_assignments)<-c('Bin','number')

#separate results into plasmid and chromosome repeats
plasmid_repeats<-repeat_assignments[as.numeric(as.character(repeat_assignments$Bin))>=1,]
chromosome_repeats<-repeat_assignments[as.numeric(as.character(repeat_assignments$Bin))==0,]

if (nrow(plasmid_repeats)==0) {
  print("No repeats associated with plasmids were found")
  bins_data$Prob_Chromosome <- round(bins_data$Prob_Chromosome,2)
  bins_data$Prob_Plasmid <- round(bins_data$Prob_Plasmid,2)
  bins_data$coverage <- round(bins_data$coverage,2)
  bins_data$Bin<-as.character(bins_data$Bin)
  full_info_assigned<-bins_data
} else {
  print("We found repeated elements associated to plasmid predictions")
 #Get all the repeat nodes 
 pl_nodes <- clean_repeats[clean_repeats$number %in% plasmid_repeats$number,] # Selecting only contigs predicted as plasmid-derived 
 raw_number <- str_split_fixed(string = pl_nodes$Contig_name, pattern = '_', n = 2)[,1]
 pl_nodes$number <- gsub(pattern = 'S', replacement = '', x = raw_number)

 #Get all the information from the plasmid nodes (Lenght, coverage, bin number, etc)
 full_info_assigned <- merge(pl_nodes, plasmid_repeats, by = 'number')

 #Add information to the results file
 full_info_assigned$Contig_length <- NULL
 full_info_assigned$Prob_Chromosome <- round(full_info_assigned$Prob_Chromosome,2)
 full_info_assigned$Prob_Plasmid <- round(full_info_assigned$Prob_Plasmid,2)
 full_info_assigned$coverage <- round(full_info_assigned$coverage,2)
 full_info_assigned$Bin<-as.numeric(as.character(full_info_assigned$Bin))
 full_info_assigned$Prediction<-'Repeat'

 bins_data$Prob_Chromosome <- round(bins_data$Prob_Chromosome,2)
 bins_data$Prob_Plasmid <- round(bins_data$Prob_Plasmid,2)
 bins_data$coverage <- round(bins_data$coverage,2)
 bins_data$Bin<-as.character(bins_data$Bin)
 full_info_assigned<-rbind(full_info_assigned,bins_data)
}

#Create the fasta files
assembly_nodes <- readDNAStringSet(filepath = path_nodes)
df_nodes <- data.frame(Contig_name = names(assembly_nodes), Sequence = paste(assembly_nodes))
df_nodes <- merge(df_nodes, full_info_assigned, by = 'Contig_name')

#Write fasta files
for(component in unique(df_nodes$Bin))
{
  nodes_component <- subset(df_nodes, df_nodes$Bin == component)
  component_complete_name <- paste(snakemake@params[["sample"]], 'bin', sep = '_')
  filename <- paste('results/', component_complete_name, sep = '')
  filename <- paste(filename,component, sep = '_')
  filename <- paste(filename,'.fasta',sep = '')
  suppressWarnings(write.fasta(sequences = as.list(nodes_component$Sequence), names = nodes_component$Contig_name, file.out = filename))
}

results_summary<-df_nodes[,c(3,9)]
colnames(results_summary)[2] <- 'Bin'

suppressWarnings(write.table(x = full_info_assigned, file = snakemake@output[["results"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
suppressWarnings(write.table(x = results_summary, file = snakemake@output[["components"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))

#format chromosome repeats and print
chromosome_repeats<-chromosome_repeats[,c(2,1)]
chromosome_repeats$Bin<-'Chromosome'
suppressWarnings(write.table(x = chromosome_repeats, file = snakemake@output[["chromosome_repeats"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))





