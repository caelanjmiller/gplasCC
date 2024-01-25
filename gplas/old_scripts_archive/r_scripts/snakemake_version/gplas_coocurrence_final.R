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
path_init_nodes <- snakemake@input[["initialize_nodes"]]
path_isolated_nodes <- snakemake@input[["isolated_nodes"]]
path_cov_variation <- snakemake@input[["coverage"]]
input_solutions <- snakemake@input[["solutions"]]
classifier <- snakemake@params[["classifier"]]
threshold <- snakemake@params[["threshold"]]
iterations <- snakemake@params[["iterations"]]
modularity_threshold <- snakemake@params[["modularity_threshold"]]




links <- read.table(file = path_links, header = TRUE)
graph_contigs <- read.table(file = path_graph_contigs, header = TRUE)

small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

repeats_graph <- read.table(file = path_graph_repeats)
repeats <- repeats_graph
repeats$number <- gsub(pattern = '\\+',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)
repeats$number <- gsub(pattern = '\\-',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)

initialize_nodes <- read.table(file = path_init_nodes, header = TRUE)
initialize_nodes <- initialize_nodes[,1]

clean_pred <- read.table(file = path_prediction, header = TRUE)


max_variation <- read.table(file = path_cov_variation, header = TRUE)
max_variation <- as.numeric(max_variation[1,1])
max_variation_small <- max_variation*5.0

max_nodes <- max(count.fields(input_solutions, sep = ","))

solutions <- read.table(input_solutions ,sep=",",fill=TRUE,col.names=1:max_nodes)

all_nodes <- NULL
  for(solution in 1:nrow(solutions))
  {
    iteration <- solutions[solution,]
    iteration <- t(iteration)
    nodes <- as.character(iteration[,1])
    nodes <- nodes[nodes != '']
    all_nodes <- append(x = all_nodes, values = nodes, after = length(all_nodes))
  }

unique_nodes <- unique(all_nodes)
unique_nodes <- unique_nodes[unique_nodes != ""]
unique_nodes <- as.character(na.omit(unique_nodes))

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

starting_nodes <- subset(unique_nodes, unique_nodes %in% unique(solutions[,1]))

number_iterations <- as.numeric(iterations)

total_pairs <- NULL

scalar1 <- function(x) {x / sqrt(sum(x^2))}

search_solutions <- as.data.frame(co_ocurrence)

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


## Looking back in the solutions if there are circular graphs

circular_sequences <- NULL
for(solution in 1:nrow(solutions))
{
  iteration <- solutions[solution,]
  iteration <- t(iteration)
  nodes <- as.character(iteration[,1])
  nodes <- nodes[nodes != '']
  if(length(nodes) > 1 & nodes[1] == nodes[length(nodes)])
  {
    
    circular_sequences <- rbind(circular_sequences, c(nodes[1],nodes[length(nodes)]))
  }
}

if(is.null(circular_sequences) == FALSE)
{
  
no_duplicated <- circular_sequences[!duplicated(circular_sequences),]
for(combination in 1:nrow(no_duplicated))
{
  combi <- no_duplicated[combination,]
  total_ocurrences <- subset(circular_sequences, circular_sequences[,2] == combi[2])
  if(nrow(total_ocurrences) == number_iterations)
  {
    df_test <- data.frame(Starting_node = combi[1],
                          Connecting_node = combi[2],
                          weight = nrow(total_ocurrences))
    
    total_pairs <- rbind(total_pairs, df_test)
  }
  
}
}

total_pairs$Starting_node <- gsub(pattern = '\\+', replacement = '', x = total_pairs$Starting_node)
total_pairs$Starting_node <- gsub(pattern = '\\-', replacement = '', x = total_pairs$Starting_node)

total_pairs$Connecting_node <- gsub(pattern = '\\+', replacement = '', x = total_pairs$Connecting_node)
total_pairs$Connecting_node <- gsub(pattern = '\\-', replacement = '', x = total_pairs$Connecting_node)

total_pairs <- subset(total_pairs,total_pairs$weight > 1)

complete_node_info <- NULL

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

initial_nodes <- gsub(pattern = '\\+', replacement = '', x = starting_nodes)
initial_nodes <- gsub(pattern = '\\-', replacement = '', x = initial_nodes)

total_pairs$Starting_node <- as.character(total_pairs$Starting_node)
total_pairs$Connecting_node <- as.character(total_pairs$Connecting_node)

total_pairs <- subset(total_pairs, total_pairs$Starting_node %in% initial_nodes)
total_pairs <- subset(total_pairs, total_pairs$Connecting_node %in% initial_nodes)


weight_counting <- NULL

if (nrow(total_pairs) == 0) {
  print ("gplas couldn't find any walks connecting plasmid-predicted nodes. Plasmid nodes will be classified as Unbinned. If this is unexpected, please assemble your genome with different parameters or with a different tool and re-run gplas.")
  pl_unassigned <- subset(clean_pred, clean_pred$Prob_Plasmid > as.numeric(as.character(threshold))) # Selecting only contigs predicted as plasmid-derived 
  
  raw_number <- str_split_fixed(string = pl_unassigned$Contig_name, pattern = '_', n = 2)[,1]
  pl_unassigned$number <- gsub(pattern = 'S', replacement = '', x = raw_number)
  pl_unassigned$Component <- 'Unbinned'
  
  pl_unassigned$Contig_length <- NULL
  pl_unassigned$Prob_Chromosome <- round(pl_unassigned$Prob_Chromosome,2)
  pl_unassigned$Prob_Plasmid <- round(pl_unassigned$Prob_Plasmid,2)
  pl_unassigned$coverage <- round(pl_unassigned$coverage,2)
  
  assembly_nodes <- readDNAStringSet(filepath = path_nodes)
  df_nodes <- data.frame(Contig_name = names(assembly_nodes), Sequence = paste(assembly_nodes))
  
  df_nodes <- merge(df_nodes, pl_unassigned, by = 'Contig_name')
  
  for(component in unique(df_nodes$Component))
  {
    nodes_component <- subset(df_nodes, df_nodes$Component == component)
    component_complete_name <- paste(snakemake@params[["sample"]], 'bin', sep = '_')
    filename <- paste('results/', component_complete_name, sep = '')
    filename <- paste(filename,component, sep = '_')
    filename <- paste(filename,'.fasta',sep = '')
    suppressWarnings(write.fasta(sequences = as.list(nodes_component$Sequence), names = nodes_component$Contig_name, file.out = filename))
    
  }
  
  colnames(pl_unassigned)[8] <- 'Bin'
  results_subgraph<-pl_unassigned[,c(5,8)]
  #colnames(results_subgraph)[2] <- 'Bin'
  
  suppressWarnings(write.table(x = pl_unassigned, file = snakemake@output[["results"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
  suppressWarnings(write.table(x = results_subgraph, file = snakemake@output[["components"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
  
  png(filename=snakemake@output[["plot_graph"]],width = 700, height = 700)
  plot.new()
  dev.off()
  quit(status=0)
} else {
for(row in 1:nrow(total_pairs))
{
  initial_node <- as.numeric(total_pairs[row,1])
  connecting_node <- as.numeric(total_pairs[row,2])
  
  raw_count <- as.numeric(total_pairs[row,3])
  
  if(initial_node < connecting_node)
  {
    pair <- paste(initial_node,connecting_node, sep = '-')
    df_count <- data.frame(Pair = pair,
               Count = raw_count)
  }
  else
  {
    pair <- paste(connecting_node,initial_node, sep = '-')
    df_count <- data.frame(Pair = pair,
               Count = raw_count)
  }
  
  weight_counting <- rbind(weight_counting, df_count)
}
}  
  
single_edge_counting <- weight_counting %>%
  group_by(Pair) %>%
  summarize(Weight = sum(Count))



weight_graph <- data.frame(From_to = as.character(str_split_fixed(string = single_edge_counting$Pair, pattern = '-', n = 2)[,1]),
                           To_from = as.character(str_split_fixed(string = single_edge_counting$Pair, pattern = '-', n = 2)[,2]),
                           weight = single_edge_counting$Weight)

total_scaled_weight <- NULL

full_graph_info <- NULL

for(node in unique(weight_graph$From_to))
{
  df_node <<- subset(weight_graph, weight_graph$From_to == node)
  scaled_weight <- scalar1(df_node$weight)
  df_node$scaled_weight <- scaled_weight
  full_graph_info <- rbind(full_graph_info, df_node)
}

full_graph_info$weight <- full_graph_info$scaled_weight

weight_graph$width <- full_graph_info$scaled_weight*5

graph_pairs <- graph_from_data_frame(weight_graph, directed = FALSE)

### Important: Only for viz purposes 

viz_graph <- weight_graph
viz_graph$weight <- NULL
graph_viz <- graph_from_data_frame(viz_graph, directed = FALSE)
graph_viz <- igraph::simplify(graph_viz, remove.multiple=FALSE)

# Simplifying the graph 

no_loops_graph <- igraph::simplify(graph_pairs, remove.multiple=FALSE)


# Starting taking the decision whether to perform a partition or not 


components_graph <- decompose(no_loops_graph, mode = c("weak"), max.comps = NA,
                              min.vertices = 2)

partitioning_components <- function(graph)
{
  
  # Spinglass algorithm
  graph_spin <- cluster_spinglass(graph)
  
  
  # Walktrap algorithm 
  graph_walktrap <- cluster_walktrap(graph)
  
  # Leading eigen values 
  
  graph_eigen <- cluster_leading_eigen(graph)
  
  
  # Louvain method 
  
  graph_louvain <- cluster_louvain(graph)
  
  
  # Community detection based on propagating labels 
  
  graph_propag <- cluster_label_prop(graph)
  
  partition_info <<- data.frame(Algorithm = c('Walktrap','Leading-eigen','Louvain'),
                                Modularity = c(modularity(graph_walktrap), 
                                               modularity(graph_eigen), 
                                               modularity(graph_louvain)))
  
}

complete_partition_info <- NULL


info_comp_member <- components(no_loops_graph)$membership
original_components <- unique(info_comp_member)
info_comp_size <- components(no_loops_graph)$csize

node_and_component <- data.frame(Node = names(info_comp_member),
                                 Original_component = as.character(info_comp_member))


information_components <- data.frame(Original_component = original_components,
                                     Size = info_comp_size)

full_info_components <- merge(node_and_component, information_components, by = 'Original_component')


if(length(components_graph) >= 1)
{
  
for(component in 1:length(components_graph))
{
  subgraph <- components_graph[[component]]
  partitioning_components(subgraph)
  first_node <- names(V(subgraph))[1]
  
  
  info_first_node <- subset(full_info_components, full_info_components$Node == first_node)
  
  partition_info$Original_component <- info_first_node$Original_component
  complete_partition_info <- rbind(complete_partition_info, partition_info)
  
}

complete_partition_info$Modularity <- round(complete_partition_info$Modularity,2)

modularity_threshold <- as.numeric(modularity_threshold)

complete_partition_info$Decision <- ifelse(complete_partition_info$Modularity >= modularity_threshold, 'Split', 'No_split')
}

singletons_component <- which((components(no_loops_graph)$csize == 1) == TRUE)


if(length(singletons_component) > 0)
{
  df_singletons <- data.frame(Algorithm = 'Independent-single_component', 
                              Modularity = 0, 
                              Original_component = singletons_component,
                              Decision = 'No_split')
  
  complete_partition_info <- rbind(complete_partition_info, df_singletons)
}

complete_partition_info

contigs_membership <- NULL
internal_component <- 1

complete_partition_info <- complete_partition_info[order(complete_partition_info$Original_component),]

for(component in unique(complete_partition_info$Original_component))
{
  decision_comp <- subset(complete_partition_info, complete_partition_info$Original_component == component)
  
  
  if(decision_comp$Algorithm[1] != 'Independent-single_component')
  {
    split_decision <- nrow(subset(decision_comp, decision_comp$Decision == 'Split'))
    no_split_decision <- nrow(subset(decision_comp, decision_comp$Decision == 'No_split'))
    if(split_decision >= no_split_decision)
    {
      algorithm_to_split <- subset(decision_comp, decision_comp$Modularity == max(decision_comp$Modularity))
      algorithm <- as.character(algorithm_to_split$Algorithm[1])
      
      if(algorithm == 'Walktrap')
      {
        graph_walktrap <- walktrap.community(components_graph[[internal_component]])
        spl_membership <- graph_walktrap$membership
        spl_names <- graph_walktrap$names
      }
      if(algorithm == 'Leading-eigen')
      {
        graph_eigen <- cluster_leading_eigen(components_graph[[internal_component]])
        spl_membership <- graph_eigen$membership
        spl_names <- graph_eigen$names
      }
      if(algorithm == 'Louvain')
      {
        graph_louvain <- cluster_louvain(components_graph[[internal_component]])
        spl_membership <- graph_louvain$membership
        spl_names <- graph_louvain$names
      }
      if(algorithm == 'Propagating-labels')
      {
        graph_propag <- cluster_label_prop(components_graph[[internal_component]])
        spl_membership <- graph_propag$membership
        spl_names <- graph_propag$names
      }
      
      info_decision_component <- data.frame(Algorithm = algorithm,
                                            Original_component = component,
                                            Bin = spl_membership, 
                                            Contig = spl_names)
      
      
      contigs_membership <- rbind(contigs_membership, info_decision_component)
      
    }
    else
    {
      nodes_component <- subset(full_info_components, full_info_components$Original_component == component)
      info_decision_component <- data.frame(Algorithm = 'Not_split_component',
                                            Original_component = component,
                                            Bin = 1, 
                                            Contig = nodes_component$Node)
      contigs_membership <- rbind(contigs_membership, info_decision_component)
    }
    internal_component <- internal_component + 1
  }
  else
  {
    nodes_component <- subset(full_info_components, full_info_components$Original_component == component)
    
    info_decision_component <- data.frame(Algorithm = 'Independent-single_component',
                                          Original_component = component,
                                          Bin = 1, 
                                          Contig = nodes_component$Node)
    contigs_membership <- rbind(contigs_membership, info_decision_component)
  }
  
}

contigs_membership$Cluster <- paste(contigs_membership$Original_component, contigs_membership$Bin, sep = '-')

contigs_membership <- contigs_membership[order(contigs_membership$Cluster),]

contigs_membership$Final_cluster <- contigs_membership$Cluster

for(number in 1: length(unique(contigs_membership$Final_cluster)))
{
  cluster_name <- unique(contigs_membership$Cluster)[number]
  contigs_membership$Final_cluster <- gsub(pattern = cluster_name, replacement = number, x = contigs_membership$Final_cluster)
}


set_colors <- c("lightblue","#d49f36","#507f2d","#84b67c","#a06fda","#df462a","#5a51dc","#5b83db","#c76c2d","#4f49a3","#552095","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

contigs_membership$Color <- set_colors[as.numeric(contigs_membership$Final_cluster)]

order_contigs <- names(V(graph_viz))

contigs_membership <- contigs_membership[match(order_contigs, as.character(contigs_membership$Contig)),]

V(graph_viz)$color <- contigs_membership$Color

png(filename=snakemake@output[["plot_graph"]],width = 700, height = 700)
plot.igraph(graph_viz)
dev.off()


results_subgraph <- data.frame(number = contigs_membership$Contig,
                               Component = contigs_membership$Final_cluster)
 
pl_nodes <- subset(clean_pred, clean_pred$Prob_Plasmid >= as.numeric(as.character(threshold))) # Selecting only contigs predicted as plasmid-derived 

raw_number <- str_split_fixed(string = pl_nodes$Contig_name, pattern = '_', n = 2)[,1]
pl_nodes$number <- gsub(pattern = 'S', replacement = '', x = raw_number)


pl_notassigned <- subset(pl_nodes,! pl_nodes$number %in% results_subgraph$number)

pl_repeats <- subset(pl_notassigned, pl_notassigned$number %in% repeats$number)

pl_unassigned <- subset(pl_notassigned,! pl_notassigned$number %in% repeats$number)

#If there are isolated nodes, remove them from the unbinned and create a new category
isolated_nodes <- read.table(file = path_isolated_nodes, header = TRUE)
pl_isolated <- subset(pl_unassigned, pl_notassigned$number %in% isolated_nodes$number)
pl_unassigned <- subset(pl_unassigned,! pl_unassigned$number %in% isolated_nodes$number)


pl_assigned <- subset(pl_nodes, pl_nodes$number %in% results_subgraph$number)
full_info_assigned <- merge(pl_assigned, results_subgraph, by = 'number')

if(nrow(pl_unassigned) >= 1)
{
  pl_unassigned$Component <- 'Unbinned'
  full_info_assigned <- rbind(full_info_assigned, pl_unassigned)
}

#Add isolated nodes category
if(nrow(pl_isolated) >= 1)
{
  isolated_nr<-1
  while (isolated_nr<= nrow(pl_isolated)) {
    isolated_identification<-paste('Isolated',as.character(isolated_nr), sep='_') 
    pl_isolated$Component[isolated_nr] <- isolated_identification
    isolated_nr<-isolated_nr+1
  }
  full_info_assigned <- rbind(full_info_assigned, pl_isolated)
}

if(nrow(pl_repeats) >= 1)
{
  pl_repeats$Component <- 'Repeat_like'
  full_info_assigned <- rbind(full_info_assigned, pl_repeats)
}

full_info_assigned$Contig_length <- NULL
full_info_assigned$Prob_Chromosome <- round(full_info_assigned$Prob_Chromosome,2)
full_info_assigned$Prob_Plasmid <- round(full_info_assigned$Prob_Plasmid,2)
full_info_assigned$coverage <- round(full_info_assigned$coverage,2)

assembly_nodes <- readDNAStringSet(filepath = path_nodes)
df_nodes <- data.frame(Contig_name = names(assembly_nodes), Sequence = paste(assembly_nodes))

df_nodes <- merge(df_nodes, full_info_assigned, by = 'Contig_name')

for(component in unique(df_nodes$Component))
{
  nodes_component <- subset(df_nodes, df_nodes$Component == component)
  component_complete_name <- paste(snakemake@params[["sample"]], 'bin', sep = '_')
  filename <- paste('results/', component_complete_name, sep = '')
  filename <- paste(filename,component, sep = '_')
  filename <- paste(filename,'.fasta',sep = '')
  suppressWarnings(write.fasta(sequences = as.list(nodes_component$Sequence), names = nodes_component$Contig_name, file.out = filename))

}

colnames(full_info_assigned)[8] <- 'Bin'
colnames(results_subgraph)[2] <- 'Bin'

suppressWarnings(write.table(x = full_info_assigned, file = snakemake@output[["results"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))
suppressWarnings(write.table(x = results_subgraph, file = snakemake@output[["components"]], append = TRUE, row.names = FALSE, quote = FALSE, col.names = TRUE))





