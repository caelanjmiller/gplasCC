#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 

suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(cooccur))
suppressMessages(library(ggrepel))

# Inputs

path_nodes <- snakemake@input[["nodes"]]
path_reference <- snakemake@input[["reference"]]
path_links <- snakemake@input[["clean_links"]]
path_prediction <- snakemake@input[["clean_prediction"]]
path_graph_contigs <- snakemake@input[["graph_contigs"]]
path_graph_repeats <- snakemake@input[["graph_repeats"]]
path_init_nodes <- snakemake@input[["initialize_nodes"]]
path_cov_variation <- snakemake@input[["coverage"]]
input_solutions <- snakemake@input[["solutions"]]
path_alignments <- snakemake@input[["alignments"]]
path_components <- snakemake@input[["components"]]

classifier <- snakemake@params[["classifier"]]
number_iterations <- snakemake@params[["iterations"]]
species <- snakemake@params[["species"]]
sample <- snakemake@params[["name"]]


# # Checking the evaluation 
# 
# path_nodes <- '/home/sergi/gplas/gplas_input/SAMN10819804_pflow_raw_nodes.fasta'
# path_reference <- '/home/sergi/gplas/independent_set/SAMN10819804/assembly.fasta'
# path_links <- '/home/sergi/gplas/coverage/SAMN10819804_pflow_clean_links.tab'
# path_prediction <- '/home/sergi/gplas/plasflow_prediction/SAMN10819804_pflow_plasmid_prediction.tab'
# path_graph_contigs <- '/home/sergi/gplas/coverage/SAMN10819804_pflow_graph_contigs.tab'
# path_graph_repeats <- '/home/sergi/gplas/coverage/SAMN10819804_pflow_repeats_graph.tab'
# path_init_nodes <- '/home/sergi/gplas/coverage/SAMN10819804_pflow_initialize_nodes.tab'
# path_cov_variation <- '/home/sergi/gplas/coverage/SAMN10819804_pflow_estimation.txt'
# input_solutions <- '/home/sergi/gplas/walks/SAMN10819804_pflow_solutions.csv'
# path_alignments <- '/home/sergi/gplas/evaluation/SAMN10819804_pflow_alignment_test.txt'
# path_components <- '/home/sergi/gplas/results/SAMN10819804_pflow_components.tab'
# 
# classifier <- 'plasflow'
# number_iterations <- 20
# species <- 'unknown'
# sample <- 'SAMN10819804'
# 


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

raw_nodes <- readDNAStringSet(filepath = path_nodes, format="fasta")

raw_contig_names <- names(raw_nodes)


# Checking if the assembly was generated using Velvet or SPAdes

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
                          Contig_length = length,
                          length = length,
                          coverage = coverage,
                          Contig_name = raw_contig_names)


contig_info$length <- as.numeric(as.character(contig_info$length))
contig_info$Contig_length <- as.numeric(as.character(contig_info$Contig_length)) # Converting the column length into a numeric column
contig_info$coverage <- as.numeric(as.character(contig_info$coverage)) # Converting the column coverage into a coverage column


alignment_quast <- read.table(file = path_alignments, sep = ' ', header = TRUE)

raw_complete_genome <- readDNAStringSet(filepath = path_reference, format="fasta")
raw_names_complete_genome <- names(raw_complete_genome)
length_complete_genomes <- as.numeric(nchar(as.character(paste(raw_complete_genome))))

info_complete_genomes <- data.frame(Reference = raw_names_complete_genome,
                                    Reference_length = length_complete_genomes)

info_complete_genomes$Reference_number <- c(1:nrow(info_complete_genomes))

info_complete_genomes$Reference <- gsub(pattern = '\\|', replacement = '_', x = info_complete_genomes$Reference)
info_complete_genomes$Reference <- gsub(pattern = '=', replacement = '_', x = info_complete_genomes$Reference)

info_complete_genomes$Reference <- gsub(pattern = ' ', replacement = '_', x = info_complete_genomes$Reference)

alignment_quast <- merge(alignment_quast, info_complete_genomes, by = 'Reference')
alignment_quast$Contig_name <- alignment_quast$Contig
alignment_quast$Contig_number <- str_split_fixed(alignment_quast$Contig, "_", 6)[,1] # Extracting the contig number
alignment_quast$Contig_number <- gsub(pattern = 'S', replacement = '', x = alignment_quast$Contig_number)

contig_info$Contig_number <- contig_info$number

alignment_quast <- merge(alignment_quast, contig_info, by = 'Contig_number')
alignment_quast$Type <- ifelse(alignment_quast$Reference_length > 4e5, 'Chromosome','Plasmid')
alignment_quast$Type <- paste(alignment_quast$Type,alignment_quast$Reference_number, sep = '')

# Remove duplicates present in the alignment (e.g. based on the fact that plasmids are circular)

non_duplicates_quast <<- alignment_quast[! duplicated(alignment_quast),]

non_transposases_quast <<- non_duplicates_quast[! non_duplicates_quast$Contig_number %in% repeats$number,]

gold_standard <- non_transposases_quast

# Including information about the number of contigs mapping to the reference genome and the total number of base-pairs 

gold_standard <- gold_standard %>%
  group_by(Type) %>%
  mutate(sum(Contig_length),
         count = n())

colnames(gold_standard)[c(13,14)] <- c('Reference_sum_bp','Reference_sum_contigs')
evaluation_ref_genomes <- gold_standard

# Dataframe with all the stats regarding the reference genomes (Number of contigs and number of base-pairs)

stats_gold_standard <- gold_standard %>% 
  group_by(Type) %>%
  summarise(sum(Contig_length),
            count = n())# Adding the contig direction into the results (we use this information to merge our gold_standard with the prediction) 

colnames(stats_gold_standard) <- c('Type','Reference_sum_bp','Reference_sum_contigs')


# Creating a dataframe to merge the quast mapping against the prediction from gplas We need to include the directionality of the contigs 

contigs_reference_genome <- gold_standard %>%
  group_by(Type) %>%
  summarise(count = n(),
            sum(Contig_length))

colnames(contigs_reference_genome) <- c('Type','Reference_sum_contigs','Reference_sum_bp')


truth_set <- subset(gold_standard, select = c('Contig_number','Type'))
truth_set <- truth_set[order(truth_set$Type, decreasing = TRUE),]

truth_set <- truth_set[! duplicated(truth_set$Contig_number),]

count_contigs_truth <- truth_set %>%
  group_by(Type) %>%
  count()

count_contigs_truth$connections <- choose(count_contigs_truth$n, 2)

#################### Fifth part
# Merging the prediction with the ground truth 
# Merging the reference results with the prediction given

results_subgraph <- read.table(file = path_components)

colnames(results_subgraph) <- c('Contig_number','Component')

results <- merge(results_subgraph, gold_standard, by = 'Contig_number')

singletons <- initialize_nodes[! initialize_nodes %in% results_subgraph$Contig_number]

singletons_length <- subset(contig_info, contig_info$number %in% singletons)


benchmark <- results %>%
  group_by(Component, Type) %>%
  summarise(count = n(), sum_comp = sum(Contig_length))

benchmark$Component <- as.character(benchmark$Component)


singletons_info <- c('No_component',NA,length(singletons),sum(singletons_length$length))

benchmark <- rbind(as.data.frame(benchmark), singletons_info)

# We need to consider what is the reference genome predominating in a particular bin, e.g. completeness is based on this assumption whereas precision (or purity of the bin) is reference-independent

sort_benchmark <- benchmark[order(benchmark$Component),]
sort_benchmark <- sort_benchmark[order(as.numeric(sort_benchmark$count), decreasing = TRUE),]

total_component_info <- sort_benchmark

reference_component <- sort_benchmark[! duplicated(sort_benchmark$Component),]
reference_component <- reference_component[order(reference_component$Component, decreasing = FALSE),]

reference_component$connections <- choose(as.numeric(as.character((reference_component$count))), 2)

completeness_df <- NULL

for(component in unique(reference_component$Component))
{
  comp_info <- subset(reference_component, reference_component$Component == component)
  total_component <- subset(total_component_info, total_component_info$Component == component)

  total_component_contigs <- sum(as.numeric(total_component$count))
  total_component_bp <- sum(as.numeric(total_component$sum_comp))
  
  total_component_connections <- choose(total_component_contigs, 2)
  
  chr_component <- subset(total_component, total_component$Type == 'Chromosome1')
  
  if(nrow(chr_component) > 0)
  {
    chr_contigs <- chr_component$count
    chr_bp <- chr_component$sum_comp
  }
  else
  {
    chr_contigs <- 0
    chr_bp <- 0
  }
  if(component == 'No_component')
  {
    completeness_info <- data.frame(sample = sample,
                                    species = species,
                                    classifier = classifier, 
                                    factor = 1.0,
                                    number_iterations = number_iterations,
                                    variation = max_variation,
                                    component = component,
                                    total_contigs = comp_info$count,
                                    total_bp = comp_info$sum_comp,
                                    total_connections = 0,
                                    reference = 0,
                                    contigs_ref_component = 0,
                                    bp_ref_component = 0,
                                    connections_ref_component = 0,
                                    contigs_reference = 0,
                                    bp_reference = 0,
                                    connections_reference = 0,
                                    completeness_contigs = 0,
                                    completeness_bp = 0,
                                    completeness_connections = 0,
                                    chr_contigs = 0,
                                    chr_bp = 0)
  }
  
  if(component != 'No_component')
  {
    
    count_truth <- subset(contigs_reference_genome$Reference_sum_contigs, contigs_reference_genome$Type == comp_info$Type)
    
    count_truth_bp <- subset(contigs_reference_genome$Reference_sum_bp, contigs_reference_genome$Type == comp_info$Type)
    
    count_truth_connections <- choose(count_truth,2)
    
    completeness_contigs <- as.numeric(as.character(comp_info$count))/count_truth
    completeness_bp <- as.numeric(as.character(comp_info$sum_comp))/count_truth_bp
    
    completeness_connections <- comp_info$connections/count_truth_connections
    
    completeness_info <- data.frame(sample = sample,
                                    species = species,
                                    classifier = classifier, 
                                    factor = 1.0,
                                    number_iterations = number_iterations,
                                    variation = max_variation,
                                    component = component,
                                    total_contigs = as.numeric(as.character(total_component_contigs)),
                                    total_bp = as.numeric(as.character(total_component_bp)),
                                    total_connections = total_component_connections,
                                    reference = comp_info$Type,
                                    contigs_ref_component = comp_info$count,
                                    bp_ref_component = comp_info$sum_comp,
                                    connections_ref_component = comp_info$connections,
                                    contigs_reference = count_truth,
                                    bp_reference = count_truth_bp,
                                    connections_reference = count_truth_connections,
                                    completeness_contigs = completeness_contigs,
                                    completeness_bp = completeness_bp,
                                    completeness_connections = completeness_connections,
                                    chr_contigs = chr_contigs,
                                    chr_bp = chr_bp)
    
  }
  suppressWarnings(completeness_df <- rbind(completeness_df, completeness_info))
}

metrics_components <- completeness_df
metrics_components$contigs_ref_component <- as.numeric(as.character(metrics_components$contigs_ref_component))
metrics_components$contigs_reference <- as.numeric(as.character(metrics_components$contigs_reference))
metrics_components$total_contigs <- as.numeric(as.character(metrics_components$total_contigs))
metrics_components$Precision_contigs <- metrics_components$contigs_ref_component/metrics_components$total_contigs

metrics_components$bp_ref_component <- as.numeric(as.character(metrics_components$bp_ref_component))
metrics_components$bp_reference <- as.numeric(as.character(metrics_components$bp_reference))
metrics_components$total_bp <- as.numeric(as.character(metrics_components$total_bp))
metrics_components$Precision_bp <- metrics_components$bp_ref_component/metrics_components$total_bp



write.table(x = metrics_components, 
            file = snakemake@output[["metrics"]],
            append = TRUE, 
            row.names = FALSE, 
            quote = FALSE, 
            col.names = TRUE)
# 
# truth_set <- subset(truth_set, !truth_set$Contig_number %in% c('representation',''))
# 
# all_truth <- NULL
# 
# for(seed in truth_set$Contig_number)
# {
#   contig_seed <<- subset(truth_set, truth_set$Contig_number == seed)
#   other_contigs <<- subset(truth_set, truth_set$Contig_number != seed)
#   
#   evaluation <- data.frame(Contig_number = seed,
#                            Pair_contig = other_contigs$Contig_number,
#                            Contig_number_type = contig_seed$Type,
#                            Pair_contig_type = other_contigs$Type)
#   
#   
#   evaluation$Truth <- ifelse(as.character(evaluation$Contig_number_type) == as.character(evaluation$Pair_contig_type), 'Same','Different')
#   
#   all_truth <- rbind(all_truth, evaluation)
# }
# 
# collection_pairs <- NULL
# 
# for(seed in 1:nrow(all_truth))
# {
#   info_row <- all_truth[seed,]
#   first_node <- as.numeric(as.character(info_row$Contig_number))
#   second_node <- as.numeric(as.character(info_row$Pair_contig))
#   
#   if(second_node < first_node)
#   {
#     contig_pair <- paste(second_node,first_node, sep = '-')
#   }
#   else
#   {
#     contig_pair <- paste(first_node,second_node, sep = '-')
#   }
#   
#   collection_pairs <- append(collection_pairs, values = contig_pair, after = length(collection_pairs))
#   
# }
# 
# all_truth$Contig_pair <- collection_pairs
# 
# all_truth <- all_truth[!duplicated(all_truth$Contig_pair),]
# 
# ## Extracting the number of connections from each plasmid and chromosome
# 
# true_connections <- subset(all_truth, all_truth$Contig_number_type == all_truth$Pair_contig_type)
# 
# count_true_connections <- true_connections %>%
#   group_by(Contig_number_type) %>%
#   count()
# 
# all_evaluation <- NULL
# 
# for(seed in results_subgraph$Contig_number)
# {
#   contig_seed <- subset(results_subgraph, results_subgraph$Contig_number == seed)
#   other_contigs <- subset(results_subgraph, results_subgraph$Contig_number != seed)
#   
#   evaluation <- data.frame(Contig_number = seed,
#                            Pair_contig = other_contigs$Contig_number,
#                            Contig_number_component = contig_seed$Component,
#                            Pair_contig_component = other_contigs$Component)
#   
#   
#   evaluation$Prediction <- ifelse(as.character(evaluation$Contig_number_component) == as.character(evaluation$Pair_contig_component), 'Same','Different')
#   
#   all_evaluation <- rbind(all_evaluation, evaluation)
# }
# 
# all_evaluation$Contig_pair <- paste(all_evaluation$Contig_number, all_evaluation$Pair_contig, sep = '-' )
# 
# 
# # Creating a confusion matrix 
# 
# eval_truth <<- merge(all_evaluation, all_truth, by = 'Contig_pair')
# 
# eval_truth$Evaluation <- ifelse(eval_truth$Truth == eval_truth$Prediction, 'Positive', 'Negative')
# 
# positive_eval_truth <- subset(eval_truth, eval_truth$Evaluation == 'Positive')
# positive_eval_truth$Evaluation <- ifelse(positive_eval_truth$Truth == 'Same', 'True_positive','True_negative')
# 
# negative_eval_truth <- subset(eval_truth, eval_truth$Evaluation == 'Negative')
# negative_eval_truth$Evaluation <- ifelse(negative_eval_truth$Truth == 'Different', 'False_positive','False_negative')
# 
# eval_truth <- rbind(positive_eval_truth, negative_eval_truth)
# 
# purity_components <- eval_truth %>% 
#   group_by(eval_truth$Contig_number_component, eval_truth$Evaluation) %>%
#   count(Pair_contig_component)
# 
# purity_components <- subset(purity_components, purity_components$`eval_truth$Evaluation` %in% c('True_positive','False_positive'))
# 
# 
# 
# precision_components <- NULL
# 
# for(component in unique(purity_components$Pair_contig_component))
# {
#   component_df <- subset(purity_components, purity_components$Pair_contig_component == component)
#   test <<- subset(component_df, component_df$`eval_truth$Evaluation` == 'True_positive')
#   
#   tp_component <- subset(component_df$n, component_df$`eval_truth$Evaluation` == 'True_positive')
#   fp_component <- subset(component_df$n, component_df$`eval_truth$Evaluation` == 'False_positive')
#   
#   if(identical(fp_component, integer(0)))
#   {
#     fp_component <- 0
#   }
#   
#   if(identical(tp_component, integer(0)))
#   {
#     tp_component <- 0
#   }
#   
#   precision <- tp_component/(tp_component+fp_component)
#   
#   record_results <- data.frame( sample = sample,
#                                 species = species,
#                                 classifier = classifier, 
#                                 factor = 1.0,
#                                 number_iterations = number_iterations,
#                                 variation = max_variation,
#                                 component = component, 
#                                 tp_component = tp_component,
#                                 fp_component = fp_component,
#                                 precision = precision)
#   
#   
#   precision_components <- rbind(precision_components, record_results)
#   
# }
# 
# write.table(x = precision_components, 
#             file = snakemake@output[["precision"]],
#             append = TRUE, 
#             row.names = FALSE, 
#             quote = FALSE, 
#             col.names = TRUE)
# 
