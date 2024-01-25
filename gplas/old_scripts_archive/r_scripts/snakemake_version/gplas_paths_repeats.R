#!/usr/bin/Rscript
# Getting the arguments out of SnakeMake

# Libraries required to generate the gplas output 
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(tidyverse))
suppressMessages(library(spatstat))
suppressMessages(library(cooccur))
suppressMessages(library(ggrepel))

# Inputs
path_nodes <- snakemake@input[["nodes"]]
path_links <- snakemake@input[["clean_links"]]
path_cov_variation <- snakemake@input[["coverage"]]
path_graph_contigs <- snakemake@input[["graph_contigs"]]
path_graph_repeats <- snakemake@input[["graph_repeats"]]
path_prediction <- snakemake@input[["clean_prediction"]]
path_init_nodes <- snakemake@input[["repeat_nodes"]]

bold_sd_coverage <- snakemake@params[["bold_sd_coverage"]]
classifier <- snakemake@params[["classifier"]]
number_iterations <- snakemake@params[["iterations"]]
filtering_threshold <- as.numeric(as.character(snakemake@params[["filt_gplas"]]))

links <- read.table(file = path_links, header = TRUE)
graph_contigs <- read.table(file = path_graph_contigs, header = TRUE)

small_contigs <- subset(graph_contigs, graph_contigs$length < 500)

clean_pred <- read.table(file = path_prediction, header = TRUE)

repeats_graph <- read.table(file = path_graph_repeats, header=TRUE) ###CHECK ABOUT THE HEADER
repeats <- repeats_graph
repeats$number <- gsub(pattern = '\\+',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)
repeats$number <- gsub(pattern = '\\-',replacement = '',x = repeats$number) # Removing directionality (to match the numbers present in mlplasmids prediction)

#export the repeat number
repeats_export<-as.data.frame(unique(repeats$number))

initialize_nodes <- read.table(file = path_init_nodes, header = TRUE)
if( dim(initialize_nodes)[1] == 0){
stop("There are no repeats in your genome. Repeat recovery step will not run")
}

initialize_nodes <- initialize_nodes[,1]

max_variation <- read.table(file = path_cov_variation, header = TRUE)
max_variation <- as.numeric(max_variation[1,1]*as.numeric(bold_sd_coverage))

#Outputs
output_path <- snakemake@output[["solutions"]]
output_connections <- snakemake@output[["connections"]]

plasmid_graph <- function(nodes = nodes, links = links, output_path, classifier, verbose = TRUE, number_iterations = number_iterations, number_nodes = 20, initial_seed, max_variation = max_variation, prob_small_repeats = prob_small_repeats, filtering_threshold = filtering_threshold, direction, classification)
{
  initial_links <- links #Links to start the function 
  
  for(iterations in 1:number_iterations) # Number of times we repeat this process
  {
    path <- NULL # We declare an empty vector
    path <- append(x = path, values = initial_seed, after = length(path)) # We add the initial seed to the path, first element in the vector
    
    if(verbose == TRUE)
    {
      print('INITIALISING THE PATH') # Adding extra information, useful to debug and observe the behaviour of the function 
    }
    
    seed <- initial_seed
    unitig_seed_classification<-classification
    links <- initial_links

    record_connections <- NULL
    record_connections <- data.frame(factor = 1.0,
                                     number_iterations = number_iterations,
                                     iteration = iterations,
                                     elongation = 0, 
                                     first_node = initial_seed, 
                                     ingoing_node = seed, 
                                     outgoing_node = seed, 
                                     Probability_pl_chr = 1,
                                     Probability_cov = 1,
                                     Probability = 1,
                                     Probability_freq = 1,
                                     Veredict = 'selected')
    
    write.table(x = record_connections, 
                file = output_connections, 
                append = TRUE, 
                row.names = FALSE, 
                quote = FALSE, 
                col.names = FALSE)
    
    #################################### Coverage of the current path #################################

    # Extracting the info from the intial node.  
    info_path <- graph_contigs[graph_contigs$number %in% path,]

    length_path <- sum(info_path$length) # Length of the path
    info_path$contribution <- info_path$length/length_path # Shorter contigs should have a lower contribution to the coverage. k-mer coverage on these contigs fluctuates drastically
    
    path_mean <- weighted.mean(x = info_path$coverage, w = info_path$contribution) # Coverage of the current path 
    
    initial_path_mean <- path_mean
    
    ##################### Elongating the path ###########################################
    
    for(elongation in 1:number_nodes)
    {
      evaluation <- 'nope' # Default integer assigned to the evaluation, if different then exit the loop
      if(verbose == TRUE)
      {
        print(paste(elongation,'elongation:','possible connections'))
      }
      
      #get all the connections (degree 1) from the repeat
      current_links <- subset(links, links$V1 == seed) # Consider the last element present in our path and observe all the possible links
     
      if(verbose == TRUE)
      {
        print("#====CURRENT LINKS")
        print(current_links)
      }
      
      if(nrow(current_links) == 0) #if there are not connections, finish the path
      {

        if(verbose == TRUE)
        {
          print(path, quote = FALSE) 
        }
        output <- paste(path, collapse = ',')
        output <- paste(output,classification, unitig_seed_classification, path_mean, collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed # There are no connections possible from this contig
        unitig_seed_classification <- classification
        break # Exiting the loop 
      }
      
      list_connections <- unique(current_links$V3) # All the possible connections 

      # We do not allow that a node which is not a repeat appears more than 1 time in any solution but we exclude the initial seed from this consideration 
      
      if(length(path) != 1 | classification!='Plasmid') # If the path has more than one element  
      {
        remove_nodes <- path[2:length(path)]
      
        remove_nodes <- gsub(x = remove_nodes, pattern = '\\+', replacement = '')
        remove_nodes <- gsub(x = remove_nodes, pattern = '\\-', replacement = '')
        
        positive_remove_nodes <- paste(remove_nodes, '+', sep = '')
        negative_remove_nodes <- paste(remove_nodes, '-', sep = '')
        
        ommit_nodes <- c(positive_remove_nodes, negative_remove_nodes)
        
        first_node <- path[1]
        
        if(direction == 'forward')
        {
          first_node_to_exclude <- gsub(pattern = '\\+', replacement = '-', x = first_node)
        }
        
        if(direction == 'reverse')
        {
          first_node_to_exclude <- gsub(pattern = '-', replacement = '+', x = first_node)
        }
        
        ommit_nodes <- c(ommit_nodes, first_node_to_exclude)
        
        # We need to remove directionality from the path to avoid paths e.g. 54-,161-,54+
        
        ommit_nodes <- ommit_nodes[! ommit_nodes %in% repeats_graph$number] # Repeats can occurr more than one time 
        list_connections <- list_connections[! list_connections %in% ommit_nodes] # Avoiding loops within the solution
      }
      
      possible_connections <- as.character(list_connections) # Connections that can take place  

      total_connections <- length(possible_connections) # Number of connections
      base_probabilities <- rep(0, total_connections) # Generating a vector with the connection probabilities
      
      #if there are no possible connections. We finish the walk
      if(length(possible_connections) < 1) 
      {
        output <- paste(path,  collapse = ',')
        output <- paste(output,classification, unitig_seed_classification,path_mean,collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed
        unitig_seed_classification<-classification
        break 
      }
      
      ################################# Probabilities of the connections #####################################################
      
      # We generate a dataframe with the connection probabilities
      prob_df <- data.frame(number = possible_connections,
                            prob = base_probabilities)
      
      if(verbose == TRUE)
      {
        print(prob_df)
      }
      
      # Replacing the base probabilities with the probabilities generated by mlplasmids 
      prob_df$number <- gsub(pattern = '\\+',replacement = '',x = prob_df$number)
      prob_df$number <- gsub(pattern = '\\-',replacement = '',x = prob_df$number)
      prob_df$number <- as.character(prob_df$number)
      
      subset_probs <- clean_pred[which(clean_pred$number %in% prob_df$number),]
      
      
      if (nrow(subset_probs)>=1){
        if (classification == 'Plasmid' | unitig_seed_classification=='Plasmid') {
          final_probs <- subset_probs$Prob_Plasmid[match(prob_df$number,subset_probs$number)]
        } else if (classification=='Chromosome' | unitig_seed_classification=='Chromosome') {
          final_probs <- subset_probs$Prob_Chromosome[match(prob_df$number,subset_probs$number)]
        } else {
          final_probs<-rep(NA,nrow(prob_df))
        }
      } else {
        final_probs<-rep(NA,nrow(prob_df))
      }

      # Short contigs do not have a reliable probability, we assign them a predefined probability (passed via the argument 'prob_small_repeats')
      final_probs[is.na(final_probs)] <- prob_small_repeats # OVERLAP!  
      
      # Transposases are also corner-cases for the machine-learning algorithm, we follow the same principle as done with short contigs 
      index_repeats <- which(prob_df$number %in% repeats$number)
      final_probs[index_repeats] <- prob_small_repeats
      
      if(verbose == TRUE)
      {
        print(final_probs)
      }
      
      # Inner loop to evaluate all the sensible connections
      times_sampling <- 0
      
      
      # From the connections we sample a possible connection based on the vector with the probabilities
      
      record_connections <- NULL
      
      record_connections <- data.frame(factor = 1.0,
                                       number_iterations = number_iterations,
                                       iteration = iterations,
                                       elongation = elongation, 
                                       first_node = initial_seed, 
                                       ingoing_node = seed, 
                                       outgoing_node = list_connections, 
                                       Probability_pl_chr = final_probs,
                                       Probability_cov = 0)

      cov_connections_info <- graph_contigs[graph_contigs$number %in% record_connections$outgoing_node,]
      
      #provide a probability of coverage connection of 1.
      cov_connections_info$Probability_cov <- 1
      
      #record this probability in the 'record_connections' dataframe.
      record_connections$Probability_cov <- cov_connections_info$Probability_cov[match(record_connections$outgoing_node, cov_connections_info$number)]
      
      #assign probabilities of cov of 0.5 to repeats and small contigs 
      record_connections$Probability_cov[which(record_connections$outgoing_node %in% repeats_graph$number)] <- prob_small_repeats
      record_connections$Probability_cov[which(record_connections$outgoing_node %in% small_contigs$number)] <- prob_small_repeats
      
      #Calculate the gplas score
      record_connections$Probability <- record_connections$Probability_pl_chr * record_connections$Probability_cov
      
      #get the corrected probabilty, which is based on all the potential connections
      record_connections$Probability_freq <- record_connections$Probability/(sum(record_connections$Probability))
      
      #start with all potential connections classified as non-selected.
      record_connections$Veredict <- 'non-selected'
      
      #If there are no potential connections with a probability superior to the threshold, write the path and stop the loop. 
      if(length(which(record_connections$Probability >= filtering_threshold)) == 0)
      {
        output <- paste(path,  collapse = ',')
        output <- paste(output,classification, unitig_seed_classification,path_mean,collapse = ',')
        write.table(x = output, 
                    file = output_path, 
                    append = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE,
                    col.names = FALSE)
        
        path <- initial_seed
        unitig_seed_classification<-classification
        
        write.table(x = record_connections, 
                    file = output_connections, 
                    append = TRUE, 
                    row.names = FALSE, 
                    quote = FALSE, 
                    col.names = FALSE)
        
        break 
      }

      filter_connections <- subset(record_connections, record_connections$Probability >= filtering_threshold)

      #randomly chose one of the potential connections and mark this as selected
      random_connection <- sample(x = filter_connections$outgoing_node, size = 1, prob = filter_connections$Probability) # Choose one connection 
      record_connections$Veredict[which(record_connections$outgoing_node == random_connection)] <- 'selected'
      
      write.table(x = record_connections, 
                  file = output_connections, 
                  append = TRUE, 
                  row.names = FALSE, 
                  quote = FALSE, 
                  col.names = FALSE)
      
      #If the selected connection is the same as the initial seed, stop the walk. This is because we have a circular path.
      if(random_connection == path[1] & classification=='Plasmid')
      {
        path <- append(values = as.character(random_connection), after = length(path), x = path)
        output <- paste(path,  collapse = ',')
        output <- paste(output,classification, unitig_seed_classification,path_mean,collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed
        unitig_seed_classification<-classification
        break 
      }
      
      #Add selected random connection to the path
      path <- append(values = as.character(random_connection), after = length(path), x = path)
      #Use the last added node as seed
      seed <- as.character(random_connection) 
      
      #updating the info_path variable.
      info_path <- graph_contigs[graph_contigs$number %in% path,]
      
      #get a seed clean, with no directions. This will be used to re-classify path in case that we are starting from a Repeat
      seed_clean<-as.character(seed)
      seed_clean<-gsub(x=seed_clean,pattern = '\\+',replacement = '')
      seed_clean<-gsub(x=seed_clean,pattern = '\\-',replacement = '')
      
      #In case that we've started the walks from a repeated element, create a new classification for the path (If the last seed has one). Otherwise, keep calling it a repeat.
      if ((unitig_seed_classification=='Repeat' | unitig_seed_classification=='Small')  & nrow(clean_pred[clean_pred$number %in% seed_clean,])>0){
        unitig_seed_clean<-as.character(seed)
        unitig_seed_clean<-gsub(x=unitig_seed_clean,pattern = '\\+',replacement = '')
        unitig_seed_clean<-gsub(x=unitig_seed_clean,pattern = '\\-',replacement = '')
        unitig_seed_info<-clean_pred[ clean_pred$number %in% unitig_seed_clean,]
        unitig_seed_classification<-as.character(unitig_seed_info$Prediction)
        
      } else if (unitig_seed_classification=='Repeat' & nrow(info_path[info_path$number %in% repeats_graph$number,])!=nrow(info_path)) {
        unitig_seed_classification<-'Small'
      }
      
      # We exit the function if we have reached the maximum number of nodes allowed per path
      if(length(path) == number_nodes) 
      {
        output <- paste(path,  collapse = ',')
        output <- paste(output,classification, unitig_seed_classification,path_mean,collapse = ',')
        write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
        path <- initial_seed
        unitig_seed_classification<-classification
        break 
      }
      
      if (unitig_seed_classification!='Repeat') {
        info_path <- info_path[! info_path$number %in% repeats_graph$number,] # Removing the contigs corresponding to transposases
        #calculate the new path mean coverage.
        length_path <- sum(info_path$length)
        info_path$contribution <- info_path$length/length_path
        path_mean <- weighted.mean(x = info_path$coverage, w = info_path$contribution)
          if (unitig_seed_classification=='Plasmid' | unitig_seed_classification=='Chromosome') {
            output <- paste(path,  collapse = ',')
            output <- paste(output,classification, unitig_seed_classification, path_mean,collapse = ',')
            write.table(x = output, file = output_path, append = TRUE, row.names = FALSE, quote = FALSE, col.names = FALSE)
            path <- initial_seed
            unitig_seed_classification<-classification
            break
          }
      } else {
        first_node<-info_path[ info_path$number %in% initial_seed, ]
        info_path <- info_path[! info_path$number %in% repeats_graph$number,]
        info_path<-rbind(info_path,first_node)
        length_path <- sum(info_path$length)
        info_path$contribution <- info_path$length/length_path
        path_mean <- weighted.mean(x = info_path$coverage, w = info_path$contribution)
      }
  
      if(verbose == TRUE)
      {
        print(paste('in this',elongation,'elongation the path finishes like this:'))
        print(path)
        print(paste('and the new path has the updated following coverage:',path_mean))
        cat("\n")
      }
      
    }
    
  }
}

for (seed in initialize_nodes) {
  set.seed(123)
  seed_info<-clean_pred[ clean_pred$number %in% seed, ]
  if (nrow(seed_info)!=0) {
    seed_classification <- as.character(seed_info$Prediction)
  } else {
    seed_classification <- 'Repeat'
  }
  
  positive_seed <- paste(seed, '+', sep = '')
  negative_seed <- paste(seed, '-', sep = '')
  
  plasmid_graph(direction = 'forward', nodes = nodes, links = links, output_path = output_path, initial_seed = positive_seed, number_iterations = number_iterations, verbose = FALSE, prob_small_repeats = 0.5, max_variation = max_variation, classifier = classifier, filtering_threshold = filtering_threshold, classification = seed_classification)
  plasmid_graph(direction = 'reverse', nodes = nodes, links =  links, output_path = output_path, initial_seed = negative_seed, number_iterations = number_iterations, verbose = FALSE, prob_small_repeats = 0.5, max_variation = max_variation, classifier = classifier, filtering_threshold = filtering_threshold, classification = seed_classification)
}
  




