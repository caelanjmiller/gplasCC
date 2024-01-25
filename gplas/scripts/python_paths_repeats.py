#!/usr/bin/env python3

#from collections import defaultdict
#from email.policy import default
import pandas as pd
import numpy as np
#import scipy.stats
#import statistics
#import igraph
#import logging
#from multiprocessing import Pool
#from functools import partial
#import copy
#from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys

#TODO run this script with test data where some initial nodes are classified as either plasmid/chromosome

manual_mode = True
if(manual_mode == True):
    import os
    os.chdir("C:\\Users\\oscar\\Documenten\\UU\\04_BiBc\\6.2_Research_Profile\\gplas\\test_runs\\test_sandbox\\")
    os.getcwd()
    #Inputs
    path_nodes = "gplas_input\\ecoli_raw_nodes.fasta"
    path_links = "coverage\\ecoli_clean_links.tab"
    path_prediction = "coverage\\ecoli_clean_prediction.tab"
    path_graph_contigs = "coverage\\ecoli_graph_contigs.tab"
    path_graph_repeats = "coverage\\ecoli_repeats_graph.tab"
    path_init_nodes = "coverage\\ecoli_repeat_nodes.tab"
    path_cov_variation = "coverage\\ecoli_estimation.txt"
    #Params
    bold_sd_coverage = float("5")
    classifier = "predict"
    number_iterations = int("20")
    filtering_threshold = float("0.1")
    #Outputs
    output_path = "walks\\repeats\\ecoli_solutions.tab"
    output_connections = "walks\\repeats\\ecoli_connections.tab"
else:
    #Inputs
    path_nodes = str(snakemake.input["nodes"])
    path_links = str(snakemake.input["clean_links"])
    path_prediction = str(snakemake.input["clean_prediction"])
    path_graph_contigs = str(snakemake.input["graph_contigs"])
    path_graph_repeats = str(snakemake.input["graph_repeats"])
    path_init_nodes = str(snakemake.input["repeat_nodes"])
    path_cov_variation = str(snakemake.input["coverage"])
    
    #Params
    bold_sd_coverage = float(snakemake.params["bold_sd_coverage"])
    classifier = str(snakemake.params["classifier"])
    number_iterations = int(snakemake.params["iterations"])
    filtering_threshold = float(snakemake.params["filt_gplas"])
    
    #Outputs
    output_path = str(snakemake.output["solutions"])
    output_connections = str(snakemake.output["connections"])

links = pd.read_csv(path_links, sep="\t", header=None)
clean_pred = pd.read_csv(path_prediction, sep="\t", header=0)
clean_pred.loc[:,"number"] = [str(number) for number in clean_pred.loc[:,"number"]] # convert "number" back to string

graph_contigs = pd.read_csv(path_graph_contigs, sep="\t", header=0)

#improve find a way to check both signed and unsigned nodes without making a copy of small/repeat nodes df
small_contigs = graph_contigs[graph_contigs["length"] < 500].copy()
small_contigs_signed_nodes = small_contigs.copy()
small_contigs_signed_nodes = [node for node in small_contigs_signed_nodes["number"]]
small_contigs.loc[:,"number"] = [name.replace("+","") for name in small_contigs["number"]]
small_contigs.loc[:,"number"] = [name.replace("-","") for name in small_contigs["number"]]

repeats = pd.read_csv(path_graph_repeats, sep="\t", header=0)
repeats_signed_nodes = repeats.copy()
repeats_signed_nodes = [node for node in repeats_signed_nodes["number"]]
repeats.loc[:,"number"] = [name.replace("+","") for name in repeats["number"]]
repeats.loc[:,"number"] = [name.replace("-","") for name in repeats["number"]]

#improve we dont use repeats_export?
#export the repeat number
repeats_export = pd.DataFrame(list(set(repeats.loc[:,"number"])))

initialize_nodes = pd.read_csv(path_init_nodes, sep="\t", header=None)

if(initialize_nodes.shape[0] == 0):
    sys.exit("There are no repeats in your genome. Repeat recovery step will not run")

initialize_nodes = [str(node) for node in initialize_nodes.iloc[:,0]]

#improve: just read value from file instead of converting to df and converting to float? speed diff is basically identical
max_variation = pd.read_csv(path_cov_variation, header=None)
max_variation = float(max_variation.iloc[0,0]) * bold_sd_coverage

#ASK check default values of 'nodes' & 'prob_small_repeats'; they are not defined in this script?
# "nodes" is just not defined or even used at all anywhere in the script; I removed it for now
def plasmid_graph(output_path, classifier, initial_seed, direction, prob_small_repeats, classification, links=links, verbose=True, number_iterations=number_iterations, number_nodes=20, max_variation=max_variation, filtering_threshold=filtering_threshold):
    for iteration in range(number_iterations): # Number of times we repeat this process
        path = [initial_seed] # We add the initial seed to the path, first element in the list
        seed = initial_seed
        unitig_seed_classification = classification

        
#improve make sure all instances of record_connections have the same columns, if it is used in later scripts
        record_connections = pd.DataFrame(data={"factor":1.0,
                                                "number_iterations":number_iterations,
                                                "iteration":iteration,
                                                "elongation":0,
                                                "first_node":initial_seed,
                                                "ingoing_node":seed,
                                                "outgoing_node":seed,
                                                "Probability_pl_chr":1,
                                                "Probability_cov":1,
                                                "Probability":1,
                                                "Probability_freq":1,
                                                "Verdict":"selected"}, index=[0])
        record_connections.to_csv(output_connections, sep="\t", index=False, header=False, mode="a")
        
        #################################### Coverage of the current path #################################
        
        # Extracting the info from our current path                     
        info_path = graph_contigs[[number in path for number in graph_contigs["number"]]].copy()
        length_path = sum(info_path["length"]) # Length of the path
        info_path["contribution"] = info_path["length"]/length_path # Shorter contigs should have a lower contribution to the coverage. k-mer coverage on these contigs fluctuates drastically
        path_mean = np.average(info_path["coverage"], weights=info_path["contribution"]) # Coverage of the current path

        ##################### Elongating the path ###########################################
        
        for elongation in range(number_nodes):
            if(verbose == True):
                print("seed", initial_seed, "iteration", iteration, "elongation", elongation)
            
            index = links.loc[:,0] == seed
            current_links = links.loc[index,:] # Consider the last element present in our path and observe all the possible links
            
            if(current_links.shape[0] == 0):
                if(verbose == True):
                    print("final path:", path, "\n")
                output = ",".join(path) #join path with ',' for later step in coocurrence_repeats script
                output = "\t".join([output, classification, unitig_seed_classification, str(path_mean)])
                with open(output_path, mode="a") as file:
                    file.write(output + "\n")
                path = [initial_seed] # There are no connections possible from this contig
                unitig_seed_classification = classification
                break # Exiting elongation loop
            
            list_connections = list(set(current_links[2])) # All the possible unique connections 

#ASK why is classificatioin != "Plasmid" also here?; why run this if not plasmid but len = 1??; change OR to AND?
            # We do not allow that a node which is not a repeat appears more than 1 time in any solution but we exclude the initial seed from this consideration 
            if((len(path) > 1) | (classification != "Plasmid")): # If the path has more than one element        
                remove_nodes = path[1:]
              
                remove_nodes = [node.replace("+","") for node in remove_nodes]
                remove_nodes = [node.replace("-","") for node in remove_nodes]
                
                positive_remove_nodes = [node + "+" for node in remove_nodes]
                negative_remove_nodes = [node + "-" for node in remove_nodes]
                
                ommit_nodes = positive_remove_nodes + negative_remove_nodes

                first_node = path[0]
                
                if(direction == 'forward'):
                    first_node_to_exclude = first_node.replace("+","-")
                    
                elif(direction == 'reverse'):
                    first_node_to_exclude = first_node.replace("-","+")
                
                ommit_nodes.append(first_node_to_exclude)
                
                # We need to remove directionality from the path to avoid paths e.g. 54-,161-,54+
                ommit_nodes = [node for node in ommit_nodes if node not in repeats_signed_nodes]
                
                list_connections = [connection for connection in list_connections if connection not in ommit_nodes] # Avoiding loops within the solution
#improve possible_connections can just be list_connections with str()
            possible_connections = [str(connection) for connection in list_connections] # Connections that can take place
            total_connections = len(possible_connections) # Number of connections
            base_probabilities = [float(0)]*total_connections # Generating a list with the connection probabilities
            
            if(total_connections < 1):
                if(verbose == True):
                    print("final path:", path, "\n")
                output = ",".join(path) #join path with ',' for later step in coocurrence_repeats script
                output = "\t".join([output, classification, unitig_seed_classification, str(path_mean)])
                with open(output_path, mode="a") as file:
                    file.write(output + "\n")
                path = [initial_seed] # There are no connections possible from this contig
                unitig_seed_classification = classification
                break # Exiting elongation loop
                
            ################################# Probabilities of the connections #####################################################
            
            # We generate a dataframe with the connection probabilities
            prob_df = pd.DataFrame(data={"number":possible_connections,
                                         "Prob_Plasmid":base_probabilities,
                                         "Prob_Chromosome":base_probabilities})
 
            # Replacing the base probabilities with the probabilities generated by mlplasmids 
            prob_df.loc[:,"number"] = [number.replace("+","") for number in prob_df["number"]]
            prob_df.loc[:,"number"] = [number.replace("-","") for number in prob_df["number"]]

#ASK check if this makes sense as a fix
#improve there has to be a better way to extract the predictions and keep the same order as in prob_df
#ASK base probs is now 0.5 check if this is 'allowed' or if something needs to change
            if(((classification != "Plasmid") & (classification != "Chromosome")) | ((unitig_seed_classification != "Plasmid") & (unitig_seed_classification != "Chromosome"))):
                final_probs = pd.Series([float(0.5)]*total_connections)
            else:
                for number in prob_df.loc[:,"number"]:
                    index = [value == number for value in clean_pred.loc[:,"number"]]
                    if(clean_pred.loc[index,"Prob_Plasmid"].shape[0] > 0):
                        prob_df.loc[prob_df.loc[:,"number"] == number, "Prob_Plasmid"] = float(clean_pred.loc[index,"Prob_Plasmid"].values[0])
                        prob_df.loc[prob_df.loc[:,"number"] == number, "Prob_Chromosome"] = float(clean_pred.loc[index,"Prob_Chromosome"].values[0])
                
                if((classification == 'Plasmid') | (unitig_seed_classification == 'Plasmid')):
                    final_probs = prob_df.loc[:,"Prob_Plasmid"].copy()
                elif((classification=='Chromosome') | (unitig_seed_classification=='Chromosome')):
                    final_probs = prob_df.loc[:,"Prob_Chromosome"].copy()
            
            # Short contigs do not have a reliable probability, we assign them a predefined probability (passed via the argument 'prob_small_repeats')
            index = prob_df.loc[:,"number"].isin(small_contigs.loc[:,"number"])
            final_probs.loc[index] = float(prob_small_repeats) # OVERLAP!  
            
            # Transposases are also corner-cases for the machine-learning algorithm, we follow the same principle as done with short contigs 
            index = prob_df.loc[:,"number"].isin(repeats.loc[:,"number"])
            final_probs.loc[index] = float(prob_small_repeats)
            
#improve fix this mess for record_connections
            record_connections = []
            for i in range(len(final_probs)):
                record_connections.append([1.0, number_iterations, iteration, elongation, initial_seed, seed, list_connections[i], final_probs[i]])
            
            record_connections = pd.DataFrame(data=record_connections, columns=["factor","number_iterations","iteration","elongation","first_node","ingoing_node","outgoing_node","Probability_pl_chr"])
            
            index = graph_contigs.loc[:,"number"].isin(record_connections.loc[:,"outgoing_node"])
            cov_connections_info = graph_contigs.loc[index,:].copy()

            cov_connections_info.loc[:,"Probability_cov"] = float(1)
#improve find a better way to do this merge
            record_connections = record_connections.merge(cov_connections_info[["number","Probability_cov"]], left_on="outgoing_node", right_on="number")
            record_connections = record_connections.drop(columns="number")
            
            index = record_connections.loc[:,"outgoing_node"].isin(repeats_signed_nodes)
            record_connections.loc[index,"Probability_cov"] = float(prob_small_repeats)
            
            index = record_connections.loc[:,"outgoing_node"].isin(small_contigs_signed_nodes)
            record_connections.loc[index,"Probability_cov"] = float(prob_small_repeats)                        
            
            record_connections.loc[:,"Probability"] = record_connections.loc[:,"Probability_pl_chr"] * record_connections.loc[:,"Probability_cov"]
            record_connections.loc[:,"Probability_freq"] = record_connections.loc[:,"Probability"]/sum(record_connections.loc[:,"Probability"])
            record_connections.loc[:,"Verdict"] = "non-selected"
            
            if(sum(record_connections.loc[:,"Probability"] >= filtering_threshold) == 0):
                if(verbose == True):
                    print("final path:", path, "\n")
                output = ",".join(path)  #join path with ',' for later step in coocurrence_repeats script
                output = "\t".join([output, classification, unitig_seed_classification, str(path_mean)])
                with open(output_path, mode="a") as file:
                    file.write(output + "\n")
                
                record_connections.to_csv(output_connections, sep="\t", index=False, header=False, mode="a")
                path = [initial_seed] # There are no connections possible from this contig
                unitig_seed_classification = classification
                break # Exiting elongation loop
                
            # Filter step to avoid going into really bad connections
            index = list(record_connections.loc[:,"Probability"] >= filtering_threshold)
            filter_connections = record_connections.loc[index,:].copy()
            if(sum(filter_connections["Probability_freq"]) != 1): # recalculate probability frequencies after filter
                filter_connections["Probability_freq"] = filter_connections.loc[:,"Probability"]/sum(filter_connections.loc[:,"Probability"])

            random_connection = np.random.choice(filter_connections.loc[:,"outgoing_node"], size=1, p=filter_connections.loc[:,"Probability_freq"])[0] # Choose one connection 
            
            index = record_connections.loc[:,"outgoing_node"] == random_connection
            record_connections.loc[index,"Verdict"] = "selected"
            
            record_connections.to_csv(output_connections, sep="\t", index=False, header=False, mode="a")

            path.append(str(random_connection))
#ASK why do we check for plasmid here?
            if((random_connection == path[0]) & (classification == "Plasmid")):
                if(verbose == True):
                    print("final path:", path, "\n")
                output = ",".join(path) #join path with ',' for later step in coocurrence_repeats script
                output = "\t".join([output, classification, unitig_seed_classification, str(path_mean)])
                with open(output_path, mode="a") as file:
                    file.write(output + "\n")
                path = [initial_seed]
                unitig_seed_classification = classification
                break # Exiting elongation loop
            
            seed = str(random_connection)
            
            #updating the info_path variable.
            index = graph_contigs.loc[:,"number"].isin(path)
            info_path = graph_contigs.loc[index,:].copy()
            
            #get a seed clean, with no directions. This will be used to re-classify path in case that we are starting from a Repeat
            seed_clean = seed
            seed_clean = seed_clean.replace("+","")
            seed_clean = seed_clean.replace("-","")
      
            #In case that we've started the walks from a repeated element, create a new classification for the path (If the last seed has one). Otherwise, keep calling it a repeat.
            if(((unitig_seed_classification == 'Repeat') | (unitig_seed_classification == 'Small')) & (sum(clean_pred.loc[:,"number"].isin([seed_clean])) > 0)):
                unitig_seed_info = clean_pred.loc[clean_pred.loc[:,"number"].isin([seed_clean]),:]
                unitig_seed_classification = unitig_seed_info.loc[:,"Prediction"].values[0]
            elif((unitig_seed_classification == 'Repeat') & (sum(info_path.loc[:,"number"].isin(repeats_signed_nodes)) != info_path.shape[0])):
                unitig_seed_classification = 'Small'
      
            # We exit the function if we have reached the maximum number of nodes allowed per path
            if(len(path) >= number_nodes): 
                if(verbose == True):
                    print("final path:", path, "\n")
                output = ",".join(path) #join path with ',' for later step in coocurrence_repeats script
                output = "\t".join([output, classification, unitig_seed_classification, str(path_mean)])
                with open(output_path, mode="a") as file:
                    file.write(output + "\n")
                path = [initial_seed]
                unitig_seed_classification = classification
                break # Exiting elongation loop
            
            if(unitig_seed_classification != 'Repeat'):
                index = ~info_path.loc[:,"number"].isin(repeats_signed_nodes)
                info_path = info_path.loc[index,:]
                length_path = sum(info_path.loc[:,"length"])
                info_path.loc[:,"contribution"] = info_path.loc[:,"length"] / length_path
                path_mean = np.average(info_path.loc[:,"coverage"], weights=info_path.loc[:,"contribution"]) # Coverage of the current path
                
                if((unitig_seed_classification == 'Plasmid') | (unitig_seed_classification == 'Chromosome')):
                    if(verbose == True):
                        print("final path:", path, "\n")
                    output = ",".join(path) #join path with ',' for later step in coocurrence_repeats script
                    output = "\t".join([output, classification, unitig_seed_classification, str(path_mean)])
                    with open(output_path, mode="a") as file:
                        file.write(output + "\n")
                    path = [initial_seed]
                    unitig_seed_classification = classification
                    break # Exiting elongation loop
            else:
                first_node = info_path.loc[info_path.loc[:,"number"].isin([initial_seed])].copy()
                index = ~info_path.loc[:,"number"].isin(repeats_signed_nodes)
                info_path = info_path.loc[index,:]
                #improve do we need to use ignore_index=True here?
                info_path = pd.concat([info_path, first_node])                
                length_path = sum(info_path.loc[:,"length"])
                info_path.loc[:,"contribution"] = info_path.loc[:,"length"] / length_path
                path_mean = np.average(info_path.loc[:,"coverage"], weights=info_path.loc[:,"contribution"]) # Coverage of the current path
                

for seed in initialize_nodes:
    np.random.seed(123)
    seed_info = clean_pred.loc[clean_pred.loc[:,"number"].isin([seed]),:]
    if(seed_info.shape[0] != 0):
        seed_classification = seed_info.loc[:,"Prediction"]
    else:
        seed_classification = "Repeat"
    
    positive_seed = seed + "+"
    negative_seed = seed + "-"
    
    plasmid_graph(output_path=output_path, classifier=classifier, initial_seed=positive_seed, direction="forward", prob_small_repeats=0.5, classification = seed_classification, links=links, verbose=False, number_iterations=number_iterations, number_nodes=100, max_variation=max_variation, filtering_threshold=filtering_threshold)
    plasmid_graph(output_path=output_path, classifier=classifier, initial_seed=negative_seed, direction="reverse", prob_small_repeats=0.5, classification = seed_classification, links=links, verbose=False, number_iterations=number_iterations, number_nodes=100, max_variation=max_variation, filtering_threshold=filtering_threshold)
