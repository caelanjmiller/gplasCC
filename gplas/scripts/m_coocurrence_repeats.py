#!/usr/bin/env python3

#from collections import defaultdict
#from email.policy import default
import pandas as pd
import numpy as np
#import scipy.stats
#import statistics
import igraph as ig
#import logging
#from multiprocessing import Pool
#from functools import partial
#import copy
from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys

#improve we have a lot of loops with "for row in range(solutions.shape[0]):"
## can we possibly merge some of them?

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
    input_solutions = "walks\\repeats\\ecoli_solutions.tab"
    path_bins = "results\\ecoli_results_no_repeats.tab"
    clean_repeats_path = "coverage\\ecoli_clean_repeats.tab"
    #Params
    classifier = "predict"
    threshold = float("0.95")
    number_iterations = int("20")
    modularity_threshold = float("0.2")
    bold_sd_coverage = float("2")
    sample_name = "ecoli"
    #Outputs
    output_dir = "results\\"
    output_results = "results\\ecoli_results.tab"
    output_components = "results\\ecoli_bins.tab"
    output_chromosomes = "results\\ecoli_chromosome_repeats.tab"
else:
    #Inputs
    path_nodes = str(snakemake.input["nodes"])
    path_links = str(snakemake.input["clean_links"])
    path_prediction = str(snakemake.input["clean_prediction"])
    path_graph_contigs = str(snakemake.input["graph_contigs"])
    path_graph_repeats = str(snakemake.input["graph_repeats"])
    path_init_nodes = str(snakemake.input["repeat_nodes"])
    path_cov_variation = str(snakemake.input["coverage"])
    input_solutions = str(snakemake.input["solutions_repeat"])
    path_bins = str(snakemake.input["bins"])
    clean_repeats_path = str(snakemake.input["clean_repeats"])
    #Params
    classifier = str(snakemake.params["classifier"])
    threshold = float(snakemake.params["threshold"])
    number_iterations = int(snakemake.params["iterations"])
    modularity_threshold = float(snakemake.params["modularity_threshold"])
    #improve rename this to repeat_sd_coverage
    bold_sd_coverage = float(snakemake.params["bold_sd_coverage"])
    sample_name = str(snakemake.params["sample"])
    #Outputs
    output_dir = "results\\"
    output_results = str(snakemake.output["results"])
    output_components = str(snakemake.output["components"])
    output_chromosomes = str(snakemake.output["chromosome_repeats"])

links = pd.read_csv(path_links, sep="\t", header=None)
clean_pred = pd.read_csv(path_prediction, sep="\t", header=0)
clean_pred.loc[:,"number"] = [str(number) for number in clean_pred.loc[:,"number"]] # convert "number" back to string

graph_contigs = pd.read_csv(path_graph_contigs, sep="\t", header=0)

#improve find a way to check both signed and unsigned nodes without making a copy of small/repeat nodes df
##what part of small/repeats / signed/unsigned is used in this/each script?
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

initialize_nodes = pd.read_csv(path_init_nodes, sep="\t", header=None)
initialize_nodes = [str(node) for node in initialize_nodes.iloc[:,0]]

#improve do we use max_variation in this script at all??
#improve: just read value from file instead of converting to df and converting to float? speed diff is basically identical
max_variation = pd.read_csv(path_cov_variation, header=None)
max_variation = float(max_variation.iloc[0,0]) * bold_sd_coverage
#improve do we use this? / it is hardcoded as 5x
max_variation_small = max_variation * 5

solutions = pd.read_csv(input_solutions, sep="\t", header=None, names=["walks", "initial_classification", "unitig_classification", "path_coverage"])
max_nodes = max([walk.count(",")+1 for walk in solutions.loc[:,"walks"]])
steps = ["".join(["step_", str(step)]) for step in range(max_nodes)]
solutions[steps] = solutions.loc[:,"walks"].str.split(',', expand=True)

solutions = solutions.drop(columns="walks")
col_order = solutions.columns.to_list()
col_order = col_order[3:] + col_order[:3]
solutions = solutions[col_order]

#get the last node of each solution
last_nodes = []
for row in range(solutions.shape[0]):
    last_node = [node for node in solutions.iloc[row,0:max_nodes].dropna()][-1]
    last_nodes.append(last_node)

last_nodes_signless = [node.replace("+","") for node in last_nodes]
last_nodes_signless = [node.replace("-","") for node in last_nodes_signless]

solutions.loc[:,"last_nodes"] = last_nodes
solutions.loc[:,"last_nodes_signless"] = last_nodes_signless

#====Merge with bin data=============
#Import information from the bins and match it with total pairs
bins_data = pd.read_csv(path_bins, sep="\t", header=0)
bins_data = bins_data.astype({"Prob_Chromosome":float,
                              "Prob_Plasmid":float,
                              "Prediction":str,
                              "Contig_name":str,
                              "number":str,
                              "length":int,
                              "coverage":float,
                              "Bin":str})
#merge it
solutions = solutions.merge(bins_data.loc[:,["number","Bin"]], how="left", left_on="last_nodes_signless", right_on="number")
solutions = solutions.drop(columns="number")
#assign "C" to the chromosome and -1 to repeats
#ASK repeats are not assigned -1 in R?
index = solutions.loc[:,"Bin"].isna()
solutions.loc[index,"Bin"] = "C"

#remove connections to unbinned unitigs
index = solutions.loc[:,"Bin"] != "Unbinned"
solutions = solutions.loc[index,:]

#remove connections to repeats only
index = solutions.loc[:,"unitig_classification"] != "Repeat"
solutions = solutions.loc[index,:]

#get the initial node without a symbol
initial_nodes = [node.replace("+","") for node in solutions.loc[:,"step_0"]]
initial_nodes = [node.replace("-","") for node in initial_nodes]
solutions.loc[:,"initial_node"] = initial_nodes

#combine inital nodes with bin in a single variable
solutions.loc[:,"repeat_bin"] = ["-".join([solutions.loc[row,"initial_node"], solutions.loc[row,"Bin"]]) for row in range(solutions.shape[0])]

#keep cases in which the last node and second node are the same (repeat directly connected to a unitig)
solutions.loc[:,"keep"] = ["yes" if ((solutions.loc[row,"step_1"] == solutions.loc[row,"last_nodes"]) &
                                     (solutions.loc[row,"unitig_classification"] != "Repeat")) else "no" for row in range(solutions.shape[0])]

index = solutions.loc[:,"keep"] == "yes"
final_solutions = solutions.loc[index,:]

#----Analyze the remainig cases to check if we have the same bin upstream and downstream from the repeat--
#separate into positive and negative solutions
index = [solutions.loc[row,"step_1"][-1] == "+" for row in range(solutions.shape[0])]
positive_solutions = solutions.loc[index,:]
index = [solutions.loc[row,"step_1"][-1] == "-" for row in range(solutions.shape[0])]
negative_solutions = solutions.loc[index,:]

#Keep only the walks if they lead to same bin, or to the chromosome in both directions
index = ((positive_solutions.loc[:,"repeat_bin"].isin(negative_solutions.loc[:,"repeat_bin"])) &
         (positive_solutions.loc[:,"keep"] == "no"))
positive_solutions_keep = positive_solutions.loc[index,:]
index = ((negative_solutions.loc[:,"repeat_bin"].isin(positive_solutions.loc[:,"repeat_bin"])) &
         (negative_solutions.loc[:,"keep"] == "no"))
negative_solutions_keep = negative_solutions.loc[index,:]

## Filter only the valid walks
final_walks = final_solutions.iloc[:,:max_nodes]
final_positive_walks = positive_solutions_keep.iloc[:,:max_nodes]
final_negative_walks = negative_solutions_keep.iloc[:,:max_nodes]

solutions = pd.concat([final_walks, final_positive_walks, final_negative_walks], ignore_index=True)

#create a list with all the nodes that appear in the plasmid walks.
all_nodes = []
for row in range(solutions.shape[0]):
    nodes = [node for node in solutions.loc[row,:].dropna()]
    all_nodes.extend(nodes)

#get unique set of nodes
unique_nodes = list(set(all_nodes))

#CREATE A CO-OCURRENCE MATRIX
##Column names are the nodes included in plasmid walks.
##Each row is a new walk
##Assign True if the node is present in walk and False if node is not present
co_ocurrence = []
for row in range(solutions.shape[0]):
    walk = [node for node in solutions.loc[row,:].dropna()]
    presence_absence = [node in walk for node in unique_nodes]
    co_ocurrence.append(presence_absence)
co_ocurrence = pd.DataFrame(co_ocurrence, columns=unique_nodes)

starting_nodes = [node for node in unique_nodes if node in solutions.loc[:,"step_0"].values]

#improve move scalar function to top of script
def scalar1(x):
    denominator = (sum([value*value for value in x]))**0.5
    scaled_x = [value/denominator for value in x]
    return scaled_x

#create a dataframe for co-ocurrence frequency (in network format)
#Start_node, Connecting_node, nr_occurences
total_pairs = []
#Get the number of times that two nodes co-ocuur in every walk
for node in starting_nodes:
    index_col = [col_name == node for col_name in unique_nodes] # select the column of the target node
    index_walks = co_ocurrence.loc[:,index_col].values # select all walks where the target node is present
    walks = co_ocurrence.loc[index_walks,:].copy()
    col_sums = [sum(walks.loc[:,col]) for col in walks] # count how often any node is present in all selected walks
    for index_connecting_node in range(len(col_sums)): # save coocurrence data for every node/connecting_node
        connecting_node = unique_nodes[index_connecting_node]
        if(connecting_node != node):
            weight = col_sums[index_connecting_node]
            total_pairs.append([node, connecting_node, weight])

"""
#====Find circular sequences
circular_sequences = []
#Extract the walks in which start-node and end-node are the same.
for row in range(solutions.shape[0]):
    walk = [node for node in solutions.loc[row,:].dropna()]
    if((len(walk) > 1) & (walk[0] == walk[-1])):
        circular_sequences.append([walk[0], walk[-1]])

#check if the number of circular walks starting from each node equals the number of iterations.
#if this is the case, add the circular walk to total_pairs
if(len(circular_sequences) > 0):
    no_duplicated = [list(unique_walk) for unique_walk in set(tuple(walk) for walk in circular_sequences)] # remove duplicate entries

    for combination in range(len(no_duplicated)):
        combi = no_duplicated[combination]
        total_ocurrences = sum([walk[1] == combi[1] for walk in circular_sequences])
        if(total_ocurrences == number_iterations):
            total_pairs.append([combi[0], combi[1], total_ocurrences])
"""

total_pairs = pd.DataFrame(total_pairs, columns=["Starting_node", "Connecting_node", "weight"])

total_pairs.loc[:,"Starting_node"] = [node.replace("+","") for node in total_pairs.loc[:,"Starting_node"]]
total_pairs.loc[:,"Starting_node"] = [node.replace("-","") for node in total_pairs.loc[:,"Starting_node"]]
total_pairs.loc[:,"Connecting_node"] = [node.replace("+","") for node in total_pairs.loc[:,"Connecting_node"]]
total_pairs.loc[:,"Connecting_node"] = [node.replace("-","") for node in total_pairs.loc[:,"Connecting_node"]]

#Filter-out cases of no-coocurrence
#ASK filter is on weight > 1 why not weight > 0??
index = [weight > 1 for weight in total_pairs.loc[:,"weight"]]
total_pairs = total_pairs.loc[index,:]

#Scale weights
complete_node_info = pd.DataFrame()
for node in list(set(total_pairs.loc[:,"Starting_node"])):
    index = total_pairs.loc[:,"Starting_node"] == node
    first_node = total_pairs.loc[index,:]
    particular_node = []

    for connecting_node in list(set(first_node.loc[:,"Connecting_node"])):
        index = first_node.loc[:,"Connecting_node"] == connecting_node
        first_second_nodes = first_node.loc[index,:]
        total_weight = sum(first_second_nodes.loc[:,"weight"])
        particular_node.append([node, connecting_node, total_weight])

    particular_node = pd.DataFrame(particular_node, columns=["Starting_node", "Connecting_node", "weight"])
    particular_node.loc[:,"scaled_weight"] = scalar1(particular_node.loc[:,"weight"]) # add a column with scaled weights using scalar1()
    complete_node_info = pd.concat([complete_node_info, particular_node], ignore_index=True)

total_pairs = complete_node_info

initial_nodes = [node.replace("+","") for node in starting_nodes]
initial_nodes = [node.replace("-","") for node in initial_nodes]

#get a list of clean untigs
clean_unitigs = clean_pred.loc[:,"number"]

#Filter out connected repeated elements. Keep only connections from starting nodes (repeats) unitigs.
index = ((total_pairs.loc[:,"Starting_node"].isin(initial_nodes)) &
         (total_pairs.loc[:,"Connecting_node"].isin(clean_unitigs)))
total_pairs = total_pairs.loc[index,:]

#Import information from the bins and match it with total pairs
#improve move bins_coverage code block to better location? (down)
bins_coverage = []
for current_bin in sorted(list(set(bins_data.loc[:,"Bin"]))):
    index = bins_data.loc[:,"Bin"] == current_bin
    bins_subset = bins_data.loc[index,:]
    mean_coverage = round(np.mean(bins_subset.loc[:,"coverage"]), 2)
    bins_coverage.append([current_bin, mean_coverage])

bins_coverage = pd.DataFrame(data=bins_coverage, columns=["Bin", "bin_coverage"])

total_pairs = total_pairs.merge(bins_data.loc[:,["number", "Bin"]], how="left", left_on="Connecting_node", right_on="number")
total_pairs = total_pairs.drop(columns="number")
#ASSIGN "C" to the chromosome
index = total_pairs.loc[:,"Bin"].isna()
total_pairs.loc[index,"Bin"] = "C"

total_pairs = total_pairs.loc[:,["Starting_node", "Bin", "weight", "scaled_weight"]]

#===Reformat the dataframe to obtain the totality of co-courences
#First check if we actually have co-ocurrence of unitigs.
if((total_pairs.shape[0] > 0) & (total_pairs.shape[1] > 0)):
    total_pairs.loc[:,"Pair"] = ["-".join([total_pairs.loc[row,"Bin"], total_pairs.loc[row,"Starting_node"]]) for row in range(total_pairs.shape[0])]
else:
    print("gplas couldn't find any walks connecting repeats to plasmid-nodes.")
    #ASK why is it randomly status=1 in the final script?
    sys.exit(1)

single_edge_counting = []
for pair in sorted(list(set(total_pairs.loc[:,"Pair"]))):
    index = total_pairs.loc[:,"Pair"] == pair
    pairs_subset = total_pairs.loc[index,:]
    sum_weight = sum(pairs_subset.loc[:,"weight"])
    single_edge_counting.append([pair, sum_weight])

pairs = [pair[0].split("-") for pair in single_edge_counting]

weight_graph = pd.DataFrame(data={"From_to":[pair[1] for pair in pairs],
                                  "To_from":[pair[0] for pair in pairs],
                                  "weight":[weight[1] for weight in single_edge_counting]})

total_scaled_weight = []
full_graph_info = pd.DataFrame()
#Get the data from coverages. Repeats and bins
clean_repeats = pd.read_csv(clean_repeats_path, sep="\t", header=0)
clean_repeats = clean_repeats.astype({"number":str,
                                      "coverage":float})

weight_graph = weight_graph.merge(clean_repeats.loc[:,["number", "coverage"]], how="left", left_on="From_to", right_on="number")
#improve find a more elegant/efficient way to properly merge dataframes
weight_graph = weight_graph.drop(columns="number")

weight_graph = weight_graph.merge(bins_coverage, how="left", left_on="To_from", right_on="Bin")
weight_graph = weight_graph.drop(columns="Bin")

#assign a coverage of 1 to chromosomal unitigs
index = weight_graph.loc[:,"bin_coverage"].isna()
weight_graph.loc[index,"bin_coverage"] = float(1)

#Explore if the combination of bins proposed by the algorithm is plausible based on coverage
repeat_assignments = pd.DataFrame()
#loop through each of the repeats
for node in sorted(list(set(weight_graph.loc[:,"From_to"]))):
    index = weight_graph.loc[:,"From_to"] == node
    df_node = weight_graph.loc[index,:].copy()
    #create a rank of the most likely connections, based on the co-ocurence count
    df_node.loc[:,"rank"] = [sorted(df_node.loc[:,"weight"], reverse=True).index(weight) for weight in df_node.loc[:,"weight"]]

    rank = 0 #start from the highest ranking (0)
    accumulated_cov = 0
    while rank <= max(df_node.loc[:,"rank"]):
        index = df_node.loc[:,"rank"] == rank
        repeat_bin = df_node.loc[index,:]
        if(repeat_bin.loc[:,"coverage"].values[0] + accumulated_cov >= repeat_bin.loc[:,"bin_coverage"].values[0]):
            accumulated_cov += repeat_bin.loc[:,"bin_coverage"].values[0]
            repeat_assignments = pd.concat([repeat_assignments, repeat_bin], ignore_index=True)
            rank += 1
        else:
            rank += 1

#ASK this whole rank stuff was just for fun and is immediatly removed again????????????????
repeat_assignments = repeat_assignments.loc[:,["To_from","From_to"]].rename(columns={"To_from":"Bin",
                                                                                     "From_to":"number"})
#separate results into plasmid and chromosome repeats
index = repeat_assignments.loc[:,"Bin"] == "C"
chromosome_repeats = repeat_assignments.loc[index,:]

index = repeat_assignments.loc[:,"Bin"] != "C"
plasmid_repeats = repeat_assignments.loc[index,:]

if(plasmid_repeats.shape[0] == 0):
    print("No repeats associated with plasmids were found")
    bins_data.loc[:,"Prob_Chromosome"] = round(bins_data.loc[:,"Prob_Chromosome"], 2)
    bins_data.loc[:,"Prob_Plasmid"] = round(bins_data.loc[:,"Prob_Plasmid"], 2)
    bins_data.loc[:,"coverage"] = round(bins_data.loc[:,"coverage"], 2)
    full_info_assigned = bins_data
else:
    print("We found repeated elements associated to plasmid predictions")
    #Get all the repeat nodes
    index = clean_repeats.loc[:,"number"].isin(plasmid_repeats.loc[:,"number"]) # Selecting only contigs predicted as plasmid-derived
    pl_nodes = clean_repeats.loc[index,:]

    #Get all the information from the plasmid nodes (Length, coverage, bin number, etc)
    full_info_assigned = pl_nodes.merge(plasmid_repeats, on="number")

    #Add information to the results file
    full_info_assigned = full_info_assigned.drop(columns="Contig_length")
    full_info_assigned.loc[:,"Prob_Chromosome"] = round(full_info_assigned.loc[:,"Prob_Chromosome"], 2)
    full_info_assigned.loc[:,"Prob_Plasmid"] = round(full_info_assigned.loc[:,"Prob_Plasmid"], 2)
    full_info_assigned.loc[:,"coverage"] = round(full_info_assigned.loc[:,"coverage"], 2)
    full_info_assigned.loc[:,"Prediction"] = "Repeat"

    bins_data.loc[:,"Prob_Chromosome"] = round(bins_data.loc[:,"Prob_Chromosome"], 2)
    bins_data.loc[:,"Prob_Plasmid"] = round(bins_data.loc[:,"Prob_Plasmid"], 2)
    bins_data.loc[:,"coverage"] = round(bins_data.loc[:,"coverage"], 2)

    full_info_assigned = pd.concat([full_info_assigned, bins_data], ignore_index=True)

#Create the fasta files
with open(path_nodes) as file:
    raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]

df_nodes = pd.DataFrame(data=raw_nodes, columns=["Contig_name", "Sequence"])
df_nodes = df_nodes.merge(full_info_assigned, on="Contig_name")

#Write fasta files
for component in set(df_nodes.loc[:,"Bin"]): #improve also use sorted() everytime we use a set()
    index = df_nodes.loc[:,"Bin"] == component
    nodes_component = df_nodes.loc[index,:]
    component_complete_name = "_".join([sample_name, "bin", str(component)])
    filename = "".join([output_dir, component_complete_name, ".fasta"])

    with open(filename, mode="w") as file:
        for contig in range(nodes_component.shape[0]):
            file.write(">" + nodes_component.iloc[contig,0] + "\n" + nodes_component.iloc[contig,1] + "\n")

results_summary = df_nodes.loc[:,["number", "Bin"]]

full_info_assigned.to_csv(output_results, sep="\t", index=False, header=True, mode="a")
results_summary.to_csv(output_components, sep="\t", index=False, header=True, mode="a")

#improve change the column order of ecoli_results_no_repeats to match the output of R?
##order node order in 'ecoli_bins_no_repeats' / ecoli_results_no_repeats & co.
### bins_no_repeats is not ordered but results_no_repeats is?

#format chromosome repeats and print
chromosome_repeats = chromosome_repeats.loc[:,["number", "Bin"]] #improve change col order from the start instead of now
chromosome_repeats.loc[:,"Bin"] = "Chromosome" #improve call it "Chromosome" from the start instead of "C"
chromosome_repeats.to_csv(output_chromosomes, sep="\t", index=False, header=True, mode="a")
