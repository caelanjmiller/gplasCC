#!/usr/bin/env python3

#from collections import defaultdict
#from email.policy import default
import pandas as pd
import numpy as np
import scipy.stats
import statistics
#import igraph
#import logging
#from multiprocessing import Pool
#from functools import partial
#import copy
from Bio.SeqIO.FastaIO import SimpleFastaParser

manual_mode = True
if(manual_mode == True):
    #Inputs
    path_nodes = "gplas_input\\ecoli_raw_nodes.fasta"
    path_links = "gplas_input\\ecoli_raw_links.txt"
    path_prediction = "test_ecoli_plasmid_prediction.tab"
    #Params
    classifier = "predict"
    threshold = float("0.95")
    #Outputs
    output_graph_contigs = "coverage\\ecoli_graph_contigs.tab"
    output_clean_links = "coverage\\ecoli_clean_links.tab"
    output_graph_repeats = "coverage\\ecoli_repeats_graph.tab"
    output_clean_prediction = "coverage\\ecoli_clean_prediction.tab"
    output_isolated_nodes = "coverage\\ecoli_isolated_nodes.tab"
    output_clean_repeats = "coverage\\ecoli_clean_repeats.tab"
    output_initialize_nodes = "coverage\\ecoli_initialize_nodes.tab"
    output_repeat_nodes = "coverage\\ecoli_repeat_nodes.tab"
    output_cov_estimate = "coverage\\ecoli_estimation.txt"
else:
    #Inputs
    path_nodes = str(snakemake.input["nodes"])
    path_links = str(snakemake.input["links"])
    path_prediction = str(snakemake.input["prediction"])
    #Params
    classifier = str(snakemake.params["classifier"])
    threshold = float(snakemake.params["threshold"])
    #Outputs
    output_graph_contigs = str(snakemake.output["graph_contigs"])
    output_clean_links = str(snakemake.output["clean_links"])
    output_graph_repeats = str(snakemake.output["graph_repeats"])
    output_clean_prediction = str(snakemake.output["clean_prediction"])
    output_isolated_nodes = str(snakemake.output["isolated_nodes"])
    output_clean_repeats = str(snakemake.output["clean_repeats"])
    output_initialize_nodes = str(snakemake.output["initialize_nodes"])
    output_repeat_nodes = str(snakemake.output["repeat_nodes"])
    output_cov_estimate = str(snakemake.output["coverage"])

with open(path_nodes) as file:
    raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]

raw_contig_names = [str(entry[0]) for entry in raw_nodes]

kc_check = sum([name.count("KC") for name in raw_contig_names])

if(kc_check == len(raw_contig_names)):
  lengths = [len(entry[1]) for entry in raw_nodes]
  kc_counts = [name.split(":", maxsplit=4)[2] for name in raw_contig_names]
  kc_counts = [int(name.replace("_","")) for name in kc_counts]
  kc_coverage = [kc/length for kc, length in zip(kc_counts, lengths)]
  coverage = [coverage/statistics.median(kc_coverage) for coverage in kc_coverage]

#improve range 1-len(names)? + cast as string
raw_number = [name.split("_")[0] for name in raw_contig_names]
number = [name.replace("S","") for name in raw_number]

if(kc_check != len(raw_contig_names)):
    raw_lengths = [name.split(":")[2] for name in raw_contig_names]
    lengths = [int(name.replace("_dp","")) for name in raw_lengths]
    coverage = [float(name.split(":")[4]) for name in raw_contig_names]

contig_info = pd.DataFrame(data={"number":number,
                                 "length":lengths,
                                 "coverage":coverage,
                                 "Contig_name":raw_contig_names})

graph_pos_contigs = pd.DataFrame(data={"number":[digit+"+" for digit in number],
                                       "length":lengths,
                                       "coverage":coverage,
                                       "Contig_name":raw_contig_names})

graph_neg_contigs = pd.DataFrame(data={"number":[digit+"-" for digit in number],
                                       "length":lengths,
                                       "coverage":coverage,
                                       "Contig_name":raw_contig_names})
#improve find a way to do this without using concat; just use dicts as intermediary instead of full dataframes?
#improve do we need to use ignore_index=True here?
graph_contigs = pd.concat([graph_pos_contigs, graph_neg_contigs])

graph_contigs.to_csv(output_graph_contigs, sep="\t", index=False)

#improve do we ever use small_contigs?/ limit is hardcoded
small_contigs = graph_contigs[graph_contigs["length"] < 500]

raw_links = pd.read_csv(path_links, sep="\t", header=None)
raw_links.rename({0:'L',1:'first_node',2:'first_sign',3:'second_node',4:'second_sign',5:'OM'}, axis='columns', inplace=True)

links = []
for index, row in raw_links.iterrows():
    reverse_first_sign = "+" if row["first_sign"] == "-" else "-"
    reverse_second_sign = "+" if row["second_sign"] == "-" else "-"

    clean_first = str(row["first_node"])+str(row["first_sign"])
    clean_second = str(row["second_node"])+str(row["second_sign"])
    info_forward = [clean_first, "to", clean_second]
    links.append(info_forward)

    clean_rev_first = str(row["second_node"])+str(reverse_second_sign)
    clean_rev_second = str(row["first_node"])+str(reverse_first_sign)
    info_reverse = [clean_rev_first, "to", clean_rev_second]
    links.append(info_reverse)

links = pd.DataFrame(links)
links.to_csv(output_clean_links, sep="\t", index=False, header=False)

unique_nodes = list(set(links[0]))
#improve skip outdegree_info as empty one and make it in one go with node_info; DONT USE CONCAT
outdegree_info = pd.DataFrame()
for node in unique_nodes:
    repeat_links = links[links[0] == node]
    unique_links = repeat_links.drop_duplicates()
    node_info = pd.DataFrame(data={"number":node,
                                   "outdegree":len(repeat_links[2])}, index = [0])
    #improve do we need to use ignore_index=True here?
    outdegree_info = pd.concat([outdegree_info, node_info])
    
#improve skip indegree_info as empty one and make it in one go with node_info; DONT USE CONCAT
indegree_info = pd.DataFrame()
for node in unique_nodes:
    repeat_links = links[links[2] == node]
    unique_links = repeat_links.drop_duplicates()
    node_info = pd.DataFrame(data={"number":node,
                                   "indegree":len(repeat_links[0])}, index = [0])
    #improve do we need to use ignore_index=True here?
    indegree_info = pd.concat([indegree_info, node_info])

repeat_info = pd.merge(outdegree_info, indegree_info, on="number")
repeats = repeat_info[(repeat_info["indegree"] > 1) | (repeat_info["outdegree"] > 1)]

repeats.to_csv(output_graph_repeats, sep="\t", index=False)

repeats.loc[:,"number"] = [name.replace("+","") for name in repeats["number"]]
repeats.loc[:,"number"] = [name.replace("-","") for name in repeats["number"]]

pred = pd.read_table(path_prediction, sep="\t", header=0)

#improve only check for plasflow if all others are the same, which is plasmidCC?
if classifier == "mlplasmids":
    clean_pred = pred

if classifier == "predict":
    clean_pred = pred

#improve remove code, we're not using plasflow anymore; plasflow code: Rscript lines 154-200
'''
if classifier == "plasflow":
    pred = pred[pred["contig_length"] > 1000]
    for index, contig in pred.iterrows():
        val = max(float(contig.iloc[:,5:32]))
        column = contig[]
        contig.loc[:,"label"] = 0
        contig.loc[:,"prob"] = 0
        contig.loc[:,"Prediction"] = 0
'''
raw_number = [name.split("_", maxsplit=1)[0] for name in clean_pred["Contig_name"]]
clean_pred.loc[:,"number"] = [number.replace("S","") for number in raw_number]
clean_pred = pd.merge(clean_pred, contig_info, on=["Contig_name","number"])

final_prediction = clean_pred[[number not in list(repeats["number"]) for number in clean_pred["number"]]]

final_prediction.to_csv(output_clean_prediction, sep="\t", index=False)

unique_nodes_signless = links[0]
unique_nodes_signless = [node.replace("+","") for node in unique_nodes_signless]
unique_nodes_signless = [node.replace("-","") for node in unique_nodes_signless]

#TODO run with test data that contains isolated nodes
isolated_nodes = contig_info[[number not in unique_nodes_signless for number in contig_info["number"]]]
isolated_nodes = pd.merge(clean_pred, isolated_nodes["Contig_name"], on="Contig_name")
isolated_nodes = isolated_nodes[isolated_nodes["Prob_Plasmid"] >= threshold]

isolated_nodes.to_csv(output_isolated_nodes, sep="\t", index=False)

repeats_final = clean_pred[[number in list(repeats["number"]) for number in clean_pred["number"]]]

repeats_final.to_csv(output_clean_repeats, sep="\t", index=False)

pl_nodes = final_prediction[final_prediction["Prob_Plasmid"] >= threshold]
pl_nodes = pl_nodes[pl_nodes["Contig_length"] > 500]
pl_nodes = pl_nodes[[number not in list(repeats["number"]) for number in pl_nodes["number"]]]
#improve is sorting nessecary here? it is already sorted on number which is based on length
#and gets unsorted again(?) when transferred to initialize_nodes
pl_nodes.sort_values(by="length", axis=0, ascending=False ,inplace=True)

initialize_nodes = list(set(pl_nodes["number"]))

with open(output_initialize_nodes, "w") as file:
    file.write("\n".join(initialize_nodes))

repeats_nodes = list(set(repeats_final["number"]))

with open(output_repeat_nodes, "w") as file:
    file.write("\n".join(repeats_nodes))

chr_contigs = clean_pred[(clean_pred["Prob_Chromosome"] > 0.7) & (clean_pred["Contig_length"] > 1000)]
cov_estimation = chr_contigs[[number not in list(repeats["number"]) for number in chr_contigs["number"]]]

#improve cov calculation with mad is different from Rscript
#same answer if scale=1 but different for scale=1.4826 (R default = 1.4826)
#do we even use cov/max_variation/max_variation_small?
sd_estimation = statistics.stdev(cov_estimation["coverage"])
cov = scipy.stats.median_abs_deviation(cov_estimation["coverage"], scale=1.4826)
max_variation = cov*1.0
max_variation_small = cov*5.0

with open(output_cov_estimate, "w") as file:
    file.write(str(sd_estimation))

#improve remove prints
print(cov)
print(sd_estimation)
print(max_variation)
print(max_variation_small)

