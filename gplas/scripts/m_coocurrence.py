import pandas as pd
import numpy as np
import igraph as ig
from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys

#improve we have a lot of loops with "for row in range(solutions.shape[0]):"
## can we possibly merge some of them?
def calculate_coocurrence(sample, number_iterations, pred_threshold, mod_threshold, mode="normal"):
    if mode == "normal":
        subdir = "normal_mode/"
    elif mode == "unbinned":
        subdir = ""
    #Inputs
    path_nodes = f"gplas_input/{sample}_raw_nodes.fasta"
    path_links = f"coverage/{sample}_clean_links.tab"
    path_prediction = f"coverage/{sample}_clean_prediction.tab"
    path_graph_contigs = f"coverage/{sample}_graph_contigs.tab"
    path_graph_repeats = f"coverage/{sample}_repeats_graph.tab"
    path_init_nodes = f"coverage/{sample}_initialize_nodes.tab"
    path_isolated_nodes = f"coverage/{sample}_isolated_nodes.tab"
    path_cov_variation = f"coverage/{sample}_estimation.txt"
    input_solutions = f"walks/{subdir}{sample}_solutions.tab"
    #Params
    #improve can't we get rid of #Params and just use the function arguments?
    number_iterations = int(number_iterations)
    pred_threshold = float(pred_threshold)
    modularity_threshold = float(mod_threshold)
    sample = str(sample)
    #Outputs
    output_dir = f"results/{subdir}"
    output_results = f"results/{subdir}{sample}_results_no_repeats.tab"
    output_components = f"results/{subdir}{sample}_bins_no_repeats.tab"
    output_png = f"results/{subdir}{sample}_plasmidome_network.png"

    links = pd.read_csv(path_links, sep="\t", header=None)
    clean_pred = pd.read_csv(path_prediction, sep="\t", header=0)
    clean_pred = clean_pred.astype({"Prob_Chromosome":float,
                                    "Prob_Plasmid":float,
                                    "Prediction":str,
                                    "Contig_name":str,
                                    "Contig_length":int,
                                    "number":str,
                                    "length":int,
                                    "coverage":float})
    
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
    
    initialize_nodes = pd.read_csv(path_init_nodes, sep="\t", header=None)
    initialize_nodes = [str(node) for node in initialize_nodes.iloc[:,0]]
    
    #improve do we use max_variation in this script at all??
    #improve: just read value from file instead of converting to df and converting to float? speed diff is basically identical
    max_variation = pd.read_csv(path_cov_variation, header=None)
    max_variation = float(max_variation.iloc[0,0])
    #improve do we use this? / it is hardcoded as 5x
    max_variation_small = max_variation * 5
    
    with open(input_solutions) as file:
        max_nodes = max([line.count("\t")+1 for line in file])
    
    solutions = pd.read_csv(input_solutions, sep="\t", header=None, names=range(max_nodes))
    
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
    
    starting_nodes = [node for node in unique_nodes if node in solutions.loc[:,0].values]
    
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
    
    #Filter out repeated elements. Keep only connections of unitigs.
    index = ((total_pairs.loc[:,"Starting_node"].isin(initial_nodes)) &
             (total_pairs.loc[:,"Connecting_node"].isin(initial_nodes)))
    total_pairs = total_pairs.loc[index,:]
    
    #===Reformat the dataframe to obtain the totality of co-courences
    weight_counting = []
    #First check if we actually have co-ocurrence of unitigs.
    if((total_pairs.shape[0] > 0) & (total_pairs.shape[1] > 0)):
        for row in range(total_pairs.shape[0]):
            initial_node = int(total_pairs.iloc[row,0])
            connecting_node = int(total_pairs.iloc[row,1])
            raw_count = int(total_pairs.iloc[row,2])
            #improve move the int() casting to only the if statement and keep is as (default)/str for the rest
            if(initial_node < connecting_node):
                pair = "-".join([str(initial_node), str(connecting_node)])
            else:
                pair = "-".join([str(connecting_node), str(initial_node)])
    
            weight_counting.append([pair, raw_count])
    else:
        print("gplas couldn't find any walks connecting plasmid-predicted nodes. Plasmid nodes will be classified as Unbinned. If this is unexpected, please assemble your genome with different parameters or with a different tool and re-run gplas.")
        index = clean_pred.loc[:,"Prob_Plasmid"] > pred_threshold
        pl_unbinned = clean_pred.loc[index,:].copy()
        pl_unbinned.loc[:,"Component"] = "Unbinned"
        pl_unbinned = pl_unbinned.drop(columns="Contig_length")
        pl_unbinned.loc[:,"Prob_Chromosome"] = round(pl_unbinned.loc[:,"Prob_Chromosome"], 2)
        pl_unbinned.loc[:,"Prob_Plasmid"] = round(pl_unbinned.loc[:,"Prob_Plasmid"], 2)
        pl_unbinned.loc[:,"coverage"] = round(pl_unbinned.loc[:,"coverage"], 2)
    
        with open(path_nodes) as file:
            raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]
    
        df_nodes = pd.DataFrame(data=raw_nodes, columns=["Contig_name", "Sequence"])
        df_nodes = df_nodes.merge(pl_unbinned, on="Contig_name")
    
        for component in set(df_nodes.loc[:,"Component"]):
            index = df_nodes.loc[:,"Component"] == component
            nodes_component = df_nodes.loc[index,:]
            component_complete_name = "_".join([sample, "bin", str(component)])
            filename = "".join([output_dir, component_complete_name, ".fasta"])
    
            with open(filename, mode="w") as file:
                for contig in range(nodes_component.shape[0]):
                    file.write(">" + nodes_component.iloc[contig,0] + "\n" + nodes_component.iloc[contig,1] + "\n")
    #improve name the column 'Bin' from the start instead of 'Component' and then renaming
    #can we use lowercase 'bin'?
        pl_unbinned = pl_unbinned.rename(columns={'Component':'Bin'})
        results_subgraph = pl_unbinned.loc[:,["number","Bin"]]
    
        pl_unbinned.to_csv(output_results, sep="\t", index=False, header=True, mode="w")
        results_subgraph.to_csv(output_components, sep="\t", index=False, header=True, mode="w")
    #improve remove this empty png?? or keep for consistency by always producing 'a' plot?
        ig.plot(None, target=output_png, bbox=(700,700))
        return
    
    weight_counting = pd.DataFrame(weight_counting, columns=["Pair","Count"])
    
    unique_pairs = list(set(weight_counting.loc[:,"Pair"]))
    sum_weights = []
    for pair in unique_pairs:
        index = list(weight_counting.loc[:,"Pair"] == pair)
        sum_weights.append(sum(weight_counting.loc[index,"Count"]))
    
    #improve remove random comma in .split()?
    pairs = [pair.split("-",) for pair in unique_pairs]
    
    weight_graph = pd.DataFrame(data={"From_to":[pair[0] for pair in pairs],
                                      "To_from":[pair[1] for pair in pairs],
                                      "weight":sum_weights})
    
    total_scaled_weight = []
    full_graph_info = pd.DataFrame()
    
    for node in list(set(weight_graph.loc[:,"From_to"])):
        index = weight_graph.loc[:,"From_to"] == node
        df_node = weight_graph.loc[index,:].copy()
        df_node.loc[:,"scaled_weight"] = scalar1(df_node.loc[:,"weight"])
        #improve do we need to use ignore_index=True here?
        full_graph_info = pd.concat([full_graph_info, df_node])
    
    full_graph_info.loc[:,"width"] = full_graph_info.loc[:,"scaled_weight"].copy() * 5
    #improve get rid of this merge and other df shenanigans; current problem is the lost order when taking unique nodes list(set(from_to))
    weight_graph = weight_graph.merge(full_graph_info, on=["From_to", "To_from", "weight"])
    
    """ # improve remove
    # change datatypes to numerical
    weight_graph = weight_graph.astype({"From_to":int,
                                        "To_from":int,
                                        "weight":int,
                                        "scaled_weight":float,
                                        "width":float})
    """
    
    graph_pairs = ig.Graph.DataFrame(weight_graph, directed=False, use_vids=False)
    
    # Simplifying the graph - Removing self-loops from the graph
    no_loops_graph = graph_pairs.simplify(multiple=False)
    
    #===Analyze if the plasmidome network can be partitioned into different sub-networks.
    #improve move components_graph down for better flow
    #Get clusters based on connectivity alone.
    components_graph = no_loops_graph.decompose(mode="weak", minelements=2)
    
    #Create function to run three community detection algorithms on the plasmidome network.
    #This fuction creates a dataframe as an output that contains the different algorithms names and the global modularity values of the resulting network.
    def partitioning_components(graph):
        # Walktrap algorithm
        dendrogram_walktrap = graph.community_walktrap()
        clusters_walktrap = dendrogram_walktrap.as_clustering()
        modularity_walktrap = clusters_walktrap.modularity
    
        # Leading eigen values
        clusters_eigen = graph.community_leading_eigenvector()
        modularity_eigen = clusters_eigen.modularity
    
        # Louvain method
        clusters_louvain = graph.community_multilevel()
        modularity_louvain = clusters_louvain.modularity
    
        # Spinglass algorithm
        ##graph_spin = graph.community_spinglass()
    
        # Community detection based on propagating labels
        ##graph_propag = graph.community_label_propagation()
    
        partition_info = pd.DataFrame(data={"Algorithm":["Walktrap", "Leading-eigen", "Louvain"],
                                            "Modularity":[modularity_walktrap, modularity_eigen, modularity_louvain]})
        return partition_info
    
    #Create a data-frame to hold the results from the different clusters
    #improve relocate complete_partition_info
    complete_partition_info = pd.DataFrame()
    
    #Get a table that contains each node and the component it belongs.
    components = no_loops_graph.connected_components()
    info_comp_member = components.membership
    original_components = sorted(list(set(info_comp_member)))
    info_comp_size = components.sizes()
    
    node_and_component = pd.DataFrame(data={"Node":[no_loops_graph.vs[index]["name"] for index in range(len(no_loops_graph.vs))],
                                            "Original_component":info_comp_member})
    
    information_components = pd.DataFrame(data={"Original_component":original_components,
                                                "Size":info_comp_size})
    
    full_info_components = node_and_component.merge(information_components, on="Original_component")
    
    #Analyze if each components needs to be sub-divided
    #improve this is always true?? when is len(components_graph) ever 0?
    if(len(components_graph) >= 1):
        #Loop through each component and run the different community detection algorithms.
        #As output we get a table with the different modularity values of each community detection algorithm for each sub-graph.
        for component in range(len(components_graph)):
            subgraph = components_graph[component]
            partition_info = partitioning_components(subgraph)
            first_node = subgraph.vs[0]["name"]
            index = full_info_components.loc[:,"Node"] == first_node
            info_first_node = full_info_components.loc[index,:]
    
            partition_info.loc[:,"Original_component"] = info_first_node.loc[:,"Original_component"][0] #improve this line causes an error (KeyError: 0) when -n is left empty and there aren't yet any gplas directories created before running the gplas command. But only sometimes, most of the time it just works??????
            complete_partition_info = pd.concat([complete_partition_info, partition_info], ignore_index=True)
    
        complete_partition_info.loc[:,"Modularity"] = round(complete_partition_info.loc[:,"Modularity"], 2)
    
        #Get the decision of splitting or not. if the modularity value is bigger than the threshold, split. Otherwise, don't.
        complete_partition_info.loc[:,"Decision"] = ["Split" if score >= modularity_threshold else "No_split" for score in complete_partition_info.loc[:,"Modularity"]]
    
    #Add singletons components, namely Single-node components from the plasmidome network.
    index = [info_comp_size[component] == 1 for component in range(len(info_comp_size))]
    #improve use np.where
    singletons_component = [i for i, x in enumerate(index) if x]
    
    #Add the singletons to the rest of the data.
    if(len(singletons_component) > 0):
        df_singletons = pd.DataFrame(data={"Algorithm":"Independent_single_component",
                                           "Modularity":0,
                                           "Original_component":singletons_component,
                                           "Decision":"No_split"})
        complete_partition_info = pd.concat([complete_partition_info, df_singletons], ignore_index=True)
    
    complete_partition_info.sort_values(by="Original_component", axis=0, ascending=True, inplace=True, ignore_index=True)
    
    #When suitable, get the results from the network partitioning algorithm.
    contigs_membership = []
    #improve can we replace internal_component with just 'component' from the for loop?
    #improve internal_component is iterated with +=1 (only if not singleton??) but it is never reset to 0?
    internal_component = 0
    
    #For determining partition, get the algorithm that provides the biggest modularity value.
    for component in sorted(list(set(complete_partition_info.loc[:,"Original_component"]))):
        index = complete_partition_info.loc[:,"Original_component"] == component
        decision_comp = complete_partition_info.loc[index,:]
    
        if(decision_comp.iloc[0,0] != "Independent_single_component"):
            split_decision = sum(decision_comp.loc[:,"Decision"].str.count("Split"))
            no_split_decision = sum(decision_comp.loc[:,"Decision"].str.count("No_split"))
    
            if(split_decision >= no_split_decision):
                index = decision_comp.loc[:,"Modularity"] == max(decision_comp.loc[:,"Modularity"])
                algorithm_to_split = decision_comp.loc[index,:]
                algorithm = algorithm_to_split.iloc[0,0]
    
                graph_component = components_graph[internal_component]
                spl_names = graph_component.vs["name"]
    
                if(algorithm == "Walktrap"):
                    dendrogram_walktrap = graph_component.community_walktrap()
                    clusters_walktrap = dendrogram_walktrap.as_clustering()
                    spl_membership = clusters_walktrap.membership
    
                elif(algorithm == "Leading-eigen"):
                    clusters_eigen = graph_component.community_leading_eigenvector()
                    spl_membership = clusters_eigen.membership
    
                elif(algorithm == "Louvain"):
                    clusters_louvain = graph_component.community_multilevel()
                    spl_membership = clusters_louvain.membership
    
                """ # improve remove this if we removed prop labels from the earlier function
                elif(algorithm == "Propagating-labels")
                    graph_propag <- cluster_label_prop(components_graph[[internal_component]])
                    spl_membership <- graph_propag$membership
                    spl_names <- graph_propag$names
                """
    
                for node in range(len(spl_names)):
                    contigs_membership.append([algorithm, component, spl_membership[node], spl_names[node]])
    
            else:
                index = full_info_components.loc[:,"Original_component"] == component
                nodes_component = full_info_components.loc[index,:]
    
                for node in range(nodes_component.shape[0]): #TODO is this 0 from membership always 0?
                    contigs_membership.append(["Not_split_component", component, 0, nodes_component.iloc[node,0]])
    
            internal_component += 1
    
        else:
            index = full_info_components.loc[:,"Original_component"] == component
            nodes_component = full_info_components.loc[index,:]
    
            for node in range(nodes_component.shape[0]):
                contigs_membership.append(["Independent_single_component", component, 0, nodes_component.iloc[node,0]])
    
    contigs_membership = pd.DataFrame(data=contigs_membership, columns=["Algorithm","Original_component","Bin","Contig"])
    contigs_membership.loc[:,"Cluster"] = ["-".join([str(contigs_membership.loc[row,"Original_component"]), str(contigs_membership.loc[row,"Bin"])]) for row in range(contigs_membership.shape[0])]
    contigs_membership.sort_values(by="Cluster", axis=0, ascending=True, inplace=True, ignore_index=True)
    contigs_membership.loc[:,"Final_cluster"] = contigs_membership.loc[:,"Cluster"].copy()
    
    unique_clusters = sorted(list(set(contigs_membership.loc[:,"Cluster"])))
    for number in range(len(unique_clusters)):
        cluster_name = unique_clusters[number]
        contigs_membership.loc[:,"Final_cluster"] = [cluster.replace(cluster_name, str(number)) for cluster in contigs_membership.loc[:,"Final_cluster"]]
    
    ### Important: Only for visualisation purposes
    """ #improve remove this
    graph_viz = no_loops_graph.copy()
    #improve do we even need viz at all? cant we just add the plotting data to no_loops_graph?
    """
    
    
    #improve program fails if there are more than 31 clusters and it runs out of hardcoded colors?
    set_colors = ["#add8e6","#d49f36","#507f2d","#84b67c","#a06fda","#df462a","#5a51dc","#5b83db","#c76c2d","#4f49a3","#552095","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975"]
    contigs_membership.loc[:,"Color"] = [set_colors[int(cluster)] for cluster in contigs_membership.loc[:,"Final_cluster"]]
    order_contigs = pd.DataFrame(data=no_loops_graph.vs["name"], columns=["Contig"])
    #improve we do this merge to rearrange contigs_membership. replace with np.match?
    contigs_membership = order_contigs.merge(contigs_membership, on="Contig")
    no_loops_graph.vs["color"] = contigs_membership.loc[:,"Color"]
    
    ig.plot(no_loops_graph,
            target=output_png,
            bbox=(700,700),
            margin=50,
            vertex_size=40,
            vertex_label=no_loops_graph.vs["name"],
            edge_width=no_loops_graph.es["width"],
            edge_color="grey")
    
    #Get the final membership data after partitioning
    results_subgraph = pd.DataFrame(data={"number":contigs_membership.loc[:,"Contig"],
                                          "Component":contigs_membership.loc[:,"Final_cluster"]}) #improve cluster+1 if preferred
    
    #Get all the plasmid nodes
    index = clean_pred.loc[:,"Prob_Plasmid"] >= pred_threshold
    pl_nodes = clean_pred.loc[index,:] # Selecting only contigs predicted as plasmid-derived
    
    """ # improve this is already done in coverage.py?? can we just skip it?
    raw_numbers = [number.split("_", 2)[0] for number in pl_nodes.loc[:,"Contig_name"]]
    pl_nodes.loc[:,"number"] = [number.replace("S", "") for number in raw_numbers]
    """
    
    #Get all the not-assigned nodes (Unbinned and repeats)
    index = ~pl_nodes.loc[:,"number"].isin(results_subgraph.loc[:,"number"])
    pl_notassigned = pl_nodes.loc[index,:]
    #Get the repeats
    index = pl_notassigned.loc[:,"number"].isin(repeats.loc[:,"number"])
    pl_repeats = pl_notassigned.loc[index,:]
    #Get the Unbinned nodes
    index = ~pl_notassigned.loc[:,"number"].isin(repeats.loc[:,"number"])
    pl_unbinned = pl_notassigned.loc[index,:]
    
    #TODO test with data containing isolated nodes
    #If there are isolated nodes, remove them from the unbinned and create a new category
    isolated_nodes = pd.read_csv(path_isolated_nodes, sep="\t", header=0)
    #ASK mistake in R?? pl_isolated is pl_unbinned, filtered with the index of pl_notassigned; probem is copied and present in .py
    index = pl_notassigned.loc[:,"number"].isin(isolated_nodes.loc[:,"number"])
    pl_isolated = pl_unbinned.loc[index,:]
    index = ~pl_unbinned.loc[:,"number"].isin(isolated_nodes.loc[:,"number"])
    pl_unbinned = pl_unbinned.loc[index,:]
    
    #Get the assigned nodes based on the final membership algorithm
    index = pl_nodes.loc[:,"number"].isin(results_subgraph.loc[:,"number"])
    pl_assigned = pl_nodes.loc[index,:]
    
    #Get all the information from the plasmid nodes (Lenght, coverage, bin number, etc)
    full_info_assigned = pl_assigned.merge(results_subgraph, on="number")
    
    #Add unbinned category
    if(pl_unbinned.shape[0] >= 1):
        pl_unbinned.loc[:,"Component"] = "Unbinned"
        full_info_assigned = pd.concat([full_info_assigned, pl_unbinned], ignore_index=True)
    
    #Add isolated nodes category
    if(pl_isolated.shape[0] >= 1):
        isolated_nr = 0
    
        while isolated_nr <= pl_isolated.shape[0]:
            isolated_identification = "_".join(["Isolated", str(isolated_nr)])
            pl_isolated.loc["Component"][isolated_nr] = isolated_identification
            isolated_nr += 1
    
        full_info_assigned = pd.concat([full_info_assigned, pl_isolated], ignore_index=True)
    
    #Add repeat-like category
    if(pl_repeats.shape[0] >= 1):
        pl_repeats.loc[:,"Component"] = "Repeat_like"
        full_info_assigned = pd.concat([full_info_assigned, pl_repeats], ignore_index=True)
    
    #Add information to the results file
    full_info_assigned = full_info_assigned.drop(columns="Contig_length")
    full_info_assigned.loc[:,"Prob_Chromosome"] = round(full_info_assigned.loc[:,"Prob_Chromosome"], 2)
    full_info_assigned.loc[:,"Prob_Plasmid"] = round(full_info_assigned.loc[:,"Prob_Plasmid"], 2)
    full_info_assigned.loc[:,"coverage"] = round(full_info_assigned.loc[:,"coverage"], 2)
    
    #Create the fasta files
    with open(path_nodes) as file:
        raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]
    
    df_nodes = pd.DataFrame(data=raw_nodes, columns=["Contig_name", "Sequence"])
    df_nodes = df_nodes.merge(full_info_assigned, on="Contig_name")
    
    #Write fasta files
    for component in set(df_nodes.loc[:,"Component"]):
        index = df_nodes.loc[:,"Component"] == component
        nodes_component = df_nodes.loc[index,:]
        component_complete_name = "_".join([sample, "bin", str(component)])
        filename = "".join([output_dir, component_complete_name, ".fasta"])
    
        with open(filename, mode="w") as file:
            for contig in range(nodes_component.shape[0]):
                file.write(">" + nodes_component.iloc[contig,0] + "\n" + nodes_component.iloc[contig,1] + "\n")
    
    #improve start the column as 'Bin' instead of changing at the end; like in coocurrence_repeats.py
    full_info_assigned = full_info_assigned.rename(columns={"Component":"Bin"})
    results_subgraph = results_subgraph.rename(columns={"Component":"Bin"})
    
    full_info_assigned.to_csv(output_results, sep='\t', index=False, header=True, mode='w')
    results_subgraph.to_csv(output_components, sep='\t', index=False, header=True, mode='w')
    
    return
#improve change the column order of ecoli_results_no_repeats to match the output of R?
##order node order in 'ecoli_bins_no_repeats' / ecoli_results_no_repeats & co.
### bins_no_repeats is not ordered but results_no_repeats is?
