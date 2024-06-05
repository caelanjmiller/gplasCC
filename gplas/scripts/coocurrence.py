import pandas as pd
import igraph as ig
from Bio.SeqIO.FastaIO import SimpleFastaParser
import warnings


def scalar1(x):
    denominator = (sum([value*value for value in x]))**0.5
    scaled_x = [value/denominator for value in x]
    return scaled_x


#Create function to run three community detection algorithms on the plasmidome network.
#This fuction creates a dataframe as an output that contains the different algorithms names and the global modularity values of the resulting network.
def partitioning_components(graph):
    # Walktrap algorithm
    dendrogram_walktrap = graph.community_walktrap()
    clusters_walktrap = dendrogram_walktrap.as_clustering()
    modularity_walktrap = clusters_walktrap.modularity

    # Leading eigen values
    # TODO find a better way to fix non-converging eigenvalue error
    # For some samples there is a ~10% chance for the eigenvalue to not converge and raise an error message
    # We use this code to catch this error and/or a runtime warning and simply try the algorithm again up to 5 times
    # It is not a pretty fix but it seems to work for now
    warnings.filterwarnings("error") # allow warnings to be caught as exceptions
    for attempt in range(5):
        try:
            clusters_eigen = graph.community_leading_eigenvector()
            break
        except (ig._igraph.InternalError, RuntimeWarning):
            pass
    warnings.resetwarnings() # return warnings to normal behavior
    modularity_eigen = clusters_eigen.modularity

    # Louvain method
    clusters_louvain = graph.community_multilevel()
    modularity_louvain = clusters_louvain.modularity

    partition_info = pd.DataFrame(data={'Algorithm':['Walktrap', 'Leading-eigen', 'Louvain'],
                                        'Modularity':[modularity_walktrap, modularity_eigen, modularity_louvain]})
    return partition_info


#TODO we have multiple loops with "for row in range(solutions.shape[0]):" can we possibly merge some of them?
def calculate_coocurrence(sample, number_iterations, pred_threshold, modularity_threshold, mode='normal'):
    if mode == 'normal':
        subdir = 'normal_mode/'
    elif mode == 'unbinned':
        subdir = ''
    #Inputs
    path_nodes = f"gplas_input/{sample}_raw_nodes.fasta"
    path_prediction = f"coverage/{sample}_clean_prediction.tab"
    path_graph_repeats = f"coverage/{sample}_repeats_graph.tab"
    path_isolated_nodes = f"coverage/{sample}_isolated_nodes.tab"
    input_solutions = f"walks/{subdir}{sample}_solutions.tab"

    #Outputs
    output_dir = f"results/{subdir}"
    output_results = f"results/{subdir}{sample}_results_no_repeats.tab"
    output_components = f"results/{subdir}{sample}_bins_no_repeats.tab"
    output_png = f"results/{subdir}{sample}_plasmidome_network.png"

    clean_pred = pd.read_csv(path_prediction, sep='\t', header=0)
    clean_pred = clean_pred.astype({'Prob_Chromosome':float,
                                    'Prob_Plasmid':float,
                                    'Prediction':str,
                                    'Contig_name':str,
                                    'Contig_length':int,
                                    'number':str,
                                    'length':int,
                                    'coverage':float})

    repeats = pd.read_csv(path_graph_repeats, sep='\t', header=0)
    repeats.loc[:,'number'] = [name.replace('+','') for name in repeats['number']]
    repeats.loc[:,'number'] = [name.replace('-','') for name in repeats['number']]

    with open(input_solutions) as file:
        max_nodes = max([line.count('\t')+1 for line in file])

    solutions = pd.read_csv(input_solutions, sep='\t', header=None, names=range(max_nodes), engine='python')
    #TODO if the line above does not define engine='python', it produces a pandas ParserError. but only for 2 of the 214 ecoli test samples?????
    #this is not a problem anywhere else in the gplas scripts or for any other test samples
    #GCA_013602835.1_ASM1360283v1
    #GCA_013802065.1_ASM1380206v1
    #ParserError: Error tokenizing data. C error: Buffer overflow caught - possible malformed input file.
    
    #create a list with all the nodes that appear in the plasmid walks.
    all_nodes = []
    for row in range(solutions.shape[0]):
        nodes = [node for node in solutions.loc[row,:].dropna()]
        all_nodes.extend(nodes)

    #get unique set of nodes
    unique_nodes = sorted(list(set(all_nodes)))

    #TODO turn co-ocurrence matrix creation into a seperate function? moving code to functions might be useful for reusing/importing code in coocurrence_repeats
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
            if connecting_node != node:
                weight = col_sums[index_connecting_node]
                total_pairs.append([node, connecting_node, weight])

    #Find circular sequences
    circular_sequences = []
    #Extract the walks in which start-node and end-node are the same.
    for row in range(solutions.shape[0]):
        walk = [node for node in solutions.loc[row,:].dropna()]
        if len(walk) > 1 and walk[0] == walk[-1]:
            circular_sequences.append([walk[0], walk[-1]])

    #check if the number of circular walks starting from each node equals the number of iterations.
    #if this is the case, add the circular walk to total_pairs
    if len(circular_sequences) > 0:
        no_duplicated = [sorted(list(unique_walk)) for unique_walk in set(tuple(walk) for walk in circular_sequences)] # remove duplicate entries

        for combination in range(len(no_duplicated)):
            combi = no_duplicated[combination]
            total_ocurrences = sum([walk[1] == combi[1] for walk in circular_sequences])
            if total_ocurrences == number_iterations:
                total_pairs.append([combi[0], combi[1], total_ocurrences])

    total_pairs = pd.DataFrame(total_pairs, columns=['Starting_node', 'Connecting_node', 'weight'])

    total_pairs.loc[:,'Starting_node'] = [node.replace('+','') for node in total_pairs.loc[:,'Starting_node']]
    total_pairs.loc[:,'Starting_node'] = [node.replace('-','') for node in total_pairs.loc[:,'Starting_node']]
    total_pairs.loc[:,'Connecting_node'] = [node.replace('+','') for node in total_pairs.loc[:,'Connecting_node']]
    total_pairs.loc[:,'Connecting_node'] = [node.replace('-','') for node in total_pairs.loc[:,'Connecting_node']]

    #Filter-out cases of no-coocurrence
    #TODO ASK Julian: we filter on weight > 1 why not weight > 0?
    index = [weight > 1 for weight in total_pairs.loc[:,'weight']]
    total_pairs = total_pairs.loc[index,:]

    #Scale weights
    complete_node_info = pd.DataFrame(columns=['Starting_node', 'Connecting_node', 'weight'])
    for node in sorted(list(set(total_pairs.loc[:,'Starting_node']))):
        index = total_pairs.loc[:,'Starting_node'] == node
        first_node = total_pairs.loc[index,:]
        particular_node = []

        for connecting_node in sorted(list(set(first_node.loc[:,'Connecting_node']))):
            index = first_node.loc[:,'Connecting_node'] == connecting_node
            first_second_nodes = first_node.loc[index,:]
            total_weight = sum(first_second_nodes.loc[:,'weight'])
            particular_node.append([node, connecting_node, total_weight])

        particular_node = pd.DataFrame(particular_node, columns=['Starting_node', 'Connecting_node', 'weight'])
        particular_node.loc[:,'scaled_weight'] = scalar1(particular_node.loc[:,'weight']) # add a column with scaled weights using scalar1()
        complete_node_info = pd.concat([complete_node_info, particular_node], ignore_index=True)

    total_pairs = complete_node_info

    initial_nodes = [node.replace('+','') for node in starting_nodes]
    initial_nodes = [node.replace('-','') for node in initial_nodes]

    #Filter out repeated elements. Keep only connections of unitigs.
    index = ((total_pairs.loc[:,'Starting_node'].isin(initial_nodes)) &
             (total_pairs.loc[:,'Connecting_node'].isin(initial_nodes)))
    total_pairs = total_pairs.loc[index,:]

    #Reformat the dataframe to obtain the totality of co-courences
    weight_counting = []
    #First check if we actually have co-ocurrence of unitigs.
    if total_pairs.shape[0] > 0 and total_pairs.shape[1] > 0:
        for row in range(total_pairs.shape[0]):
            initial_node = total_pairs.iloc[row,0]
            connecting_node = total_pairs.iloc[row,1]
            raw_count = int(total_pairs.iloc[row,2])
            if int(initial_node) < int(connecting_node):
                pair = f"{initial_node}-{connecting_node}"
            else:
                pair = f"{connecting_node}-{initial_node}"

            weight_counting.append([pair, raw_count])
    else:
        index = clean_pred.loc[:,'Prob_Plasmid'] > pred_threshold
        pl_unbinned = clean_pred.loc[index,:].copy()
        pl_unbinned = pl_unbinned.drop(columns='Contig_length')
        pl_unbinned['Component'] = pd.Series(dtype='object')
        pl_unbinned['Prob_Chromosome'] = pd.Series(dtype='float')
        pl_unbinned['Prob_Plasmid'] = pd.Series(dtype='float')
        pl_unbinned['coverage'] = pd.Series(dtype='float')
        if pl_unbinned.shape[0] > 0:
            pl_unbinned.loc[:,'Component'] = 'Unbinned'
            pl_unbinned.loc[:,'Prob_Chromosome'] = round(pl_unbinned.loc[:,'Prob_Chromosome'], 2)
            pl_unbinned.loc[:,'Prob_Plasmid'] = round(pl_unbinned.loc[:,'Prob_Plasmid'], 2)
            pl_unbinned.loc[:,'coverage'] = round(pl_unbinned.loc[:,'coverage'], 2)

        with open(path_nodes) as file:
            raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]

        df_nodes = pd.DataFrame(data=raw_nodes, columns=['Contig_name', 'Sequence'])
        df_nodes = df_nodes.merge(pl_unbinned, on='Contig_name')

        for component in sorted(list(set(df_nodes.loc[:,'Component']))):
            index = df_nodes.loc[:,'Component'] == component
            nodes_component = df_nodes.loc[index,:]
            component_complete_name = f"{sample}_bin_{component}"
            filename = f"{output_dir}{component_complete_name}.fasta"

            with open(filename, mode='w') as file:
                for contig in range(nodes_component.shape[0]):
                    file.write('>' + nodes_component.iloc[contig,0] + '\n' + nodes_component.iloc[contig,1] + '\n')
        #TODO if possible, name the column 'Bin' from the start instead of 'Component' and then renaming it later
        pl_unbinned = pl_unbinned.rename(columns={'Component':'Bin'})
        results_subgraph = pl_unbinned.loc[:,['number','Bin']]

        pl_unbinned.to_csv(output_results, sep='\t', index=False, header=True, mode='w')
        results_subgraph.to_csv(output_components, sep='\t', index=False, header=True, mode='w')
        return False

    weight_counting = pd.DataFrame(weight_counting, columns=['Pair','Count'])

    unique_pairs = sorted(list(set(weight_counting.loc[:,'Pair'])))
    sum_weights = []
    for pair in unique_pairs:
        index = list(weight_counting.loc[:,'Pair'] == pair)
        sum_weights.append(sum(weight_counting.loc[index,'Count']))

    pairs = [pair.split('-') for pair in unique_pairs]

    weight_graph = pd.DataFrame(data={'From_to':[pair[0] for pair in pairs],
                                      'To_from':[pair[1] for pair in pairs],
                                      'weight':sum_weights})

    full_graph_info = pd.DataFrame()

    for node in sorted(list(set(weight_graph.loc[:,'From_to']))):
        index = weight_graph.loc[:,'From_to'] == node
        df_node = weight_graph.loc[index,:].copy()
        df_node.loc[:,'scaled_weight'] = scalar1(df_node.loc[:,'weight'])
        full_graph_info = pd.concat([full_graph_info, df_node], ignore_index=True)

    full_graph_info.loc[:,'width'] = full_graph_info.loc[:,'scaled_weight'].copy() * 5

    graph_pairs = ig.Graph.DataFrame(full_graph_info, directed=False, use_vids=False)

    # Simplifying the graph - Removing self-loops from the graph
    no_loops_graph = graph_pairs.simplify(multiple=False)

    #Analyze if the plasmidome network can be partitioned into different sub-networks.
    #Get clusters based on connectivity alone.
    components_graph = no_loops_graph.decompose(mode='weak', minelements=2)
    
    #Get a table that contains each node and the component it belongs to
    components = no_loops_graph.connected_components()
    info_comp_member = components.membership
    original_components = sorted(list(set(info_comp_member)))
    info_comp_size = components.sizes()

    node_and_component = pd.DataFrame(data={'Node':[no_loops_graph.vs[index]['name'] for index in range(len(no_loops_graph.vs))],
                                            'Original_component':info_comp_member})

    information_components = pd.DataFrame(data={'Original_component':original_components,
                                                'Size':info_comp_size})

    full_info_components = node_and_component.merge(information_components, on='Original_component')

    #Analyze if each component needs to be sub-divided
    #Create a data-frame to hold the results from the different clusters
    complete_partition_info = pd.DataFrame()
    #TODO this statement is always true?? when is len(components_graph) ever 0? can we remove this statement?
    if len(components_graph) >= 1:
        #Loop through each component and run the different community detection algorithms.
        #As output we get a table with the different modularity values of each community detection algorithm for each sub-graph.
        for component in range(len(components_graph)):
            subgraph = components_graph[component]
            partition_info = partitioning_components(subgraph)
            first_node = subgraph.vs[0]['name']
            index = full_info_components.loc[:,'Node'] == first_node
            info_first_node = full_info_components.loc[index,:]

            partition_info.loc[:,'Original_component'] = info_first_node.loc[:,'Original_component'].values[0]
            complete_partition_info = pd.concat([complete_partition_info, partition_info], ignore_index=True)

        complete_partition_info.loc[:,'Modularity'] = round(complete_partition_info.loc[:,'Modularity'], 2)

        #Get the decision of splitting or not. if the modularity value is bigger than the threshold, split. Otherwise, don't.
        complete_partition_info.loc[:,'Decision'] = ['Split' if score >= modularity_threshold else 'No_split' for score in complete_partition_info.loc[:,'Modularity']]

    #Add singletons components, namely Single-node components from the plasmidome network.
    index = [info_comp_size[component] == 1 for component in range(len(info_comp_size))]
    singletons_component = [i for i, x in enumerate(index) if x]  # TODO replace this list comprehension with np.where?

    #Add the singletons to the rest of the data.
    if len(singletons_component) > 0:
        df_singletons = pd.DataFrame(data={'Algorithm':'Independent_single_component',
                                           'Modularity':0,
                                           'Original_component':singletons_component,
                                           'Decision':'No_split'})
        complete_partition_info = pd.concat([complete_partition_info, df_singletons], ignore_index=True)

    complete_partition_info.sort_values(by='Original_component', axis=0, ascending=True, inplace=True, ignore_index=True)

    #When suitable, get the results from the network partitioning algorithm.
    contigs_membership = []

    internal_component = 0

    #For determining partition, get the algorithm that provides the biggest modularity value.
    for component in sorted(list(set(complete_partition_info.loc[:,'Original_component']))):
        index = complete_partition_info.loc[:,'Original_component'] == component
        decision_comp = complete_partition_info.loc[index,:]

        if decision_comp.iloc[0,0] != 'Independent_single_component':
            split_decision = sum(decision_comp.loc[:,'Decision'].str.count('Split'))
            no_split_decision = sum(decision_comp.loc[:,'Decision'].str.count('No_split'))

            if split_decision >= no_split_decision:
                index = decision_comp.loc[:,'Modularity'] == max(decision_comp.loc[:,'Modularity'])
                algorithm_to_split = decision_comp.loc[index,:]
                algorithm = algorithm_to_split.iloc[0,0]

                graph_component = components_graph[internal_component]
                spl_names = graph_component.vs['name']

                if algorithm == 'Walktrap':
                    dendrogram_walktrap = graph_component.community_walktrap()
                    clusters_walktrap = dendrogram_walktrap.as_clustering()
                    spl_membership = clusters_walktrap.membership

                elif algorithm == 'Leading-eigen':
                    clusters_eigen = graph_component.community_leading_eigenvector()
                    spl_membership = clusters_eigen.membership

                elif algorithm == 'Louvain':
                    clusters_louvain = graph_component.community_multilevel()
                    spl_membership = clusters_louvain.membership

                for node in range(len(spl_names)):
                    contigs_membership.append([algorithm, component, spl_membership[node], spl_names[node]])

            else:
                index = full_info_components.loc[:,'Original_component'] == component
                nodes_component = full_info_components.loc[index,:]

                for node in range(nodes_component.shape[0]): #TODO check if this 0 from membership is always 0?
                    contigs_membership.append(['Not_split_component', component, 0, nodes_component.iloc[node,0]])

            internal_component += 1

        else:
            index = full_info_components.loc[:,'Original_component'] == component
            nodes_component = full_info_components.loc[index,:]

            for node in range(nodes_component.shape[0]):
                contigs_membership.append(['Independent_single_component', component, 0, nodes_component.iloc[node,0]])

    contigs_membership = pd.DataFrame(data=contigs_membership, columns=['Algorithm','Original_component','Bin','Contig'])
    contigs_membership.loc[:,'Cluster'] = ['-'.join([str(contigs_membership.loc[row,'Original_component']), str(contigs_membership.loc[row,'Bin'])]) for row in range(contigs_membership.shape[0])]
    contigs_membership.sort_values(by='Cluster', axis=0, ascending=True, inplace=True, ignore_index=True)
    contigs_membership.loc[:,'Final_cluster'] = contigs_membership.loc[:,'Cluster'].copy()

    unique_clusters = sorted(list(set(contigs_membership.loc[:,'Cluster'])))
    for number in range(len(unique_clusters)):
        cluster_name = unique_clusters[number]
        contigs_membership.loc[:,'Final_cluster'] = [cluster.replace(cluster_name, str(number)) for cluster in contigs_membership.loc[:,'Final_cluster']]

    #TODO what happens if there are more than 31 clusters and it runs out of hardcoded colors?
    set_colors = ['#add8e6','#d49f36','#507f2d','#84b67c','#a06fda','#df462a','#5a51dc','#5b83db','#c76c2d','#4f49a3','#552095','#82702d','#dd6bbb','#334c22','#d83979','#55baad','#dc4555','#62aad3','#8c3025','#417d61','#862977','#bba672','#403367','#da8a6d','#a79cd4','#71482c','#c689d0','#6b2940','#d593a7','#895c8b','#bd5975']
    contigs_membership.loc[:,'Color'] = [set_colors[int(cluster)] for cluster in contigs_membership.loc[:,'Final_cluster']]
    order_contigs = pd.DataFrame(data=no_loops_graph.vs['name'], columns=['Contig'])
    contigs_membership = order_contigs.merge(contigs_membership, on='Contig')
    no_loops_graph.vs['color'] = contigs_membership.loc[:,'Color']

    ig.plot(no_loops_graph,
            target=output_png,
            bbox=(700,700),
            margin=50,
            vertex_size=40,
            vertex_label=no_loops_graph.vs['name'],
            edge_width=no_loops_graph.es['width'],
            edge_color='grey')

    #Get the final membership data after partitioning
    results_subgraph = pd.DataFrame(data={'number':contigs_membership.loc[:,'Contig'],
                                          'Component':contigs_membership.loc[:,'Final_cluster']})

    #Get all the plasmid nodes
    index = clean_pred.loc[:,'Prob_Plasmid'] >= pred_threshold
    pl_nodes = clean_pred.loc[index,:] # Selecting only contigs predicted as plasmid-derived

    #Get all the not-assigned nodes (Unbinned and repeats)
    index = ~pl_nodes.loc[:,'number'].isin(results_subgraph.loc[:,'number'])
    pl_notassigned = pl_nodes.loc[index,:]
    #Get the repeats
    index = pl_notassigned.loc[:,'number'].isin(repeats.loc[:,'number'])
    pl_repeats = pl_notassigned.loc[index,:]
    #Get the Unbinned nodes
    index = ~pl_notassigned.loc[:,'number'].isin(repeats.loc[:,'number'])
    pl_unbinned = pl_notassigned.loc[index,:]

    #If there are isolated nodes, remove them from the unbinned and create a new category
    isolated_nodes = pd.read_csv(path_isolated_nodes, sep='\t', header=0)
    #TODO ASK Julian: this looks like a mistake in original R-script?
    #pl_isolated is pl_unbinned, filtered with the index of pl_notassigned?
    #this mistake is copied from the R-script and still present in this python script
    index = pl_notassigned.loc[:,'number'].isin(isolated_nodes.loc[:,'number'])
    pl_isolated = pl_unbinned.loc[index,:]
    index = ~pl_unbinned.loc[:,'number'].isin(isolated_nodes.loc[:,'number'])
    pl_unbinned = pl_unbinned.loc[index,:]

    #Get the assigned nodes based on the final membership algorithm
    index = pl_nodes.loc[:,'number'].isin(results_subgraph.loc[:,'number'])
    pl_assigned = pl_nodes.loc[index,:]

    #Get all the information from the plasmid nodes (Lenght, coverage, bin number, etc)
    full_info_assigned = pl_assigned.merge(results_subgraph, on='number')

    #Add unbinned category
    if pl_unbinned.shape[0] >= 1:
        pl_unbinned.loc[:,'Component'] = 'Unbinned'
        full_info_assigned = pd.concat([full_info_assigned, pl_unbinned], ignore_index=True)

    #Add isolated nodes category
    if pl_isolated.shape[0] >= 1:
        isolated_nr = 0

        while isolated_nr <= pl_isolated.shape[0]:
            isolated_identification = f"Isolated_{isolated_nr}"
            pl_isolated.loc['Component'][isolated_nr] = isolated_identification
            isolated_nr += 1

        full_info_assigned = pd.concat([full_info_assigned, pl_isolated], ignore_index=True)

    #Add repeat-like category
    if pl_repeats.shape[0] >= 1:
        pl_repeats.loc[:,'Component'] = 'Repeat_like'
        full_info_assigned = pd.concat([full_info_assigned, pl_repeats], ignore_index=True)

    #Add information to the results file
    full_info_assigned = full_info_assigned.drop(columns='Contig_length')
    full_info_assigned.loc[:,'Prob_Chromosome'] = round(full_info_assigned.loc[:,'Prob_Chromosome'], 2)
    full_info_assigned.loc[:,'Prob_Plasmid'] = round(full_info_assigned.loc[:,'Prob_Plasmid'], 2)
    full_info_assigned.loc[:,'coverage'] = round(full_info_assigned.loc[:,'coverage'], 2)

    #Create the fasta files
    with open(path_nodes) as file:
        raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]

    df_nodes = pd.DataFrame(data=raw_nodes, columns=['Contig_name', 'Sequence'])
    df_nodes = df_nodes.merge(full_info_assigned, on='Contig_name')

    #Write fasta files
    for component in sorted(list(set(df_nodes.loc[:,'Component']))):
        index = df_nodes.loc[:,'Component'] == component
        nodes_component = df_nodes.loc[index,:]
        component_complete_name = f"{sample}_bin_{component}"
        filename = f"{output_dir}{component_complete_name}.fasta"

        with open(filename, mode='w') as file:
            for contig in range(nodes_component.shape[0]):
                file.write('>' + nodes_component.iloc[contig,0] + '\n' + nodes_component.iloc[contig,1] + '\n')

    #TODO if possible start the column name as 'Bin' instead of renaming it at the end; like in coocurrence_repeats.py
    full_info_assigned = full_info_assigned.rename(columns={'Component':'Bin'})
    results_subgraph = results_subgraph.rename(columns={'Component':'Bin'})

    full_info_assigned.to_csv(output_results, sep='\t', index=False, header=True, mode='w')
    results_subgraph.to_csv(output_components, sep='\t', index=False, header=True, mode='w')

    return True
