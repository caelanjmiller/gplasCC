import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser


def scalar1(x):  # TODO can we do something like "from coocurrence.py import scalar1"? it is the same identical function
    denominator = (sum([value*value for value in x]))**0.5
    scaled_x = [value/denominator for value in x]
    return scaled_x


#TODO we have multiple loops with "for row in range(solutions.shape[0]):" can we possibly merge some of them?

def calculate_coocurrence_repeats(sample):
    #Inputs
    path_nodes = f"gplas_input/{sample}_raw_nodes.fasta"
    path_prediction = f"coverage/{sample}_clean_prediction.tab"
    input_solutions = f"walks/repeats/{sample}_solutions.tab"
    path_bins = f"results/{sample}_results_no_repeats.tab"
    clean_repeats_path = f"coverage/{sample}_clean_repeats.tab"

    #Outputs
    output_dir = "results/"
    output_results = f"results/{sample}_results.tab"
    output_components = f"results/{sample}_bins.tab"
    output_chromosomes = f"results/{sample}_chromosome_repeats.tab"

    clean_pred = pd.read_csv(path_prediction, sep='\t', header=0)
    clean_pred = clean_pred.astype({'Prob_Chromosome':float,
                                    'Prob_Plasmid':float,
                                    'Prediction':str,
                                    'Contig_name':str,
                                    'Contig_length':int,
                                    'number':str,
                                    'length':int,
                                    'coverage':float})

    solutions = pd.read_csv(input_solutions, sep='\t', header=None, names=['walks', 'initial_classification', 'unitig_classification', 'path_coverage'])
    max_nodes = max([walk.count(',')+1 for walk in solutions.loc[:,'walks']])
    steps = [''.join(['step_', str(step)]) for step in range(max_nodes)]
    stepsDf = pd.DataFrame(solutions.loc[:,'walks'].str.split(',', expand=True))
    stepsDf.columns = steps
    solutions = pd.concat([solutions, stepsDf], axis=1)

    solutions = solutions.drop(columns='walks')
    col_order = solutions.columns.to_list()
    col_order = col_order[3:] + col_order[:3]
    solutions = solutions[col_order]

    #get the last node of each solution
    last_nodes = []
    for row in range(solutions.shape[0]):
        last_node = [node for node in solutions.iloc[row,0:max_nodes].dropna()][-1]
        last_nodes.append(last_node)

    last_nodes_signless = [node.replace('+','') for node in last_nodes]
    last_nodes_signless = [node.replace('-','') for node in last_nodes_signless]

    solutions.loc[:,'last_nodes'] = last_nodes
    solutions.loc[:,'last_nodes_signless'] = last_nodes_signless

    #====Merge with bin data=============
    #Import information from the bins and match it with total pairs
    bins_data = pd.read_csv(path_bins, sep='\t', header=0)
    bins_data = bins_data.astype({'Prob_Chromosome':float,
                                  'Prob_Plasmid':float,
                                  'Prediction':str,
                                  'Contig_name':str,
                                  'number':str,
                                  'length':int,
                                  'coverage':float,
                                  'Bin':str})
    #merge it
    solutions = solutions.merge(bins_data.loc[:,['number','Bin']], how='left', left_on='last_nodes_signless', right_on='number')
    solutions = solutions.drop(columns='number')
    #assign 'C' to the chromosome and -1 to repeats
    #TODO ASK Julian: repeats are not assigned -1 in the original R-script? (i.e. comment above is lying) does it matter?
    index = solutions.loc[:,'Bin'].isna()
    solutions.loc[index,'Bin'] = 'C'

    #remove connections to unbinned unitigs
    index = solutions.loc[:,'Bin'] != 'Unbinned'
    solutions = solutions.loc[index,:]

    #remove connections to repeats only
    index = solutions.loc[:,'unitig_classification'] != 'Repeat'
    solutions = solutions.loc[index,:]

    #get the initial node without a symbol
    initial_nodes = [node.replace('+','') for node in solutions.loc[:,'step_0']]
    initial_nodes = [node.replace('-','') for node in initial_nodes]
    solutions.loc[:,'initial_node'] = initial_nodes
    solutions.reset_index(inplace=True, drop=True)

    #combine inital nodes with bin in a single variable
    solutions.loc[:,'repeat_bin'] = ['-'.join([solutions.loc[row,'initial_node'], solutions.loc[row,'Bin']]) for row in range(solutions.shape[0])]

    #keep cases in which the last node and second node are the same (repeat directly connected to a unitig)
    solutions.loc[:,'keep'] = ['yes' if ((solutions.loc[row,'step_1'] == solutions.loc[row,'last_nodes']) &
                                         (solutions.loc[row,'unitig_classification'] != 'Repeat')) else 'no' for row in range(solutions.shape[0])]

    index = solutions.loc[:,'keep'] == 'yes'
    final_solutions = solutions.loc[index,:]

    #----Analyze the remainig cases to check if we have the same bin upstream and downstream from the repeat--
    #separate into positive and negative solutions
    index = [solutions.loc[row,'step_1'][-1] == '+' for row in range(solutions.shape[0])]
    positive_solutions = solutions.loc[index,:]
    index = [solutions.loc[row,'step_1'][-1] == '-' for row in range(solutions.shape[0])]
    negative_solutions = solutions.loc[index,:]

    #Keep only the walks if they lead to same bin, or to the chromosome in both directions
    index = ((positive_solutions.loc[:,'repeat_bin'].isin(negative_solutions.loc[:,'repeat_bin'])) &
             (positive_solutions.loc[:,'keep'] == 'no'))
    positive_solutions_keep = positive_solutions.loc[index,:]
    index = ((negative_solutions.loc[:,'repeat_bin'].isin(positive_solutions.loc[:,'repeat_bin'])) &
             (negative_solutions.loc[:,'keep'] == 'no'))
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
    unique_nodes = sorted(list(set(all_nodes)))

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

    starting_nodes = [node for node in unique_nodes if node in solutions.loc[:,'step_0'].values]

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

    total_pairs = pd.DataFrame(total_pairs, columns=['Starting_node', 'Connecting_node', 'weight'])

    total_pairs.loc[:,'Starting_node'] = [node.replace('+','') for node in total_pairs.loc[:,'Starting_node']]
    total_pairs.loc[:,'Starting_node'] = [node.replace('-','') for node in total_pairs.loc[:,'Starting_node']]
    total_pairs.loc[:,'Connecting_node'] = [node.replace('+','') for node in total_pairs.loc[:,'Connecting_node']]
    total_pairs.loc[:,'Connecting_node'] = [node.replace('-','') for node in total_pairs.loc[:,'Connecting_node']]

    #Filter-out cases of no-coocurrence
    #TODO ASK Julian: filter is on weight > 1 why not weight > 0?
    index = [weight > 1 for weight in total_pairs.loc[:,'weight']]
    total_pairs = total_pairs.loc[index,:]

    #Scale weights
    complete_node_info = pd.DataFrame(columns=['Starting_node', 'Connecting_node', 'weight', 'scaled_weight'])
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
        if complete_node_info.shape[0] == 0:  # if statement to prevent "FutureWarning: concatenation with empty DataFrame"
            complete_node_info = particular_node
        else:
            complete_node_info = pd.concat([complete_node_info, particular_node], ignore_index=True)

    total_pairs = complete_node_info

    initial_nodes = [node.replace('+','') for node in starting_nodes]
    initial_nodes = [node.replace('-','') for node in initial_nodes]

    #get a list of clean untigs
    clean_unitigs = clean_pred.loc[:,'number']

    #Filter out connected repeated elements. Keep only connections from starting nodes (repeats) unitigs.
    index = ((total_pairs.loc[:,'Starting_node'].isin(initial_nodes)) &
             (total_pairs.loc[:,'Connecting_node'].isin(clean_unitigs)))
    total_pairs = total_pairs.loc[index,:]

    total_pairs = total_pairs.merge(bins_data.loc[:,['number', 'Bin']], how='left', left_on='Connecting_node', right_on='number')
    total_pairs = total_pairs.drop(columns='number')
    #ASSIGN 'C' to the chromosome
    index = total_pairs.loc[:,'Bin'].isna()
    total_pairs.loc[index,'Bin'] = 'C'

    total_pairs = total_pairs.loc[:,['Starting_node', 'Bin', 'weight', 'scaled_weight']]

    #===Reformat the dataframe to obtain the totality of co-courences
    #First check if we actually have co-ocurrence of unitigs.
    if total_pairs.shape[0] > 0 and total_pairs.shape[1] > 0:
        total_pairs.loc[:,'Pair'] = ['-'.join([total_pairs.loc[row,'Bin'], total_pairs.loc[row,'Starting_node']]) for row in range(total_pairs.shape[0])]
    else:
        return False

    single_edge_counting = []
    for pair in sorted(list(set(total_pairs.loc[:,'Pair']))):
        index = total_pairs.loc[:,'Pair'] == pair
        pairs_subset = total_pairs.loc[index,:]
        sum_weight = sum(pairs_subset.loc[:,'weight'])
        single_edge_counting.append([pair, sum_weight])

    pairs = [pair[0].split('-') for pair in single_edge_counting]

    weight_graph = pd.DataFrame(data={'From_to':[pair[1] for pair in pairs],
                                      'To_from':[pair[0] for pair in pairs],
                                      'weight':[weight[1] for weight in single_edge_counting]})

    #Get the data from coverages. Repeats and bins
    clean_repeats = pd.read_csv(clean_repeats_path, sep='\t', header=0)
    clean_repeats = clean_repeats.astype({'number':str,
                                          'coverage':float})

    bins_coverage = []
    for current_bin in sorted(list(set(bins_data.loc[:,'Bin']))):
        index = bins_data.loc[:,'Bin'] == current_bin
        bins_subset = bins_data.loc[index,:]
        mean_coverage = round(np.mean(bins_subset.loc[:,'coverage']), 2)
        bins_coverage.append([current_bin, mean_coverage])

    bins_coverage = pd.DataFrame(data=bins_coverage, columns=['Bin', 'bin_coverage'])

    weight_graph = weight_graph.merge(clean_repeats.loc[:,['number', 'coverage']], how='left', left_on='From_to', right_on='number')
    #TODO find a more elegant/efficient way to properly merge dataframes
    weight_graph = weight_graph.drop(columns='number')

    weight_graph = weight_graph.merge(bins_coverage, how='left', left_on='To_from', right_on='Bin')
    weight_graph = weight_graph.drop(columns='Bin')

    #assign a coverage of 1 to chromosomal unitigs
    index = weight_graph.loc[:,'bin_coverage'].isna()
    weight_graph.loc[index,'bin_coverage'] = float(1)

    #Explore if the combination of bins proposed by the algorithm is plausible based on coverage
    repeat_assignments = pd.DataFrame(columns=['From_to', 'To_from', 'weight', 'coverage', 'bin_coverage', 'rank'])
    #loop through each of the repeats
    for node in sorted(list(set(weight_graph.loc[:,'From_to']))):
        index = weight_graph.loc[:,'From_to'] == node
        df_node = weight_graph.loc[index,:].copy()
        #create a rank of the most likely connections, based on the co-ocurence count
        df_node.sort_values(by='weight', axis=0, ascending=False, ignore_index=True, inplace=True)
        df_node.loc[:,'rank'] = df_node.index.tolist()

        rank = 0 #start from the highest ranking (0)
        accumulated_cov = 0
        while rank <= max(df_node.loc[:,'rank']):
            index = df_node.loc[:,'rank'] == rank
            repeat_bin = df_node.loc[index,:]
            if repeat_bin.loc[:,'coverage'].values[0] + accumulated_cov >= repeat_bin.loc[:,'bin_coverage'].values[0]:
                accumulated_cov += repeat_bin.loc[:,'bin_coverage'].values[0]
                if repeat_assignments.shape[0] == 0: # if statement to prevent "FutureWarning: concatenation with empty DataFrame"
                    repeat_assignments = repeat_bin
                else:
                    repeat_assignments = pd.concat([repeat_assignments, repeat_bin], ignore_index=True)
            rank += 1
            
    repeat_assignments = repeat_assignments.loc[:,['To_from','From_to']].rename(columns={'To_from':'Bin',
                                                                                         'From_to':'number'})
    #separate results into plasmid and chromosome repeats
    index = repeat_assignments.loc[:,'Bin'] == 'C'
    chromosome_repeats = repeat_assignments.loc[index,:]

    index = repeat_assignments.loc[:,'Bin'] != 'C'
    plasmid_repeats = repeat_assignments.loc[index,:]

    print("Adding repeated elements to the predictions........... completed!")
    if plasmid_repeats.shape[0] == 0:
        print("gplas did not find repeated elements associated with plasmid predictions")
        bins_data.loc[:,'Prob_Chromosome'] = round(bins_data.loc[:,'Prob_Chromosome'], 2)
        bins_data.loc[:,'Prob_Plasmid'] = round(bins_data.loc[:,'Prob_Plasmid'], 2)
        bins_data.loc[:,'coverage'] = round(bins_data.loc[:,'coverage'], 2)
        full_info_assigned = bins_data
    else:
        print("gplas found repeated elements associated with plasmid predictions")
        #Get all the repeat nodes
        index = clean_repeats.loc[:,'number'].isin(plasmid_repeats.loc[:,'number']) # Selecting only contigs predicted as plasmid-derived
        pl_nodes = clean_repeats.loc[index,:]

        #Get all the information from the plasmid nodes (Length, coverage, bin number, etc)
        full_info_assigned = pl_nodes.merge(plasmid_repeats, on='number')

        #Add information to the results file
        full_info_assigned = full_info_assigned.drop(columns='Contig_length')
        full_info_assigned.loc[:,'Prob_Chromosome'] = round(full_info_assigned.loc[:,'Prob_Chromosome'], 2)
        full_info_assigned.loc[:,'Prob_Plasmid'] = round(full_info_assigned.loc[:,'Prob_Plasmid'], 2)
        full_info_assigned.loc[:,'coverage'] = round(full_info_assigned.loc[:,'coverage'], 2)
        full_info_assigned.loc[:,'Prediction'] = 'Repeat'

        bins_data.loc[:,'Prob_Chromosome'] = round(bins_data.loc[:,'Prob_Chromosome'], 2)
        bins_data.loc[:,'Prob_Plasmid'] = round(bins_data.loc[:,'Prob_Plasmid'], 2)
        bins_data.loc[:,'coverage'] = round(bins_data.loc[:,'coverage'], 2)

        full_info_assigned = pd.concat([full_info_assigned, bins_data], ignore_index=True)

    #Create the fasta files
    with open(path_nodes) as file:
        raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]

    df_nodes = pd.DataFrame(data=raw_nodes, columns=['Contig_name', 'Sequence'])
    df_nodes = df_nodes.merge(full_info_assigned, on='Contig_name')

    #Write fasta files
    for component in sorted(list(set(df_nodes.loc[:,'Bin']))):
        index = df_nodes.loc[:,'Bin'] == component
        nodes_component = df_nodes.loc[index,:]
        component_complete_name = f"{sample}_bin_{component}"
        filename = f"{output_dir}{component_complete_name}.fasta"

        with open(filename, mode='w') as file:
            for contig in range(nodes_component.shape[0]):
                file.write('>' + nodes_component.iloc[contig,0] + '\n' + nodes_component.iloc[contig,1] + '\n')

    results_summary = df_nodes.loc[:,['number', 'Bin']]

    full_info_assigned.to_csv(output_results, sep='\t', index=False, header=True, mode='w')
    results_summary.to_csv(output_components, sep='\t', index=False, header=True, mode='w')

    #format chromosome repeats and print
    chromosome_repeats = chromosome_repeats.loc[:,['number', 'Bin']] #TODO if possible change column order from the start instead of changing it at the end
    chromosome_repeats.loc[:,'Bin'] = 'Chromosome' #TODO if possible call the column name 'Chromosome' from the start instead of 'C' and renaming it at the end
    chromosome_repeats.to_csv(output_chromosomes, sep='\t', index=False, header=True, mode='w')

    return True
