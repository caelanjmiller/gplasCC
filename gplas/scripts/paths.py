import pandas as pd
import numpy as np
import scipy.stats


def generate_paths(sample, number_iterations, filtering_threshold, sd_coverage=1, mode='normal'):
    #Inputs
    path_links = f"coverage/{sample}_clean_links.tab"
    path_prediction = f"coverage/{sample}_clean_prediction.tab"
    path_graph_contigs = f"coverage/{sample}_graph_contigs.tab"
    path_graph_repeats = f"coverage/{sample}_repeats_graph.tab"
    path_init_nodes = f"coverage/{sample}_initialize_nodes.tab"
    path_cov_variation = f"coverage/{sample}_estimation.txt"
    #Params
    # TODO make these into user-tunable parameters?
    number_nodes = 100
    prob_small_repeats = 0.5
    #Outputs Normal Mode
    if mode == 'normal':
        output_path = f"walks/normal_mode/{sample}_solutions.tab"
    #Outputs Bold Mode
    elif mode == 'bold':
        output_path = f"walks/bold_mode/{sample}_solutions_bold.tab"

    links = pd.read_csv(path_links, sep='\t', header=None)
    clean_pred = pd.read_csv(path_prediction, sep='\t', header=0)
    clean_pred = clean_pred.astype({'Prob_Chromosome':float,
                                    'Prob_Plasmid':float,
                                    'Prediction':str,
                                    'Contig_name':str,
                                    'Contig_length':int,
                                    'number':str,
                                    'length':int,
                                    'coverage':float})

    graph_contigs = pd.read_csv(path_graph_contigs, sep='\t', header=0)

    small_contigs = graph_contigs[graph_contigs['length'] < 500].copy()
    small_contigs_signed_nodes = small_contigs.copy()
    small_contigs_signed_nodes = [node for node in small_contigs_signed_nodes['number']]
    small_contigs.loc[:,'number'] = [name.replace('+','') for name in small_contigs['number']]
    small_contigs.loc[:,'number'] = [name.replace('-','') for name in small_contigs['number']]

    repeats = pd.read_csv(path_graph_repeats, sep='\t', header=0)
    repeats_signed_nodes = repeats.copy()
    repeats_signed_nodes = [node for node in repeats_signed_nodes['number']]
    repeats.loc[:,'number'] = [name.replace('+','') for name in repeats['number']]
    repeats.loc[:,'number'] = [name.replace('-','') for name in repeats['number']]

    initialize_nodes = pd.read_csv(path_init_nodes, sep='\t', header=None)
    initialize_nodes = [str(node) for node in initialize_nodes.iloc[:,0]]
    
    with open(path_cov_variation, mode='r') as file:
        max_variation = float(file.readline()) * sd_coverage


    def plasmid_graph(initial_seed, output_path, links, number_iterations, number_nodes, max_variation, filtering_threshold, prob_small_repeats, direction):
        paths_list = []
        for iteration in range(number_iterations): # Number of times we repeat this process
            path = [initial_seed] # We add the initial seed to the path, first element in the list
            seed = initial_seed

            #################################### Coverage of the current path #################################

            # Extracting the info from our current path                     
            info_path = graph_contigs[[number in path for number in graph_contigs['number']]].copy()
            info_path = info_path[[number not in repeats_signed_nodes for number in info_path['number']]] # Removing the contigs corresponding to transposases
            length_path = sum(info_path['length']) # Length of the path
            info_path['contribution'] = info_path['length']/length_path # Shorter contigs should have a lower contribution to the coverage. k-mer coverage on these contigs fluctuates drastically
            path_mean = np.average(info_path['coverage'], weights=info_path['contribution']) # Coverage of the current path

            ##################### Elongating the path ###########################################

            for elongation in range(number_nodes):
                index = links.loc[:,0] == seed
                current_links = links.loc[index,:] # Consider the last element present in our path and observe all the possible links

                if(current_links.shape[0] == 0):
                    output = '\t'.join(path)
                    paths_list.append(output)
                    path = [initial_seed] # There are no connections possible from this contig
                    break # Exiting elongation loop

                list_connections = sorted(list(set(current_links[2]))) # All the possible unique connections 

                # We do not allow that a node which is not a repeat appears more than 1 time in any solution but we exclude the initial seed from this consideration 
                if(len(path) > 1): # If the path has more than one element        
                    remove_nodes = path[1:]

                    remove_nodes = [node.replace('+','') for node in remove_nodes]
                    remove_nodes = [node.replace('-','') for node in remove_nodes]

                    positive_remove_nodes = [node + '+' for node in remove_nodes]
                    negative_remove_nodes = [node + '-' for node in remove_nodes]

                    ommit_nodes = positive_remove_nodes + negative_remove_nodes

                    first_node = path[0]

                    if(direction == 'forward'):
                        first_node_to_exclude = first_node.replace('+','-')

                    elif(direction == 'reverse'):
                        first_node_to_exclude = first_node.replace('-','+')

                    ommit_nodes.append(first_node_to_exclude)

                    # We need to remove directionality from the path to avoid paths e.g. 54-,161-,54+
                    ommit_nodes = [node for node in ommit_nodes if node not in repeats_signed_nodes]

                    list_connections = [str(connection) for connection in list_connections if connection not in ommit_nodes] # Avoiding loops within the solution

                total_connections = len(list_connections) # Number of connections
                base_probabilities = [float(0)]*total_connections # Generating a list with the connection probabilities

                if(total_connections < 1):
                    output = '\t'.join(path)
                    paths_list.append(output)
                    path = [initial_seed] # There are no connections possible from this contig
                    break # Exiting elongation loop

                ################################# Probabilities of the connections #####################################################

                # We generate a dataframe with the connection probabilities
                prob_df = pd.DataFrame(data={'number':list_connections,
                                             'Prob_Plasmid':base_probabilities})

                # Replacing the base probabilities with the probabilities generated by mlplasmids 
                prob_df.loc[:,'number'] = [number.replace('+','') for number in prob_df['number']]
                prob_df.loc[:,'number'] = [number.replace('-','') for number in prob_df['number']]

                for number in prob_df.loc[:,'number']:
                    matchID = clean_pred.index[clean_pred['number'] == number]
                    if any(matchID):
                        index = prob_df.loc[:,'number'] == number
                        prob_df.loc[index, 'Prob_Plasmid'] = float(clean_pred.loc[matchID, 'Prob_Plasmid'].values[0])

                # Short contigs do not have a reliable probability, we assign them a predefined probability (passed via the argument 'prob_small_repeats')
                index = prob_df.loc[:,'number'].isin(small_contigs.loc[:,'number'])
                prob_df.loc[index, 'Prob_Plasmid'] = float(prob_small_repeats) # OVERLAP!  

                # Transposases are also corner-cases for the machine-learning algorithm, we follow the same principle as done with short contigs 
                index = prob_df.loc[:,'number'].isin(repeats.loc[:,'number'])
                prob_df.loc[index, 'Prob_Plasmid'] = float(prob_small_repeats)

                final_probs = list(prob_df.loc[:,'Prob_Plasmid'])

                record_connections = []
                for i in range(len(final_probs)):
                    record_connections.append([1.0, number_iterations, iteration, elongation, initial_seed, seed, list_connections[i], final_probs[i]])

                record_connections = pd.DataFrame(data=record_connections, columns=['factor','number_iterations','iteration','elongation','first_node','ingoing_node','outgoing_node','Probability_pl_chr'])

                index = graph_contigs.loc[:,'number'].isin(record_connections.loc[:,'outgoing_node'])
                cov_connections_info = graph_contigs.loc[index,:].copy()

                up_cutoff = cov_connections_info.loc[:,'coverage'].copy() + max_variation            
                down_cutoff = cov_connections_info.loc[:,'coverage'].copy() - max_variation

                up_threshold = scipy.stats.norm.cdf(up_cutoff, loc=path_mean, scale=max_variation)
                down_threshold = scipy.stats.norm.cdf(down_cutoff, loc=path_mean, scale=max_variation)

                # Simple test
                window = up_threshold - down_threshold
                cov_connections_info.loc[:,'Probability_cov'] = abs(window)
                #TODO find a better way to do this merge
                record_connections = record_connections.merge(cov_connections_info[['number','Probability_cov']], left_on='outgoing_node', right_on='number')
                record_connections = record_connections.drop(columns='number')

                index = record_connections.loc[:,'outgoing_node'].isin(repeats_signed_nodes)
                record_connections.loc[index,'Probability_cov'] = float(prob_small_repeats)

                index = record_connections.loc[:,'outgoing_node'].isin(small_contigs_signed_nodes)
                record_connections.loc[index,'Probability_cov'] = float(prob_small_repeats)                        

                record_connections.loc[:,'Probability'] = record_connections.loc[:,'Probability_pl_chr'] * record_connections.loc[:,'Probability_cov']
                record_connections.loc[:,'Probability_freq'] = record_connections.loc[:,'Probability']/sum(record_connections.loc[:,'Probability'])
                record_connections.loc[:,'Verdict'] = 'non-selected'

                if(sum(record_connections.loc[:,'Probability'] >= filtering_threshold) == 0):
                    output = '\t'.join(path)
                    paths_list.append(output)
                    path = [initial_seed] # There are no connections possible from this contig
                    break # Exiting elongation loop

                # Filter step to avoid going into really bad connections
                index = list(record_connections.loc[:,'Probability'] >= filtering_threshold)
                filter_connections = record_connections.loc[index,:].copy()
                if(sum(filter_connections['Probability_freq']) != 1): # recalculate probability frequencies after filter
                    filter_connections['Probability_freq'] = filter_connections.loc[:,'Probability']/sum(filter_connections.loc[:,'Probability'])

                random_connection = np.random.choice(filter_connections.loc[:,'outgoing_node'], size=1, p=filter_connections.loc[:,'Probability_freq'])[0] # Choose one connection 
                path.append(str(random_connection))

                if(random_connection == path[0]):
                    output = '\t'.join(path)
                    paths_list.append(output)
                    path = [initial_seed] 
                    break # Exiting elongation loop

                seed = str(random_connection)

                if(len(path) >= number_nodes): # We only exit the function if we have reached the maximum number of nodes allowed per path
                    output = '\t'.join(path)
                    paths_list.append(output)
                    path = [initial_seed] 
                    break # Exiting elongation loop

                index = graph_contigs.loc[:,'number'].isin(path)
                info_path = graph_contigs.loc[index,:].copy()

                index = ~info_path.loc[:,'number'].isin(repeats_signed_nodes)
                info_path = info_path.loc[index,:]

                length_path = sum(info_path.loc[:,'length'])
                info_path.loc[:,'contribution'] = info_path.loc[:,'length'] / length_path
                path_mean = np.average(info_path.loc[:,'coverage'], weights=info_path.loc[:,'contribution']) # Coverage of the current path

        return paths_list


    final_paths = []
    for seed in initialize_nodes:
        np.random.seed(123)
        positive_seed = seed + '+'
        negative_seed = seed + '-'
        final_paths.extend(plasmid_graph(positive_seed, output_path, links, number_iterations, number_nodes, max_variation, filtering_threshold, prob_small_repeats, direction='forward'))
        final_paths.extend(plasmid_graph(negative_seed, output_path, links, number_iterations, number_nodes, max_variation, filtering_threshold, prob_small_repeats, direction='reverse'))

    with open(output_path, mode='w') as outfile:
        for path in final_paths:
            outfile.write(path + '\n')

    return
