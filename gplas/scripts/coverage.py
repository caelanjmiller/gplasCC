import pandas as pd
import statistics
from Bio.SeqIO.FastaIO import SimpleFastaParser

def coverage(sample, path_prediction, pred_threshold):
    #Inputs    
    path_nodes = f"gplas_input/{sample}_raw_nodes.fasta"
    path_links = f"gplas_input/{sample}_raw_links.txt"
    
    #Outputs
    output_graph_contigs = f"coverage/{sample}_graph_contigs.tab"
    output_clean_links = f"coverage/{sample}_clean_links.tab"
    output_graph_repeats = f"coverage/{sample}_repeats_graph.tab"
    output_clean_prediction = f"coverage/{sample}_clean_prediction.tab"
    output_isolated_nodes = f"coverage/{sample}_isolated_nodes.tab"
    output_clean_repeats = f"coverage/{sample}_clean_repeats.tab"
    output_initialize_nodes = f"coverage/{sample}_initialize_nodes.tab"
    output_repeat_nodes = f"coverage/{sample}_repeat_nodes.tab"
    output_cov_estimate = f"coverage/{sample}_estimation.txt"
        
    with open(path_nodes) as file:
        raw_nodes = [[str(values[0]), str(values[1])] for values in SimpleFastaParser(file)]
    
    raw_contig_names = [str(entry[0]) for entry in raw_nodes]
    
    kc_check = sum([name.count('KC') for name in raw_contig_names])
    
    if kc_check == len(raw_contig_names):
        lengths = [len(entry[1]) for entry in raw_nodes]
        kc_counts = [name.split(':', maxsplit=4)[2] for name in raw_contig_names]
        kc_counts = [int(name.replace('_','')) for name in kc_counts]
        kc_coverage = [kc/length for kc, length in zip(kc_counts, lengths)]
        coverage = [coverage/statistics.median(kc_coverage) for coverage in kc_coverage]
    else:
        raw_lengths = [name.split(':')[2] for name in raw_contig_names]
        lengths = [int(name.replace('_dp','')) for name in raw_lengths]
        coverage = [float(name.split(':')[4]) for name in raw_contig_names]

    raw_number = [name.split('_')[0] for name in raw_contig_names]
    number = [name.replace('S','') for name in raw_number]

    contig_info = pd.DataFrame(data={'number':number,
                                     'length':lengths,
                                     'coverage':coverage,
                                     'Contig_name':raw_contig_names})

    graph_pos_contigs = pd.DataFrame(data={'number':[digit+'+' for digit in number],
                                           'length':lengths,
                                           'coverage':coverage,
                                           'Contig_name':raw_contig_names})

    graph_neg_contigs = pd.DataFrame(data={'number':[digit+'-' for digit in number],
                                           'length':lengths,
                                           'coverage':coverage,
                                           'Contig_name':raw_contig_names})

    graph_contigs = pd.concat([graph_pos_contigs, graph_neg_contigs], ignore_index=True)
    graph_contigs.to_csv(output_graph_contigs, sep='\t', index=False, mode='w')

    raw_links = pd.read_csv(path_links, sep='\t', header=None)
    raw_links.rename({0:'L',1:'first_node',2:'first_sign',3:'second_node',4:'second_sign',5:'OM'}, axis='columns', inplace=True)

    links = []
    for index, row in raw_links.iterrows():
        reverse_first_sign = '+' if row['first_sign'] == '-' else '-'
        reverse_second_sign = '+' if row['second_sign'] == '-' else '-'

        clean_first = str(row['first_node'])+str(row['first_sign'])
        clean_second = str(row['second_node'])+str(row['second_sign'])
        info_forward = [clean_first, 'to', clean_second]
        links.append(info_forward)

        clean_rev_first = str(row['second_node'])+str(reverse_second_sign)
        clean_rev_second = str(row['first_node'])+str(reverse_first_sign)
        info_reverse = [clean_rev_first, 'to', clean_rev_second]
        links.append(info_reverse)

    links = pd.DataFrame(links)
    links.to_csv(output_clean_links, sep='\t', index=False, header=False, mode='w')

    unique_nodes = sorted(list(set(links[0])))

    outdegree_info = []
    for node in unique_nodes:
        repeat_links = links[links[0] == node]
        unique_links = repeat_links.drop_duplicates()
        outdegree_info.append([node, len(unique_links[2])])

    outdegree_info = pd.DataFrame(outdegree_info, columns=['number', 'outdegree'])

    indegree_info = []
    for node in unique_nodes:
        repeat_links = links[links[2] == node]
        unique_links = repeat_links.drop_duplicates()
        indegree_info.append([node, len(unique_links[0])])

    indegree_info = pd.DataFrame(indegree_info, columns=['number', 'indegree'])

    repeat_info = pd.merge(outdegree_info, indegree_info, on='number')
    repeats = repeat_info[(repeat_info['indegree'] > 1) | (repeat_info['outdegree'] > 1)]

    repeats.to_csv(output_graph_repeats, sep='\t', index=False, mode='w')

    repeats.loc[:,'number'] = [name.replace('+','') for name in repeats['number']]
    repeats.loc[:,'number'] = [name.replace('-','') for name in repeats['number']]

    pred = pd.read_table(path_prediction, sep='\t', header=0)

    raw_number = [name.split('_', maxsplit=1)[0] for name in pred['Contig_name']]
    pred.loc[:,'number'] = [number.replace('S','') for number in raw_number]
    pred = pd.merge(pred, contig_info, on=['Contig_name','number'])

    final_prediction = pred[[number not in list(repeats['number']) for number in pred['number']]]

    final_prediction.to_csv(output_clean_prediction, sep='\t', index=False, mode='w')

    unique_nodes_signless = links[0]
    unique_nodes_signless = [node.replace('+','') for node in unique_nodes_signless]
    unique_nodes_signless = [node.replace('-','') for node in unique_nodes_signless]

    isolated_nodes = contig_info[[number not in unique_nodes_signless for number in contig_info['number']]]
    isolated_nodes = pd.merge(pred, isolated_nodes['Contig_name'], on='Contig_name')
    isolated_nodes = isolated_nodes[isolated_nodes['Prob_Plasmid'] >= pred_threshold]

    isolated_nodes.to_csv(output_isolated_nodes, sep='\t', index=False, mode='w')

    repeats_final = pred[[number in list(repeats['number']) for number in pred['number']]]

    repeats_final.to_csv(output_clean_repeats, sep='\t', index=False, mode='w')

    pl_nodes = final_prediction[final_prediction['Prob_Plasmid'] >= pred_threshold]
    pl_nodes = pl_nodes[pl_nodes['Contig_length'] > 500]  # TODO ASK Julian: is this threshold needed?
    index = [number not in list(repeats['number']) for number in pl_nodes['number']]
    pl_nodes = pl_nodes.loc[index,:]

    initialize_nodes = sorted(list(set(pl_nodes['number'])))
    with open(output_initialize_nodes, 'w') as file:
        for node in initialize_nodes:
            file.write(node + '\n')

    repeats_nodes = sorted(list(set(repeats_final['number'])))
    with open(output_repeat_nodes, 'w') as file:
        for node in repeats_nodes:
            file.write(node + '\n')

    chr_contigs = pred[(pred['Prob_Chromosome'] > 0.7) & (pred['Contig_length'] > 1000)]  # TODO this is hardcoded, both values could/should become parameters?
    cov_estimation = chr_contigs[[number not in list(repeats['number']) for number in chr_contigs['number']]]
    sd_estimation = statistics.stdev(cov_estimation['coverage'])

    with open(output_cov_estimate, 'w') as file:
        file.write(str(sd_estimation))

    return
