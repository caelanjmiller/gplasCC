import shutil


def extract_nodes(sample, infile, minlen):
    output_links = f"gplas_input/{sample}_raw_links.txt"
    output_nodes = f"gplas_input/{sample}_raw_nodes.fasta"
    output_contigs = f"gplas_input/{sample}_contigs.fasta"

    with open(infile,'r') as graph, open(output_links,'w') as links, open(output_nodes,'w') as nodes, open(output_contigs,'w') as contigs:
        for line in graph:
            line = line.rstrip()
            if line[0] == 'S':
                cols = line.split('\t')
                number = cols[0] + str(cols[1])
                sequence = cols[2]
                if 'LN' in str(cols[3]): # Unicycler assembly
                    information = '_'.join(cols[3:])
                elif 'KC' in str(cols[3]): # Spades assembly
                    information = cols[3]
                else: # Empty sequence field (0 length), or other error in gfa format
                      # TODO exit the tool if the assembly is not from Unicycler or Spades? right now it would just skip all nodes and probably run into an error later in the workflow
                    continue # Skip node and continue to the next
                nodes.write(f">{number}_{information}\n{sequence}\n")
                if len(sequence) >= minlen:
                    contigs.write(f">{number}_{information}\n{sequence}\n")
            elif line[0] == 'L':
                links.write(f"{line}\n")


def extract_unbinned_solutions(sample):
    normal_results = f"results/normal_mode/{sample}_results_no_repeats.tab"
    bold_walks = f"walks/bold_mode/{sample}_solutions_bold.tab"
    unbinned_walks = f"walks/unbinned_nodes/{sample}_solutions_unbinned.tab"   
    #Get unbinned nodes
    unbinned_nodes = []
    with open(normal_results,'r') as file:
        for line in file:
            line = line.rstrip()
            cols = line.split('\t')
            component = cols[7]
            if component == 'Unbinned':
                number = str(cols[4])
                unbinned_nodes.append(number)
    #Select bold walks that start with unbinned nodes
    with open(bold_walks,'r') as infile, open(unbinned_walks,'w') as outfile:
        for line in infile:
            first_node = str(line.split('\t', 1)[0])
            first_node_unsigned = first_node.replace('+','')
            first_node_unsigned = first_node_unsigned.replace('-','')
            if first_node_unsigned in unbinned_nodes:
                outfile.write(line)
    #Combine solutions from bold and normal mode
    normal_walks = f"walks/normal_mode/{sample}_solutions.tab"
    combined_walks = f"walks/{sample}_solutions.tab"
    shutil.copyfile(normal_walks, combined_walks)
    with open(unbinned_walks,'r') as infile, open(combined_walks,'a') as outfile:
        for line in infile:
            outfile.write(line)
