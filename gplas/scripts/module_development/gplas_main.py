#!/usr/bin/env python3
# -*- coding: utf-8 -*-

if __name__ == "__main__":

    #from cProfile import run
    import shutil
    import linecache
    import glob
    import os
    import sys
    import argparse
    #import time
    import subprocess
    #from pathlib import Path
    #from .version import version as VERSION
    VERSION="1.0.0"
    ##os.chdir("C:/Users/oscar/Documenten/UU/04_BiBc/6.2_Research_Profile/gplas/gplas-2-python/gplas/scripts/module_development")
    ##os.getcwd()
    
    from m_check_independent_prediction_format import check_prediction
    from m_coverage import coverage

    ##workflow layout
    #0_mkdirs
    #1_awk_links
    #2_awk_nodes_extract_nodes
    #3_awk_nodes_filter
    #4_awk_nodes_rename

    #coverage(sample=sample, path_prediction=path_prediction, classifier=classifier, threshold=threshold)

    
    #m_paths
    #m_paths_bold
    #m_coocurrence
    
    #5_extract_unbinned_nodes
    #6_combine_solutions
    
    #m_coocurrence_final
    #m_paths_repeats
    #m_coocurrence_repeats
    
    #remove_intermediate_files
    
    
    
    # Directories
    pkgdir = os.path.dirname(__file__)
    #scriptdir = f"{pkgdir}/scripts/module_development"
    scriptdir = pkgdir
    
    #******************************#
    #*                            *#
    #* Command line parsing       *#
    #*                            *#
    #******************************#

    class PriorityPrinting(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if option_string == "-h" or option_string == "--help":
                parser.print_help()
            elif option_string == "-v" or option_string == "--version":
                print(f"gplas version {VERSION}")
            parser.exit()

    #create a function to pass float ranges
    parser = argparse.ArgumentParser(description='gplas: A tool for binning plasmid-predicted contigs into individual predictions.',formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=False)
    parser.register('action','printing',PriorityPrinting)
    parser.add_argument('-i', '--input', type=str, required=True, help="Path to the graph file in GFA (.gfa) format, used to extract nodes and links")
    parser.add_argument('-c', '--classifier', type=str, required=True, choices=['extract','predict'], help="Select to extract nodes from the assembly graph or to predict individual plasmids.")
    parser.add_argument('-n', '--name', type=str, default='unnamed', help="Output name used in the gplas files")
    parser.add_argument('-P', '--prediction', help="Path to the binary classification file")
    parser.add_argument('-k', '--keep', action='store_true', help="Keep intermediary files")
    parser.add_argument('-t', '--threshold_prediction', type=float, default=0.5, help="Prediction threshold for plasmid-derived sequences")
    parser.add_argument('-b', '--bold_walks',type=int, default=5, help="Coverage variance allowed for bold walks to recover unbinned plasmid-predicted nodes")
    parser.add_argument('-r', '--repeats_coverage_sd',type=int, default=2, help="Coverage variance allowed for assigning repeats to bins")
    parser.add_argument('-x', '--number_iterations',type=int, default=20,help="Number of walk iterations per starting node")
    parser.add_argument('-f', '--filt_gplas', type=float, default=0.1, help="filtering threshold to reject outgoing edges")
    parser.add_argument('-e', '--edge_threshold', type=float, default=0.1, help="Edge threshold")
    parser.add_argument('-q', '--modularity_threshold', type=float, default=0.2, help="Modularity threshold to split components in the plasmidome network")
    parser.add_argument('-l', '--length_filter', type=int, default=1000, help="Filtering threshold for sequence length")
    parser.add_argument('-h', '--help', action='printing', nargs=0, help="Prints this message")
    parser.add_argument('-v', '--version', action='printing', nargs=0, help="Prints gplas version")
    args = parser.parse_args()

    #Success Message
    def success_message():
        print ('\n')
        print(read_logo)
        print ('\n')
        print(f"""

    Congratulations! Prediction succesfully done.

    Your results are in results/

    We hope it helps your research, thanks for using gplas version {VERSION}!

    Please cite: https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
    """)

        sys.exit(0)

    def success_message_extract():
        print ('\n')
        print(read_logo)
        print ('\n')
        print(f"""

    Congratulations! Your nodes have been succesfully extracted.

    Your results are in gplas_input/{args.name}_contigs.fasta. Please, use an external tool to classify the nodes in this file, and then bin them into individual plasmids using gplas.

    We hope it helps your research, thanks for using gplas version {VERSION}!.

    Please cite: https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
    """)

        sys.exit(0)

    #Prediction not successful
    def error_message():
        print('\n')
        print(f"""
    Looks like no plasmids could be detected in your assembly graph

    Please check the file:   coverage/{args.name}_clean_prediction.tab.
    If all contigs were predicted as chromosome, gplas probably failed at the step to create random walks starting from plasmid seeds. If that's the case, probably your isolate does not carry any plasmid(s)
    If you don't see any files present at:   gplas_input/  or  coverage/  most likely the installation of gplas failed at some point


            """)
        sys.exit(0)

    def error_message_extract():
        print('\n')
        print("""
    Looks like the nodes were not extracted. Please, check above for error messages.

            """)
        sys.exit(0)

    #******************************#
    #*                            *#
    #*        Start gplas         *#
    #*                            *#
    #******************************#

    #Print messages
    with open(f'{pkgdir}/figures/logo.txt', 'r') as file:
        read_logo = file.read()
    print ('\n')
    print(read_logo)
    print("##################################################################")

    #Print chosen parameters
    print("Your results will be named:...........................", args.name)
    print("Input graph:..........................................", args.input)
    print("Classifier:...........................................", args.classifier)
    print("Threshold for predicting plasmid-derived contigs:.....", args.threshold_prediction)
    print("Number of plasmid walks created per node:.............", args.number_iterations)
    print("Threshold of gplas scores:............................", args.filt_gplas)
    print("Minimum frequency to consider an edge:................", args.edge_threshold)
    print("Modularity threshold used to partition the network:...", args.modularity_threshold)
    print("Coverage SD for bold mode:............................", args.bold_walks)
    print("Coverage SD for repeats:..............................", args.repeats_coverage_sd)
    print("Minimum sequence length:..............................", args.length_filter)
    print("##################################################################")

    """
    #3. Run analysis
    ##3.1  Extract nodes and links from the assembly graph
    print("We first need to extract the contigs from the assembly graph, these contigs are later used for your binary prediction.\n")
    extract_links_command=f'bash {scriptdir}/1_awk_links.sh {args.name} {args.input}'
    subprocess.run(extract_links_command, shell=True, text=True, executable='/bin/bash')
    #extract_nodes_command=f'bash {scriptdir}/2_awk_nodes_extract_nodes.sh'
    #subprocess.run(extract_nodes_command, shell=True, text=True, executable='/bin/bash')
    #filter_nodes_command=f'bash {scriptdir}/3_awk_nodes_filter.sh'
    #subprocess.run(filter_nodes_command, shell=True, text=True, executable='/bin/bash')
    #rename_nodes_command=f'bash {scriptdir}/4_awk_nodes_rename.sh'
    #subprocess.run(rename_nodes_command, shell=True, text=True, executable='/bin/bash')
    
    ##3.2.1 If in extract mode, exit workflow after extraction is complete
    if args.classifier == "extract":
        pass  
    
    ##3.2.2 If not in extract mode, run workflow until the _raw_nodes.fasta file is obtainted.
    else:
        
        #ddd

        #3.4 Check if the independent prediction file is correctly formatted.
        print("Checking if prediction file is correctly formatted.\n")
        #check_prediction(sample=sample, path_prediction=path_prediction)
    
        #if file is correctly formated, continue the workflow
        #ddd
        
        #3.5 Check for Unbinned contigs
        unbinned_path=f'results/normal_mode/{args.name}_bin_Unbinned.fasta'
        if os.path.exists(unbinned_path):
            ##3.5.1 run bold mode if contigs were left unbinned.
            print('\n')
            print('Some contigs were left Unbinned, running gplas in bold mode')
            print('\n')
            #ddd
            
        else:
            for file in glob.glob(f"results/normal_mode/{args.name}*"):
                shutil.copy(file, "results/")
        
    ## 3.6 ADD REPEATED ELEMENTS.
    ##3.6.1 Now add the repeats to the final bins
    repeated_elements_path=f'coverage/{args.name}_repeat_nodes.tab'
    line_content=linecache.getline(repeated_elements_path,2)
    if line_content:
        print('\n')
        print('Adding repeated elements to the predictions')
        print('\n')
        #ddd
    
    ##3.6.2 If there are not repeated elementss, just rename the results files.
    else:
        shutil.move("results/{args.name}_results_no_repeats.tab", "results/{args.name}_results.tab")
        shutil.move("results/{args.name}_bins_no_repeats.tab", "results/{args.name}_bins.tab")

    ##3.7 If the -k flag was not selected, delete intermediary files
    if args.keep==False and args.classifier!='extract':
      print("Intermediate files will be deleted. If you want to keep these files, use the -k flag")
      remove_command=f'bash {scriptdir}/remove_intermediate_files.sh -n {args.name}'
      subprocess.run(remove_command, shell=True, text=True, executable='/bin/bash')

    ##3.8 Check that output has been correctly created
    final_results_path=f'results/{args.name}_results.tab'
    raw_nodes_path=f'gplas_input/{args.name}_contigs.fasta'

    if args.classifier!='extract':
        if os.path.exists(final_results_path):
            success_message()
            sys.exit(0)
        else:
            error_message()
            sys.exit(1)
    else:
        if os.path.exists(raw_nodes_path):
            success_message_extract()
            sys.exit(0)
        else:
            error_message_extract()
            sys.exit(1)
    """

    #TODO is this ever used?
    def start():
        print("Starting gplas")
