#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#from cProfile import run
import shutil
import linecache
#import glob
import os
import sys
import argparse
#import subprocess
#from pathlib import Path
import glob
from .version import version as VERSION
#VERSION="1.0.0"
import time
##os.chdir("C:/Users/oscar/Documenten/UU/04_BiBc/6.2_Research_Profile/gplas/gplas-2-python/gplas/scripts")
##os.getcwd()

#TODO change script/function/import names
from gplas.scripts.m_node_extraction import extract_nodes
from gplas.scripts.m_node_extraction import extract_unbinned_solutions
from gplas.scripts.m_check_independent_prediction_format import check_prediction
from gplas.scripts.m_coverage import coverage
from gplas.scripts.m_paths import generate_paths
from gplas.scripts.m_coocurrence import calculate_coocurrence
from gplas.scripts.m_paths_repeats import generate_repeat_paths
from gplas.scripts.m_coocurrence_repeats import calculate_coocurrence_repeats
from gplas.scripts.m_run_plasmidcc import run_plasmidCC
from gplas.scripts import m_utils as utils

start_time = time.time()

# Directories
pkgdir = os.path.dirname(__file__)

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
        elif option_string == "--speciesopts":
            for species in utils.speciesopts:
                print(f"'{species}'")
        parser.exit()

#create a function to pass float ranges
parser = argparse.ArgumentParser(description='gplas: A tool for binning plasmid-predicted contigs into individual predictions.',formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=False)
parser.register('action','printing',PriorityPrinting)

inputgroup = parser.add_argument_group("General")
inputgroup.add_argument('-i', '--input', type=utils.is_valid_file, required=True, help="Path to the graph file in GFA (.gfa) format, used to extract nodes and links")
inputgroup.add_argument('-n', '--name', type=str, default='unnamed', help="Output name used in the gplas files")

classifiergroup = inputgroup.add_mutually_exclusive_group(required=True)
classifiergroup.add_argument('-s', '--species', type=utils.check_species, help="Choose a species database for plasmidCC classification. Use --speciesopts for a list of all supported species")
classifiergroup.add_argument('-P', '--prediction', type=utils.file_exists, help="If not using plasmidCC. Provide a path to an independent binary classification file")
classifiergroup.add_argument('--extract', action='store_true', help="extract FASTA sequences from the assembly graph to use with an external classifier")
#TODO for now we are missing the option to supply a custom centrifuge database directly via gplas

paramgroup = parser.add_argument_group("Parameters")
paramgroup.add_argument('-t', '--threshold_prediction', type=float, default=0.5, help="Prediction threshold for plasmid-derived sequences")
paramgroup.add_argument('-b', '--bold_walks', type=int, default=5, help="Coverage variance allowed for bold walks to recover unbinned plasmid-predicted nodes")
paramgroup.add_argument('-r', '--repeats_coverage_sd', type=int, default=2, help="Coverage variance allowed for assigning repeats to bins")
paramgroup.add_argument('-x', '--number_iterations', type=int, default=20,help="Number of walk iterations per starting node")
paramgroup.add_argument('-f', '--filt_gplas', type=float, default=0.1, help="filtering threshold to reject outgoing edges")
paramgroup.add_argument('-e', '--edge_threshold', type=float, default=0.1, help="Edge threshold")
paramgroup.add_argument('-q', '--modularity_threshold', type=float, default=0.2, help="Modularity threshold to split components in the plasmidome network")
paramgroup.add_argument('-l', '--length_filter', type=int, default=1000, help="Filtering threshold for sequence length")

utilgroup = parser.add_argument_group("Utility")
utilgroup.add_argument('-k', '--keep', action='store_true', help="Keep intermediary files")
#utilgroup.add_argument('--threads', type=int, default=1, help="Max number of threads to ")#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
utilgroup.add_argument('--silent', action='store_true', help="Suppress verbose print statements")
utilgroup.add_argument('--speciesopts', action='printing', nargs=0, help="Prints a list of all supported species for the -s flag")
utilgroup.add_argument('-h', '--help', action='printing', nargs=0, help="Prints this message")
utilgroup.add_argument('-v', '--version', action='printing', nargs=0, help="Prints gplas version")
args = parser.parse_args()

#Success Messages
def success_message():
    print('\n')
    print(read_logo)
    print(f"""
Congratulations! Prediction succesfully done.
Your results are in results/
We hope it helps your research, thank you for using gplas version {VERSION}

Please cite: https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
""")
    end_time = time.time()
    duration = end_time - start_time
    verbose_print(f"gplas took {round(duration,1)} seconds to run")
    sys.exit(0)

def success_message_extract():
    print('\n')
    print(read_logo)
    print(f"""
Congratulations! Your nodes have been succesfully extracted.
Your results are in gplas_input/{args.name}_contigs.fasta. Please, use an external tool to classify the nodes in this file, and then bin them into individual plasmids using gplas.
We hope it helps your research, thank you for using gplas version {VERSION}

Please cite: https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
""")
    sys.exit(0)

def verbose_print(message, end='\n'):
    if args.silent:
        return
    else:
        return print(message, end=end)

#******************************#
#*                            *#
#*        Start gplas         *#
#*                            *#
#******************************#

#Print messages
with open(f'{pkgdir}/../figures/logo.txt', 'r') as file: #TODO fix this path for final version
    read_logo = file.read()
print('\n')
print(read_logo)
print('\n')
print("##################################################################")

#Print chosen parameters
print("Your results will be named:...........................", args.name)
print("Input graph:..........................................", args.input)
print("Threshold for predicting plasmid-derived contigs:.....", args.threshold_prediction)
print("Number of plasmid walks created per node:.............", args.number_iterations)
print("Threshold of gplas scores:............................", args.filt_gplas)
print("Minimum frequency to consider an edge:................", args.edge_threshold)
print("Modularity threshold used to partition the network:...", args.modularity_threshold)
print("Coverage SD for bold mode:............................", args.bold_walks)
print("Coverage SD for repeats:..............................", args.repeats_coverage_sd)
print("Minimum sequence length:..............................", args.length_filter)
print("##################################################################\n")

##_1.0 Run analysis
#_1.1 Extract nodes and links from the assembly graph
verbose_print("Extracting contigs from the assembly graph...", end='\r')
os.makedirs("gplas_input", exist_ok=True)

extract_nodes(sample = args.name,
              infile = args.input,
              maxlen = args.length_filter)

utils.check_output(f"gplas_input/{args.name}_raw_nodes.fasta")
verbose_print("Extracting contigs from the assembly graph completed!")

#_1.2 If in extract mode, exit workflow after succesful extraction. Else continue workflow
if args.extract:
    success_message_extract() #Exits the workflow

##_2.0 Obtain correct prediction file
#_2.1 Run plasmidCC if no independent prediction file was given
##TODO fix plasmidCC FASTA input and only give it the sample_contigs.fasta as input to prevent double extracting of nodes
if args.species:
    os.makedirs("plasmidCC", exist_ok=True)
    
    run_plasmidCC(infile = args.input,
                  sample = args.name,
                  species = args.species,
                  maxlen = args.length_filter)
    
    #cleanup_centrifuge(sample = args.name)
    
    path_prediction = f"plasmidCC/{args.name}/{args.name}_gplas.tsv"
else:
    path_prediction = args.prediction

#_2.2 Check if the prediction file is correctly formatted.
verbose_print("Checking prediction file format...", end='\r')

check_prediction(sample = args.name,
                 path_prediction = path_prediction)

verbose_print("Checking prediction file format completed!")

##_3.0 Run gplas in normal mode
#_3.1 Extract nodes/links from the assembly graph
verbose_print("Calculating base coverage...", end='\r')
os.makedirs("coverage", exist_ok=True)

coverage(sample = args.name,
         path_prediction = path_prediction,
         pred_threshold = args.threshold_prediction)

verbose_print("Calculating base coverage completed!")

#_3.2 Generate random walks
verbose_print("Generating random walks in normal mode...", end='\r')
os.makedirs("walks/normal_mode", exist_ok=True)

#TODO this will append paths/connections to previous file if using the same sample name
##instead of appending to file each loop, append to list and all the way at the end write list of lists to file?
##also needs to use multiprocessing which complicates things
generate_paths(sample = args.name,
               number_iterations = args.number_iterations,
               filt_threshold = args.filt_gplas,
               mode = "normal",
               sd_coverage = 1)

verbose_print("Generating random walks in normal mode completed!")

#_3.3 Calculate coocurrence between walks
verbose_print("Calculating coocurrence of random walks...", end='\r')
os.makedirs("results/normal_mode", exist_ok=True)

#TODO coocurrence script breaks if reruning gplas with the same sample name
##see generate_paths() appending to old file; coocurrence doesnt break if you remove the previous 'solutions' file
calculate_coocurrence(sample = args.name,
                      number_iterations = args.number_iterations,
                      pred_threshold = args.threshold_prediction,
                      mod_threshold = args.modularity_threshold,
                      mode = "normal")

utils.check_output(f"results/normal_mode/{args.name}_results_no_repeats.tab")
verbose_print("Calculating coocurrence of random walks completed!")

##_4.0 Resolve unbinned contigs
#_4.1 Check for unbinned contigs
unbinned_path=f'results/normal_mode/{args.name}_bin_Unbinned.fasta'
if os.path.exists(unbinned_path):
    #_4.1.1 Run gplas in bold mode if contigs were left unbinned
    verbose_print("Some contigs were left unbinned")#improve tell user how many contigs are unbinned?
    #_4.1.1.1 Generate random walks
    verbose_print("Generating random walks in bold mode...", end='\r')
    os.makedirs("walks/bold_mode", exist_ok=True)

    generate_paths(sample = args.name,
                   number_iterations = args.number_iterations,
                   filt_threshold = args.filt_gplas,
                   mode = "bold",
                   sd_coverage = args.bold_walks)
    
    verbose_print("Generating random walks in bold mode completed!")
    
    #_4.1.1.2 Extract unbinned solutions
    verbose_print("Extracting unbinned contigs from bold walks...", end='\r')
    os.makedirs("walks/unbinned_nodes", exist_ok=True)

    extract_unbinned_solutions(sample = args.name)
    
    verbose_print("Extracting unbinned contigs from bold walks completed!")

    #_4.1.1.3 Recalculate coocurrence of walks using the combined solutions
    verbose_print("Recalculating coocurrence of random walks...", end='\r')
    
    calculate_coocurrence(sample = args.name,
                          number_iterations = args.number_iterations,
                          pred_threshold = args.threshold_prediction,
                          mod_threshold = args.modularity_threshold,
                          mode = "unbinned")

    verbose_print("Recalculating coocurrence of random walks completed!")
    
#_4.1.2 Copy files from normal mode if there were no unbinned contigs
else:
    for file in glob.glob(f"results/normal_mode/{args.name}*"):
        shutil.copy(file, "results/")

utils.check_output(f"results/{args.name}_results_no_repeats.tab")

##_5.0 Add repeated elements
#_5.1 Check for repeats
repeated_elements_path=f"coverage/{args.name}_repeat_nodes.tab"
line_content=linecache.getline(repeated_elements_path,1)
if line_content:
    #_5.1.1 Run gplas on repeated elements
    verbose_print("Adding repeated elements to the predictions...", end='\r')
    os.makedirs("walks/repeats", exist_ok=True)
    
    #_5.1.1.1 Generate random walks
    generate_repeat_paths(sample = args.name,
                          number_iterations = args.number_iterations,
                          filt_threshold = args.filt_gplas,
                          sd_coverage = args.repeats_coverage_sd)

    #_5.1.1.2 Calculate coocurrence between walks
    calculate_coocurrence_repeats(sample = args.name,
                                  number_iterations = args.number_iterations,
                                  pred_threshold = args.threshold_prediction,
                                  mod_threshold = args.modularity_threshold,
                                  sd_coverage = args.repeats_coverage_sd)
    
    verbose_print("Adding repeated elements to the predictions completed!")


#_5.1.2 If there are no repeated elements, just rename the results files.
else:
    shutil.move("results/{args.name}_results_no_repeats.tab", "results/{args.name}_results.tab")
    shutil.move("results/{args.name}_bins_no_repeats.tab", "results/{args.name}_bins.tab")

utils.check_output(f"results/{args.name}_results.tab")

##_6.0 If the -k flag was not selected, delete intermediary files
if args.keep==False and args.extract==False:
    verbose_print("Intermediate files will be deleted. If you want to keep these files, use the -k flag")
    #Coverage files
    utils.delete_file(f"coverage/{args.name}_clean_links.tab")
    utils.delete_file(f"coverage/{args.name}_clean_prediction.tab")
    utils.delete_file(f"coverage/{args.name}_clean_repeats.tab")
    utils.delete_file(f"coverage/{args.name}_estimation.txt")
    utils.delete_file(f"coverage/{args.name}_graph_contigs.tab")
    utils.delete_file(f"coverage/{args.name}_initialize_nodes.tab")
    utils.delete_file(f"coverage/{args.name}_isolated_nodes.tab")
    utils.delete_file(f"coverage/{args.name}_repeat_nodes.tab")
    utils.delete_file(f"coverage/{args.name}_repeats_graph.tab")
    #Walks normal mode
    utils.delete_file(f"walks/normal_mode/{args.name}_solutions.tab")
    #Walks bold mode + unbinned solutions
    utils.delete_file(f"walks/bold_mode/{args.name}_solutions_bold.tab")
    utils.delete_file(f"walks/unbinned_nodes/{args.name}_solutions_unbinned.tab")
    utils.delete_file(f"walks/{args.name}_solutions.tab")
    #Walks repeats
    utils.delete_file(f"walks/repeats/{args.name}_solutions.tab")
    #Results no_repeats
    utils.delete_file(f"results/{args.name}_results_no_repeats.tab")
    utils.delete_file(f"results/{args.name}_bins_no_repeats.tab")
    
    ##Delete directories if they exist and are empty
    utils.delete_empty_dir("coverage/")
    utils.delete_empty_dir("walks/normal_mode/")
    utils.delete_empty_dir("walks/bold_mode/")
    utils.delete_empty_dir("walks/unbinned_nodes/")
    utils.delete_empty_dir("walks/repeats/")
    utils.delete_empty_dir("walks/")
    
##_7.0 Show success message and exit workflow
success_message()
