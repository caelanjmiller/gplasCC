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
from pathlib import Path
from .version import version as VERSION
#VERSION="1.0.0"
import time
##os.chdir("C:/Users/oscar/Documenten/UU/04_BiBc/6.2_Research_Profile/gplas/gplas-2-python/gplas/scripts")
##os.getcwd()

#TODO change script/function/import names
from gplas.scripts.m_check_independent_prediction_format import check_prediction
from gplas.scripts.m_coverage import coverage
from gplas.scripts.m_paths import generate_paths
from gplas.scripts.m_coocurrence import calculate_coocurrence
from gplas.scripts.m_paths_repeats import generate_repeat_paths
from gplas.scripts.m_coocurrence_repeats import calculate_coocurrence_repeats
from gplas.scripts.m_run_plasmidcc import run_plasmidCC
#from plasmidCC.scripts import utils as utilsCC

start_time = time.time()

# Directories
pkgdir = os.path.dirname(__file__)
#pkgdir = Path(__file__).parent.resolve() #improve for if we want to repace os with pathlib?
#scriptdir = f"{pkgdir}/scripts/module_development"
#scriptdir = f"{pkgdir}/scripts"

#******************************#
#*                            *#
#* Command line parsing       *#
#*                            *#
#******************************#

#speciesopts needs to be manually updated if plasmidCC ever changes their species options
speciesopts = ['General','Escherichia coli','Enterococcus faecium','Enterococcus faecalis','Salmonella enterica','Staphylococcus aureus','Acinetobacter baumannii','Klebsiella pneumoniae']
def check_species(arg):
    if not arg in speciesopts:
        raise argparse.ArgumentTypeError(f"\'{arg}\' is not a recognised species" + "\nUse gplas with the --speciesopts flag for a list of all supported species")
    return arg

class PriorityPrinting(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if option_string == "-h" or option_string == "--help":
            parser.print_help()
        elif option_string == "-v" or option_string == "--version":
            print(f"gplas version {VERSION}")
        elif option_string == "--speciesopts":
            for species in speciesopts:
                print(f"\'{species}\'")
        parser.exit()

#create a function to pass float ranges
parser = argparse.ArgumentParser(description='gplas: A tool for binning plasmid-predicted contigs into individual predictions.',formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=False)
parser.register('action','printing',PriorityPrinting)

inputgroup = parser.add_argument_group("General")
inputgroup.add_argument('-i', '--input', type=str, required=True, help="Path to the graph file in GFA (.gfa) format, used to extract nodes and links")
inputgroup.add_argument('-m', '--mode', type=str, required=True, choices=['extract','predict'], help="extract: extract FASTA sequences from the assembly graph to use with an external classifier | predict: run the full gplas workflow to predict individual plasmids from the assembly graph")
inputgroup.add_argument('-n', '--name', type=str, default='unnamed', help="Output name used in the gplas files")

classifiergroup = inputgroup.add_mutually_exclusive_group(required=True)
classifiergroup.add_argument('-s', '--species', type=check_species, help="Choose a species database for plasmidCC classification. Use --speciesopts for a list of all supported species")
classifiergroup.add_argument('-P', '--prediction', type=str, help="If not using plasmidCC. Provide a path to an independent binary classification file")
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

#Success Message
def success_message():
    print('\n')
    print(read_logo)
    #print('\n')
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
    #print('\n')
    print(f"""
Congratulations! Your nodes have been succesfully extracted.
Your results are in gplas_input/{args.name}_contigs.fasta. Please, use an external tool to classify the nodes in this file, and then bin them into individual plasmids using gplas.
We hope it helps your research, thank you for using gplas version {VERSION}

Please cite: https://academic.oup.com/bioinformatics/article/36/12/3874/5818483
""")

    sys.exit(0)

#Prediction not successful
def error_message():
    print('\n')
    print(f"""
Looks like no plasmids could be detected in your assembly graph

Please check the file:   coverage/{args.name}_clean_prediction.tab
If all contigs were predicted as chromosome, gplas probably failed at the step to create random walks starting from plasmid seeds. If that's the case, probably your isolate does not carry any plasmid(s)
If you don't see any files present at:  gplas_input/  or  coverage/  most likely the installation of gplas failed at some point
""")
    sys.exit(0)

def error_message_extract():
    print('\n')
    print("""
Looks like the nodes were not extracted. Please, check above for error messages.
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

if args.prediction:
    classifier = 'external'
else:
    classifier = 'plasmidCC'

#3. Run analysis
verbose_print("Extracting contigs from the assembly graph...", end='\r')
Path("gplas_input").mkdir(parents=True, exist_ok=True)

##3.1  Extract nodes and links from the assembly graph
output_links = f"gplas_input/{args.name}_raw_links.txt"
output_nodes = f"gplas_input/{args.name}_raw_nodes.fasta"
output_contigs = f"gplas_input/{args.name}_contigs.fasta"

#TODO turn this into an extract() function
with open(args.input,"r") as graph, open(output_links,"w") as links, open(output_nodes,"w") as nodes, open(output_contigs,"w") as contigs:
    for line in graph:
        line = line.rstrip()
        if line[0] == "S":
            cols = line.split("\t")
            number = cols[0] + str(cols[1])
            sequence = cols[2]
            if len(cols) > 4: #unicycler
                information = "_".join(cols[3:])
            else: #spades
                information = cols[3]
            nodes.write(f">{number}_{information}\n{sequence}\n")
            if len(sequence) >= args.length_filter:
                contigs.write(f">{number}_{information}\n{sequence}\n")
        elif line[0] == "L":
            links.write(f"{line}\n")

#Check if output has been correctly created
#Path(f"gplas_input/{args.name}_raw_nodes.fasta").exists() #improve for if we want to replace os with pathlib?
#TODO turn this into a check_output() function
if os.path.exists(f"gplas_input/{args.name}_raw_nodes.fasta") == False:
    error_message_extract()
    sys.exit(1)    

verbose_print("Extracting contigs from the assembly graph completed!")

##3.2 If in extract mode, exit workflow after succesful extraction. Else continue workflow
if args.mode == "extract":
    success_message_extract()
    sys.exit(0)

##3.CC Run plasmidCC if no independent prediction file was given
##TODO fix plasmidCC FASTA input and only give it the sample_contigs.fasta as input to prevent double extracting of nodes
if classifier == 'plasmidCC':
    Path("plasmidCC").mkdir(parents=True, exist_ok=True)    
    run_plasmidCC(infile = args.input,
                  sample = args.name,
                  species = args.species,
                  maxlen = args.length_filter)
    #cleanup_centrifuge(sample = args.name)
    path_prediction = f"plasmidCC/{args.name}/{args.name}_gplas.tsv"
else:
    path_prediction = args.prediction

##3.3 Check if the prediction file is correctly formatted.
verbose_print("Checking prediction file format...", end='\r')
try:
    check_prediction(sample = args.name, path_prediction = path_prediction)
except Exception as e:
    print(e)
    sys.exit(1)
    
verbose_print("Checking prediction file format completed!")

##3.4 Run gplas in normal mode
##3.4.1 Extract information from the assembly graph
verbose_print("Calculating base coverage...", end='\r')
Path("coverage").mkdir(parents=True, exist_ok=True)

coverage(sample = args.name,
         path_prediction = path_prediction,
         pred_threshold = args.threshold_prediction)
verbose_print("Calculating base coverage completed!")

##3.4.2 Generate random walks
verbose_print("Generating random walks in normal mode...", end='\r')
Path("walks/normal_mode").mkdir(parents=True, exist_ok=True)

#TODO this will append paths/connections to previous file if using the same sample name
##instead of appending to file each loop, append to list and all the way at the end write list of lists to file?
##also needs to use multiprocessing which complicates things
generate_paths(sample = args.name,
               number_iterations = args.number_iterations,
               filt_threshold = args.filt_gplas,
               mode = "normal",
               sd_coverage = 1)
verbose_print("Generating random walks in normal mode completed!")

##3.4.3 Calculate coocurrence between walks
verbose_print("Calculating coocurrence of random walks...", end='\r')
Path("results/normal_mode").mkdir(parents=True, exist_ok=True)

#TODO coocurrence script breaks if reruning gplas with the same sample name
##see generate_paths() appending to old file; coocurrence doesnt break if you remove the previous 'solutions' file
calculate_coocurrence(sample = args.name,
                      number_iterations = args.number_iterations,
                      pred_threshold = args.threshold_prediction,
                      mod_threshold = args.modularity_threshold,
                      mode = "normal")

#Check if output has been correctly created
if os.path.exists(f"results/normal_mode/{args.name}_results_no_repeats.tab") == False:
    #make this also an error_message() function??
    sys.exit("ERROR: Something went wrong while running gplas in normal mode")

verbose_print("Calculating coocurrence of random walks completed!")

#3.5 Check for Unbinned contigs
unbinned_path=f'results/normal_mode/{args.name}_bin_Unbinned.fasta'
if os.path.exists(unbinned_path):
    ##3.5.1 Run bold mode if contigs were left unbinned.
    verbose_print("Some contigs were left unbinned")#improve tell user how many contigs are unbinned?
    verbose_print("Generating random walks in bold mode...", end='\r')
    Path("walks/bold_mode").mkdir(parents=True, exist_ok=True)

    generate_paths(sample = args.name,
                   number_iterations = args.number_iterations,
                   filt_threshold = args.filt_gplas,
                   mode = "bold",
                   sd_coverage = args.bold_walks)
    verbose_print("Generating random walks in bold mode completed!")
    
    ##3.5.2 Extract unbinned solutions
    verbose_print("Extracting unbinned contigs from bold walks...", end='\r')
    Path("walks/unbinned_nodes").mkdir(parents=True, exist_ok=True)
    
    normal_results = f"results/normal_mode/{args.name}_results_no_repeats.tab"
    bold_walks = f"walks/bold_mode/{args.name}_solutions_bold.tab"
    unbinned_walks = f"walks/unbinned_nodes/{args.name}_solutions_unbinned.tab"   
    #Get unbinned nodes
    #improve turn this into a extract_unbinned() / combine_solutions() function
    unbinned_nodes = []
    with open(normal_results,"r") as file:
        for line in file:
            line = line.rstrip()
            cols = line.split("\t")
            number = str(cols[4]) #improve move number= to after the if component== statement; should be (slightly) more efficient
            component = cols[7]
            if component == "Unbinned":
                unbinned_nodes.append(number)
    #Select bold walks that start with unbinned nodes
    with open(bold_walks,"r") as infile, open(unbinned_walks,"w") as outfile:
        for line in infile:
            first_node = str(line.split("\t", 1)[0])
            first_node_unsigned = first_node[:-1] #improve find a better way to remove +/- from end of number?
            if first_node_unsigned in unbinned_nodes:
                outfile.write(line)
    ##3.5.3 Combine solutions from bold and normal mode
    normal_walks = f"walks/normal_mode/{args.name}_solutions.tab"
    combined_walks = f"walks/{args.name}_solutions.tab"
    shutil.copyfile(normal_walks, combined_walks)
    with open(unbinned_walks,"r") as infile, open(combined_walks,"a") as outfile:
        for line in infile:
            outfile.write(line)
    
    verbose_print("Extracting unbinned contigs from bold walks completed!")

    ##3.5.4 Recalculate coocurrence of walks using the combined solutions
    verbose_print("Recalculating coocurrence of random walks...", end='\r')
    calculate_coocurrence(sample = args.name,
                          number_iterations = args.number_iterations,
                          pred_threshold = args.threshold_prediction,
                          mod_threshold = args.modularity_threshold,
                          mode = "unbinned")

else:
    for file in Path("results/normal_mode/").glob(f"{args.name}*"):
        shutil.copy(file, "results/")

#Check if output has been correctly created
if os.path.exists(f"results/{args.name}_results_no_repeats.tab") == False:
    #make this also an error_message() function??
    sys.exit("ERROR: Something went wrong while running gplas in bold mode")

verbose_print("Recalculating coocurrence of random walks completed!")

## 3.6 ADD REPEATED ELEMENTS.
##3.6.1 Now add the repeats to the final bins
repeated_elements_path=f"coverage/{args.name}_repeat_nodes.tab"
line_content=linecache.getline(repeated_elements_path,1)
if line_content:
    verbose_print("Adding repeated elements to the predictions...", end='\r')
    Path("walks/repeats").mkdir(parents=True, exist_ok=True)

    generate_repeat_paths(sample = args.name,
                          number_iterations = args.number_iterations,
                          filt_threshold = args.filt_gplas,
                          sd_coverage = args.repeats_coverage_sd)

    calculate_coocurrence_repeats(sample = args.name,
                                  number_iterations = args.number_iterations,
                                  pred_threshold = args.threshold_prediction,
                                  mod_threshold = args.modularity_threshold,
                                  sd_coverage = args.repeats_coverage_sd)

##3.6.2 If there are not repeated elementss, just rename the results files.
else:
    shutil.move("results/{args.name}_results_no_repeats.tab", "results/{args.name}_results.tab")
    shutil.move("results/{args.name}_bins_no_repeats.tab", "results/{args.name}_bins.tab")

#Check if output has been correctly created
if os.path.exists(f"results/{args.name}_results.tab") == False:
    #make this also an error_message() function??
    sys.exit("ERROR: Something went wrong while running gplas on the repeated elements")

verbose_print("Adding repeated elements to the predictions completed!")


##3.7 If the -k flag was not selected, delete intermediary files
if args.keep==False and args.mode!='extract':
    verbose_print("Intermediate files will be deleted. If you want to keep these files, use the -k flag")
    #improve use Path().unlink(missing_ok=True)?
    ##Delete files
    #Coverage files
    if Path(f"coverage/{args.name}_estimation.txt").exists():
        Path(f"coverage/{args.name}_clean_links.tab").unlink()
        Path(f"coverage/{args.name}_clean_prediction.tab").unlink()
        Path(f"coverage/{args.name}_clean_repeats.tab").unlink()
        Path(f"coverage/{args.name}_estimation.txt").unlink()
        Path(f"coverage/{args.name}_graph_contigs.tab").unlink()
        Path(f"coverage/{args.name}_initialize_nodes.tab").unlink()
        Path(f"coverage/{args.name}_isolated_nodes.tab").unlink()
        Path(f"coverage/{args.name}_repeat_nodes.tab").unlink()
        Path(f"coverage/{args.name}_repeats_graph.tab").unlink()
    #Walks normal mode
    if Path(f"walks/normal_mode/{args.name}_solutions.tab").exists():
        Path(f"walks/normal_mode/{args.name}_solutions.tab").unlink()
        #Path(f"walks/normal_mode/{args.name}_connections.tab").unlink()
    #Walks bold mode + unbinned solutions
    if Path(f"walks/bold_mode/{args.name}_solutions_bold.tab").exists():
        Path(f"walks/bold_mode/{args.name}_solutions_bold.tab").unlink()
        #Path(f"walks/bold_mode/{args.name}_connections_bold.tab").unlink()
        Path(f"walks/unbinned_nodes/{args.name}_solutions_unbinned.tab").unlink()
        Path(f"walks/{args.name}_solutions.tab").unlink()
    #Walks repeats
    if Path(f"walks/repeats/{args.name}_solutions.tab").exists():
        Path(f"walks/repeats/{args.name}_solutions.tab").unlink()
        #Path(f"walks/repeats/{args.name}_connections.tab").unlink()
    #Results no_repeats
    if Path(f"results/{args.name}_results_no_repeats.tab").exists():
        Path(f"results/{args.name}_results_no_repeats.tab").unlink()
        Path(f"results/{args.name}_bins_no_repeats.tab").unlink()
    
    ##Delete directories if they exist and are empty
    if Path("coverage/").exists() and not any(Path("coverage/").glob("*")):
        Path("coverage/").rmdir()

    if Path("walks/normal_mode/").exists() and not any(Path("walks/normal_mode/").glob("*")):
        Path("walks/normal_mode/").rmdir()

    if Path("walks/bold_mode/").exists() and not any(Path("walks/bold_mode/").glob("*")):
        Path("walks/bold_mode/").rmdir()

    if Path("walks/unbinned_nodes/").exists() and not any(Path("walks/unbinned_nodes/").glob("*")):
        Path("walks/unbinned_nodes/").rmdir()

    if Path("walks/repeats/").exists() and not any(Path("walks/repeats/").glob("*")):
        Path("walks/repeats/").rmdir()

    if Path("walks/").exists() and not any(Path("walks/").glob("*")):
        Path("walks/").rmdir()
    
##3.8 If there are no errors: Show success message and exit workflow
success_message()
sys.exit(0)

#TODO is this ever used?
def start():
    print("Starting gplas")
