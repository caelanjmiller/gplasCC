import shutil
import linecache
import os
import sys
import argparse
import glob
from .version import version as VERSION
#VERSION="1.0.0"
import time
from plasmidCC.scripts import utils as utilsCC
#import glob
#import subprocess
#from pathlib import Path

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
from gplas.scripts.m_run_plasmidcc import print_speciesopts
from gplas.scripts import m_utils as utils

start_time = time.time()

# Directories
pkgdir = os.path.dirname(__file__)

# Load ASCII logo
with open(f'{pkgdir}/figures/logo.txt', 'r') as file:
    read_logo = file.read()

#******************************#
#*                            *#
#* Command line parsing       *#
#*                            *#
#******************************#

class PriorityPrinting(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if option_string == '-h' or option_string == '--help':
            print(read_logo + '\n')
            parser.print_help()
        elif option_string == '-v' or option_string == '--version':
            print(f"gplas version {VERSION}")
        elif option_string == '--speciesopts':
            print_speciesopts()
        parser.exit()

#create a function to pass float ranges
parser = argparse.ArgumentParser(description="gplas: A tool for binning plasmid-predicted contigs into individual predictions", add_help=False)
parser.register('action', 'printing', PriorityPrinting)

inputgroup = parser.add_argument_group('General')
inputgroup.add_argument('-i', dest='input', type=utils.is_valid_file, required=True, help="Path to the graph file in GFA (.gfa) format, used to extract nodes and links")
#inputgroup.add_argument('-o', dest='outdir', type=utils.is_valid_dir, default=".", help="Output directory")  # TODO go through all code and append outdir where needed
inputgroup.add_argument('-n', dest='name', type=str, help="Name prefix for output files (default: input file name)")

classifiergroup = inputgroup.add_mutually_exclusive_group(required=True)
classifiergroup.add_argument('-s', dest='species', type=utils.check_species, help="Choose a species database for plasmidCC classification. Use --speciesopts for a list of all supported species")
classifiergroup.add_argument('-p', dest='custom_db_path', type=utilsCC.verify_user_db, help="Path to a custom Centrifuge database (name without file extensions)")
classifiergroup.add_argument('-P', dest='prediction', type=utils.file_exists, help="If not using plasmidCC. Provide a path to an independent binary classification file")
classifiergroup.add_argument('--extract', action='store_true', help="extract FASTA sequences from the assembly graph to use with an external classifier")

paramgroup = parser.add_argument_group('Parameters')
paramgroup.add_argument('-t', dest='threshold_prediction', type=float, default=0.5, help="Prediction threshold for plasmid-derived sequences (default: %(default)s)")
paramgroup.add_argument('-b', dest='bold_coverage_sd', type=int, default=5, help="Coverage variance allowed for bold walks to recover unbinned plasmid-predicted nodes (default: %(default)s)")
paramgroup.add_argument('-r', dest='repeats_coverage_sd', type=int, default=2, help="Coverage variance allowed for assigning repeats to bins (default: %(default)s)")
paramgroup.add_argument('-x', dest='number_iterations', type=int, default=20,help="Number of walk iterations per starting node (default: %(default)s)")
paramgroup.add_argument('-f', dest='filt_gplas', type=float, default=0.1, help="filtering threshold to reject outgoing edges (default: %(default)s)")
paramgroup.add_argument('-e', dest='edge_threshold', type=float, default=0.1, help="Edge threshold (default: %(default)s)")
paramgroup.add_argument('-q', dest='modularity_threshold', type=float, default=0.2, help="Modularity threshold to split components in the plasmidome network (default: %(default)s)")
paramgroup.add_argument('-l', dest='length_filter', type=int, default=1000, help="Filtering threshold for sequence length (default: %(default)s)")

othergroup = parser.add_argument_group('Other')
othergroup.add_argument('-k', '--keep', action='store_true', help="Keep intermediary files")
#othergroup.add_argument('--threads', type=int, default=1, help="Max number of threads to ")  #TODO add multi processing
othergroup.add_argument('--silent', action='store_true', help="Suppress verbose print statements")

infogroup = parser.add_argument_group('Info')
infogroup.add_argument('--speciesopts', action='printing', nargs=0, help="Prints a list of all supported species for the -s flag")
infogroup.add_argument('-v', '--version', action='printing', nargs=0, help="Prints gplas version")
infogroup.add_argument('-h', '--help', action='printing', nargs=0, help="Prints this message")
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
Your results are in gplas_input/{sample}_contigs.fasta. Please, use an external tool to classify the nodes in this file, and then bin them into individual plasmids using gplas.
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

# Check if the user specified a run name or to use file name data
infile = os.path.abspath(args.input)
infilename = os.path.basename(args.input)

if args.name:
    sample = args.name
else:
    sample, _ = os.path.splitext(infilename)

#Print messages
print('\n')
print(read_logo)
print('\n')

#Print chosen parameters
print("##################################################################")
print("Your results will be named:...........................", sample)
print("Input graph:..........................................", infilename)
print("Threshold for predicting plasmid-derived contigs:.....", args.threshold_prediction)
print("Number of plasmid walks created per node:.............", args.number_iterations)
print("Threshold of gplas scores:............................", args.filt_gplas)
print("Minimum frequency to consider an edge:................", args.edge_threshold)
print("Modularity threshold used to partition the network:...", args.modularity_threshold)
print("Coverage SD for bold mode:............................", args.bold_coverage_sd)
print("Coverage SD for repeats:..............................", args.repeats_coverage_sd)
print("Minimum sequence length:..............................", args.length_filter)
print("##################################################################" + '\n')

##_1.0 Run analysis
#_1.1 Extract nodes and links from the assembly graph
verbose_print("Extracting contigs from the assembly graph...", end='\r')
os.makedirs('gplas_input', exist_ok=True)

extract_nodes(sample, infile, args.length_filter)

utils.check_output(f"gplas_input/{sample}_raw_nodes.fasta")
verbose_print("Extracting contigs from the assembly graph completed!")

#_1.2 If in extract mode, exit workflow after succesful extraction. Else continue workflow
if args.extract:
    success_message_extract() #Exits the workflow

##_2.0 Obtain correct prediction file
#_2.1 Run plasmidCC if no independent prediction file was given
if args.species or args.custom_db_path:
    verbose_print("Running plasmidCC to generate prediction file..." + '\n')
    os.makedirs('plasmidCC', exist_ok=True)
    inputFASTA = f"gplas_input/{sample}_contigs.fasta"

    run_plasmidCC(inputFASTA, sample, args.length_filter, args.species, args.custom_db_path)
    utils.cleanup_centrifuge(sample)

    print('\n', end='')
    path_prediction = f"plasmidCC/{sample}/{sample}_gplas.tab"
    plasmidCC = True
else:
    path_prediction = args.prediction
    plasmidCC = False

#_2.2 Check if the prediction file is correctly formatted.
verbose_print("Checking prediction file...", end='\r')

check_prediction(sample, path_prediction, plasmidCC)

verbose_print("Valid prediction file found!")

##_3.0 Run gplas in normal mode
#_3.1 Extract nodes/links from the assembly graph
verbose_print("Calculating base coverage...", end='\r')
os.makedirs('coverage', exist_ok=True)

coverage(sample, path_prediction, args.threshold_prediction)

verbose_print("Calculating base coverage completed!")

#_3.2 Generate random walks
verbose_print("Generating random walks in normal mode...", end='\r')
os.makedirs("walks/normal_mode", exist_ok=True)

generate_paths(sample, args.number_iterations, args.filt_gplas, mode='normal')

verbose_print("Generating random walks in normal mode completed!")

#_3.3 Calculate coocurrence between walks
verbose_print("Calculating coocurrence of random walks...", end='\r')
os.makedirs("results/normal_mode", exist_ok=True)

calculate_coocurrence(sample, args.number_iterations, args.threshold_prediction, args.modularity_threshold, mode='normal')

utils.check_output(f"results/normal_mode/{sample}_results_no_repeats.tab")
verbose_print("Calculating coocurrence of random walks completed!")

##_4.0 Resolve unbinned contigs
#_4.1 Check for unbinned contigs
unbinned_path=f"results/normal_mode/{sample}_bin_Unbinned.fasta"
if os.path.exists(unbinned_path):
    #_4.1.1 Run gplas in bold mode if contigs were left unbinned
    verbose_print("Some contigs were left unbinned")  # improve tell user how many contigs are unbinned?
    #_4.1.1.1 Generate random walks
    verbose_print("Generating random walks in bold mode...", end='\r')
    os.makedirs("walks/bold_mode", exist_ok=True)

    generate_paths(sample, args.number_iterations, args.filt_gplas, args.bold_coverage_sd, mode='bold')

    verbose_print("Generating random walks in bold mode completed!")

    #_4.1.1.2 Extract unbinned solutions
    verbose_print("Extracting unbinned contigs from bold walks...", end='\r')
    os.makedirs("walks/unbinned_nodes", exist_ok=True)

    extract_unbinned_solutions(sample)

    verbose_print("Extracting unbinned contigs from bold walks completed!")

    #_4.1.1.3 Recalculate coocurrence of walks using the combined solutions
    verbose_print("Recalculating coocurrence of random walks...", end='\r')

    calculate_coocurrence(sample, args.number_iterations, args.threshold_prediction, args.modularity_threshold, mode='unbinned')

    verbose_print("Recalculating coocurrence of random walks completed!")

#_4.1.2 Copy files from normal mode if there were no unbinned contigs
else:
    for file in glob.glob(f"results/normal_mode/{sample}*"):
        shutil.copy(file, "results/")

utils.check_output(f"results/{sample}_results_no_repeats.tab")

##_5.0 Add repeated elements
#_5.1 Check for repeats
repeated_elements_path=f"coverage/{sample}_repeat_nodes.tab"
line_content=linecache.getline(repeated_elements_path,1)
if line_content:
    #_5.1.1 Run gplas on repeated elements
    verbose_print("Adding repeated elements to the predictions...", end='\r')
    os.makedirs("walks/repeats", exist_ok=True)

    #_5.1.1.1 Generate random walks
    generate_repeat_paths(sample, args.number_iterations, args.filt_gplas)

    #_5.1.1.2 Calculate coocurrence between walks
    calculate_coocurrence_repeats(sample, args.number_iterations, args.threshold_prediction, args.modularity_threshold, args.repeats_coverage_sd)

    verbose_print("Adding repeated elements to the predictions completed!")


#_5.1.2 If there are no repeated elements, just rename the results files.
else:
    shutil.move("results/{sample}_results_no_repeats.tab", "results/{sample}_results.tab")
    shutil.move("results/{sample}_bins_no_repeats.tab", "results/{sample}_bins.tab")

utils.check_output(f"results/{sample}_results.tab")

##_6.0 If the -k flag was not selected, delete intermediary files
if not args.keep:
    utils.cleanup_intermediary_files(sample)

##_7.0 Show success message and exit workflow
success_message()
