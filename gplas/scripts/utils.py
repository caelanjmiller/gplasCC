import argparse
import os
import sys
from plasmidCC.scripts.utils import speciesopts


def quit_tool(exitcode=0):
    if exitcode != 0:
        print('\n', end='')
        print("This run of gplas has ended unexpectedly. Please check above for any error messages")
        sys.exit(-1)
    else:
        sys.exit(0)


def file_exists(arg):
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(f"'{arg}' is not an existing file" + "\nPlease make sure the file exists and is spelled correctly")
    return str(arg)


def is_valid_file(arg, extensions=['gfa']):
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(f"'{arg}' is not an existing file" + "\nPlease make sure the file exists and is spelled correctly")
    _, file_extension = os.path.splitext(arg)
    if not file_extension[1:].lower() in extensions:
        raise argparse.ArgumentTypeError(f"'{arg}' is not a file of type {' or '.join(extensions)}")
    return arg


def is_valid_dir(arg):
    if not os.path.isdir(arg):
        raise argparse.ArgumentTypeError(f"'{arg}' is not an existing directory, and I am afraid to create it")
    return arg


def check_species(arg):
    if not arg in speciesopts:
        raise argparse.ArgumentTypeError(f"'{arg}' is not a recognised species" + "\nUse gplas with the --speciesopts flag for a list of all supported species")
    return arg


def check_output(path):
    if not os.path.exists(path):
        print('\n')
        print(f"Failed to create the following output: {path}")
        quit_tool(-1)


def delete_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)


def delete_empty_dir(dir_path):
    if os.path.exists(dir_path):
        if not any(os.listdir(dir_path)):
            os.rmdir(dir_path)


def cleanup_centrifuge(sample):
    delete_file(f"plasmidCC/{sample}/{sample}_plasmids.fasta")
    delete_file(f"plasmidCC/{sample}/{sample}_centrifuge_classified.txt")


def cleanup_intermediary_files(sample):
    #Coverage files
    delete_file(f"coverage/{sample}_clean_links.tab")
    delete_file(f"coverage/{sample}_clean_prediction.tab")
    delete_file(f"coverage/{sample}_clean_repeats.tab")
    delete_file(f"coverage/{sample}_estimation.txt")
    delete_file(f"coverage/{sample}_graph_contigs.tab")
    delete_file(f"coverage/{sample}_initialize_nodes.tab")
    delete_file(f"coverage/{sample}_isolated_nodes.tab")
    delete_file(f"coverage/{sample}_repeat_nodes.tab")
    delete_file(f"coverage/{sample}_repeats_graph.tab")
    #Walks normal mode
    delete_file(f"walks/normal_mode/{sample}_solutions.tab")
    #Walks bold mode + unbinned solutions
    delete_file(f"walks/bold_mode/{sample}_solutions_bold.tab")
    delete_file(f"walks/unbinned_nodes/{sample}_solutions_unbinned.tab")
    delete_file(f"walks/{sample}_solutions.tab")
    #Walks repeats
    delete_file(f"walks/repeats/{sample}_solutions.tab")
    #Results no_repeats
    delete_file(f"results/{sample}_results_no_repeats.tab")
    delete_file(f"results/{sample}_bins_no_repeats.tab")
    #Centrifuge classification
    delete_file(f"plasmidCC/{sample}/{sample}_gplas.tab")
    #Delete directories if they exist and are empty
    delete_empty_dir("coverage/")
    delete_empty_dir("walks/normal_mode/")
    delete_empty_dir("walks/bold_mode/")
    delete_empty_dir("walks/unbinned_nodes/")
    delete_empty_dir("walks/repeats/")
    delete_empty_dir("walks/")
    delete_empty_dir(f"plasmidCC/{sample}/")
    delete_empty_dir("plasmidCC/")
