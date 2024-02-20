import argparse
import os
import sys

def file_exists(arg):
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(f"'{arg}' is not an existing file" + "\nPlease make sure the file exists and is spelled correctly")
    return arg

def is_valid_file(arg, extensions=['gfa']):
    if not os.path.isfile(arg):
        raise argparse.ArgumentTypeError(f"'{arg}' is not an existing file" + "\nPlease make sure the file exists and is spelled correctly")
    _, file_extension = os.path.splitext(arg)
    if not file_extension[1:].lower() in extensions:
        raise argparse.ArgumentTypeError(f"'{arg}' is not a file of type {' or '.join(extensions)}")
    return arg

#speciesopts needs to be manually updated if plasmidCC ever changes their species options
# can we do something like import plasmidCC.speciesopts (with a try except?)?
speciesopts = ['General',
               'Escherichia_coli',
               'Enterococcus_faecium',
               'Enterococcus_faecalis',
               'Salmonella_enterica',
               'Staphylococcus_aureus',
               'Acinetobacter_baumannii',
               'Klebsiella_pneumoniae']
def check_species(arg):
    if not arg in speciesopts:
        raise argparse.ArgumentTypeError(f"'{arg}' is not a recognised species" + "\nUse gplas with the --speciesopts flag for a list of all supported species")
    return arg

def check_output(path):
    if not os.path.exists(path):
        print("\n")
        print("Something went wrong while running gplas")
        print(f"Failed to create the following output: {path}")
        sys.exit(1)

def delete_file(file_path):
    if os.path.exists(file_path):
        os.remove(file_path)

def delete_empty_dir(dir_path):
    if os.path.exists(dir_path):
        if not any(os.listdir(dir_path)):
            os.rmdir(dir_path)
