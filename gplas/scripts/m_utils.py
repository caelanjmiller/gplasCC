import argparse
import os
import sys

#def quit_tool(exitcode):
#    print("\n")
#    print("gplas has run into a problem:")
#    print(exitcode)
#    sys.exit(1)

#def is_valid_dir(arg):
#    if not os.path.isdir(arg):
#        raise argparse.ArgumentTypeError(f"'{arg}' is not an existing directory" + "\nPlease make sure the directory exists and is spelled correctly")
#    return(arg)

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
speciesopts = ['General','Escherichia coli','Enterococcus faecium','Enterococcus faecalis','Salmonella enterica','Staphylococcus aureus','Acinetobacter baumannii','Klebsiella pneumoniae']
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
