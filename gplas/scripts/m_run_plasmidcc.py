import subprocess
import shutil

#TODO does not yet support the -d(-p) options for a custom database (path) or -t for nr of threads
#TODO change -o to be less annoying with nested directories; prob needs to change in plasmidCC sourcecode; or move files and delete directories
def run_plasmidCC(infile, sample, species, maxlen):
    if shutil.which("plasmidCC"):
        cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -s \'{species}\' -l {maxlen} -D -g"
        try:
            subprocess.call(cmd, shell=True)
        except subprocess.CalledProcessError as e: #TODO improve this error handling
            print("\n")
            print("plasmidCC has run into an unexpected error!")
            raise Exception(e.returncode)
    else:
        raise Exception("plasmidCC is not installed, please verify your installation.")

def print_speciesopts():
    if shutil.which("plasmidCC"):
        cmd = "plasmidCC --speciesopts"
        try:
            subprocess.call(cmd, shell=True)
        except subprocess.CalledProcessError as e: #TODO improve this error handling
            print("\n")
            print("plasmidCC has run into an unexpected error!")
            raise Exception(e.returncode)
    else:
        raise Exception("plasmidCC is not installed, please verify your installation.")
    