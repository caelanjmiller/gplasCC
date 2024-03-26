import subprocess
import shutil

#TODO does not yet support the -d(-p) options for a custom database (path) or -t for nr of threads
#TODO change -o to be less annoying with nested directories; prob needs to change in plasmidCC sourcecode; or move files and delete directories
def run_plasmidCC(infile, sample, maxlen, species, custom_db_path):
    if shutil.which("plasmidCC"):
        if species:
            cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -s {species} -l {maxlen} -D -g -f"
        elif custom_db_path:
            cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -p {custom_db_path} -l {maxlen} -D -g -f"
        else:
            raise Exception("No input given for either species (-s) or custom_db_path (-p)")

        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e: #TODO improve this error handling; plasmidCC has its own error message already with quit_tool()
            print("\n")
            print("plasmidCC has run into an unexpected error!")
            raise Exception(e.returncode)
    else:
        raise Exception("plasmidCC is not installed, please verify your installation.")

def print_speciesopts():
    if shutil.which("plasmidCC"):
        cmd = "plasmidCC --speciesopts"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e: #TODO improve this error handling
            print("\n")
            print("plasmidCC has run into an unexpected error!")
            raise Exception(e.returncode)
    else:
        raise Exception("plasmidCC is not installed, please verify your installation.")
    