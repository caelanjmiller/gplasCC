import subprocess
from gplas.scripts.utils import quit_tool

def run_plasmidCC(infile, sample, minlen, species, custom_db_path):
        if species:
            cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -s {species} -l {minlen} -D -g -f"
        elif custom_db_path:
            cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -p {custom_db_path} -l {minlen} -D -g -f"

        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as err:
            print(f"plasmidCC returned non-zero exit status: {err.returncode}")
            quit_tool(err.returncode)


def print_speciesopts():
        cmd = "plasmidCC --speciesopts"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as err:
            print(f"plasmidCC returned non-zero exit status: {err.returncode}")
            quit_tool(err.returncode)
