import subprocess
import shutil


def run_plasmidCC(infile, sample, minlen, species, custom_db_path):
    if shutil.which('plasmidCC'):
        if species:
            cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -s {species} -l {minlen} -D -g -f"
        elif custom_db_path:
            cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -p {custom_db_path} -l {minlen} -D -g -f"

        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as err:
            raise Exception(f"plasmidCC returned non-zero exit status: {err.returncode}")
    else:
        raise Exception("plasmidCC is not installed, please verify your installation")


def print_speciesopts():
    if shutil.which('plasmidCC'):
        cmd = "plasmidCC --speciesopts"
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as err:
            raise Exception(f"plasmidCC returned non-zero exit status: {err.returncode}")
    else:
        raise Exception("plasmidCC is not installed, please verify your installation")
