import plasmidCC
import subprocess
import shutil
import os

#TODO does not yet support the -d option for custom databases or -t for nr of threads
#TODO fix -o outdir and -p db_path
def run_plasmidCC(infile, sample, maxlen):
    if shutil.which("plasmidCC"):
        cmd = f"plasmidCC -i {infile} -o plasmidCC -n {sample} -s 'Escherichia coli' -l {maxlen} -D -g"
        try:
            return(subprocess.call(cmd, shell=True))
        except subprocess.CalledProcessError as e:
            print("plasmidCC has run into an error!")
            return(e.returncode)
    else:
        print("plasmidCC is not installed, please verify your installation.")

def cleanup_centrifuge(): #function to remove all unused centrifuge output
    return