#!/usr/bin/env python3
if __name__ == "__main__":

    #from cProfile import run
    #import shutil
    #import linecache
    #import glob
    #import os
    #import sys
    #import argparse
    #import time
    #import subprocess
    #from pathlib import Path
    #from .version import version as VERSION
    ##os.chdir("C:/Users/oscar/Documenten/UU/04_BiBc/6.2_Research_Profile/gplas/gplas-2-python/gplas/scripts/module_development")
    ##os.getcwd()
    
    from m_check_independent_prediction_format import check_prediction
    from m_coverage import coverage
           
    ##arguments
    #only need sample, .gfa file, prediction file, and -flag parameters; all the rest is hardcoded like f"coverage/{sample}_clean_links"
    sample = "ecoli"
    path_assembly_graph = "test_ecoli.gfa"
    path_prediction = "test_ecoli_plasmid_prediction.tab"
    #Params
    classifier = "predict"
    threshold = float("0.95")

    ##workflow layout
    #0_mkdirs
    #1_awk_links
    #2_awk_nodes_extract_nodes
    #3_awk_nodes_filter
    #4_awk_nodes_rename
    
    check_prediction(sample=sample, path_prediction=path_prediction)

    coverage(sample=sample, path_prediction=path_prediction, classifier=classifier, threshold=threshold)

    
    #m_paths
    #m_paths_bold
    #m_coocurrence
    
    #5_extract_unbinned_nodes
    #6_combine_solutions
    
    #m_coocurrence_final
    #m_paths_repeats
    #m_coocurrence_repeats
    
    #remove_intermediate_files