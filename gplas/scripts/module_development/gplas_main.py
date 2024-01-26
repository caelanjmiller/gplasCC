#!/usr/bin/env python3
if __name__ == "__main__":
    #import sys
    #import os
    #os.chdir("C:/Users/oscar/Documenten/UU/04_BiBc/6.2_Research_Profile/gplas/gplas-2-python/gplas/scripts/module_development")
    #os.getcwd()
    
    from m_check_independent_prediction_format import check_prediction
       
    ##arguments
    sample_name = "ecoli"
    path_prediction = "test_ecoli_plasmid_prediction.tab"
    
    
    
    ##workflow layout
    #0_mkdirs
    #1_awk_links
    #2_awk_nodes_extract_nodes
    #3_awk_nodes_filter
    #4_awk_nodes_rename
    
    check_prediction(sample_name, path_prediction)
    
    #m_coverage
    #m_paths
    #m_paths_bold
    #m_coocurrence
    
    #5_extract_unbinned_nodes
    #6_combine_solutions
    
    #m_coocurrence_final
    #m_paths_repeats
    #m_coocurrence_repeats
    
    #remove_intermediate_files