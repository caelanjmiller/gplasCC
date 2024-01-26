#!/usr/bin/env python3
if __name__ == "__main__":
    #import sys
    #import os
    #os.chdir("C:/Users/oscar/Documenten/UU/04_BiBc/6.2_Research_Profile/gplas/gplas-2-python/gplas/scripts/module_development")
    #os.getcwd()
    
    from m_check_independent_prediction_format import check_prediction
    from m_coverage import coverage
       
    ##arguments
    sample_name = "ecoli"
    #Inputs    
    path_nodes = "gplas_input/ecoli_raw_nodes.fasta"
    path_links = "gplas_input/ecoli_raw_links.txt"
    path_prediction = "test_ecoli_plasmid_prediction.tab"
        
    #Params
    classifier = "predict"
    threshold = float("0.95")
    
    #Outputs
    output_graph_contigs = "coverage/ecoli_graph_contigs.tab"
    output_clean_links = "coverage/ecoli_clean_links.tab"
    output_graph_repeats = "coverage/ecoli_repeats_graph.tab"
    output_clean_prediction = "coverage/ecoli_clean_prediction.tab"
    output_isolated_nodes = "coverage/ecoli_isolated_nodes.tab"
    output_clean_repeats = "coverage/ecoli_clean_repeats.tab"
    output_initialize_nodes = "coverage/ecoli_initialize_nodes.tab"
    output_repeat_nodes = "coverage/ecoli_repeat_nodes.tab"
    output_cov_estimate = "coverage/ecoli_estimation.txt"
    
    ##workflow layout
    #0_mkdirs
    #1_awk_links
    #2_awk_nodes_extract_nodes
    #3_awk_nodes_filter
    #4_awk_nodes_rename
    
    check_prediction(sample_name, path_prediction)

    coverage(path_nodes, path_links, path_prediction, classifier, threshold, output_graph_contigs, output_clean_links, output_graph_repeats, output_clean_prediction, output_isolated_nodes, output_clean_repeats, output_initialize_nodes, output_repeat_nodes, output_cov_estimate)

    
    #m_paths
    #m_paths_bold
    #m_coocurrence
    
    #5_extract_unbinned_nodes
    #6_combine_solutions
    
    #m_coocurrence_final
    #m_paths_repeats
    #m_coocurrence_repeats
    
    #remove_intermediate_files