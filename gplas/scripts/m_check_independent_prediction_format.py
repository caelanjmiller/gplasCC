#!/usr/bin/env python3

import os
import sys
import pandas as pd
from pandas.api.types import is_integer_dtype, is_float_dtype, is_object_dtype, is_string_dtype
from Bio.SeqIO.FastaIO import SimpleFastaParser

def check_prediction(sample, path_prediction):
    
    #get a path for fasta file.
    #improve f string
    raw_nodes_path = f"gplas_input/{sample}_raw_nodes.fasta" 
    
    #Check if prediction file exists
    if(os.path.exists(path_prediction) == False):
        sys.exit("Prediction file does not exist or name is incorrect. Please, check your input for the -P argument.")
    
    #Check if fasta file exists
    if(os.path.exists(raw_nodes_path) == False):
        sys.exit("Nodes file does not exist or name is incorrect. The nodes file is produced by running the 'extract' command in gplas. Results will be located on the directory gplas_input/ and will be named as follows: $name_raw_nodes.fasta.")
    
    #load prediction file.
    prediction_file = pd.read_csv(path_prediction, sep="\t", header=0)
    
    #######################################################################################################################################
    
    #1. Check the number of columns
    if(prediction_file.shape[1] != 5):
        sys.exit("Error in prediction file format. The file should contain 5 columns, and they should be tab separated.")
    
    #######################################################################################################################################
    
    #2. Check the column names
    if(prediction_file.columns[0] != "Prob_Chromosome"):
        sys.exit("Error in prediction file format. First column should be named Prob_Chromosome (case sensitive).")
    
    if(prediction_file.columns[1] != "Prob_Plasmid"):
        sys.exit("Error in prediction file format. Second column should be named Prob_Plasmid (case sensitive).")
    
    if(prediction_file.columns[2] != "Prediction"):
        sys.exit("Error in prediction file format. Third column should be named Prediction (case sensitive).")
    
    if(prediction_file.columns[3] != "Contig_name"):
        sys.exit("Error in prediction file format. Fourth column should be named Contig_name (case sensitive).")
    
    if(prediction_file.columns[4] != "Contig_length"):
        sys.exit("Error in prediction file format. Fifth column should be named Contig_length (case sensitive).")
    
    #######################################################################################################################################
    
    #3. Check the data type of every column
    if((is_float_dtype(prediction_file.dtypes["Prob_Chromosome"]) == False) &
       (is_integer_dtype(prediction_file.dtypes["Prob_Chromosome"]) == False)):
        sys.exit("Error in prediction file format. First column should contain numeric values between 0 and 1.")
    
    if((is_float_dtype(prediction_file.dtypes["Prob_Plasmid"]) == False) &
       (is_integer_dtype(prediction_file.dtypes["Prob_Plasmid"]) == False)):
        sys.exit("Error in prediction file format. Second column should contain numeric values between 0 and 1.")
    
    if((is_object_dtype(prediction_file.dtypes["Prediction"]) == False) &
       (is_string_dtype(prediction_file.dtypes["Prediction"]) == False)):
        sys.exit("Error in prediction file format. Third column should contain data formatted as character indicating the type of prediciton (Plasmid or Chromosome). This is a case sensitive input.")
    
    if((is_object_dtype(prediction_file.dtypes["Contig_name"]) == False) &
       (is_string_dtype(prediction_file.dtypes["Contig_name"]) == False)):
        sys.exit("Error in prediction file format. Fourth column should contain data formatted as character indicating the contig names. Contig names should match exactly those provided in the FASTA file provided by the extract command.")
        
    if(is_integer_dtype(prediction_file.dtypes["Contig_length"]) == False):
        sys.exit("Error in prediction file format. Fifth column should contain integer values with lenght of contigs.")
        
    #######################################################################################################################################
    
    #4. check if plasmids exist in the prediction
    if("Plasmid" not in prediction_file["Prediction"].values):
        sys.exit("No plasmids are predicted, so gplas can't do anything")
    
    #######################################################################################################################################
    
    #5. Check if the names of the contigs in the predictions match the ones on the FASTA file
    ##5.1 get headers from fastafile.
    with open(raw_nodes_path) as file:
        fasta_headers = [str(node[0]) for node in SimpleFastaParser(file)]
    
    ##5.2 Get the headers from the prediction file
    prediction_headers = prediction_file["Contig_name"]
    
    ##5.3 See if the predictions are in the fasta headers
    #improve do we also need to check the other way around? now there can be contigs in the FASTA that are not in the prediction file
    comparison_output = [header in fasta_headers for header in prediction_headers]
    
    if(False in comparison_output):
        sys.exit("Error in contig names. Contig names in plasmid prediction file should match exactly with those in FASTA file obtained after runnning the 'extract' command")
    
    #######################################################################################################################################
    
    #6. All checks are successful
    print("Congrats! Your prediction file is correctly formatted. Now we are moving on to predicting your plasmids.")
    
    return