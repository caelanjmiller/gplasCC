import os
import pandas as pd
from pandas.api.types import is_integer_dtype, is_float_dtype, is_object_dtype, is_string_dtype
from Bio.SeqIO.FastaIO import SimpleFastaParser
from m_utils import quit_tool


def format_error(message):
    print('\n')
    print("Error in prediction file format:")
    print(message)
    quit_tool(-1)


def check_prediction(sample, path_prediction, plasmidCC):
    #get a path for fasta file.
    raw_nodes_path = f"gplas_input/{sample}_raw_nodes.fasta" 

    #Check if prediction file exists
    if not os.path.exists(path_prediction):
        print('\n')
        print("Prediction file does not exist or name is incorrect")
        if plasmidCC:
            print("Please check if something went wrong with your plasmidCC run")
        else:
            print("Please check your input for the -P argument")
        quit_tool(-1)

    #Check if fasta file exists
    if not os.path.exists(raw_nodes_path):
        print("Nodes file does not exist or name is incorrect")
        print(f"Please verify that '{raw_nodes_path}' exsists")
        quit_tool(-1)

    #load prediction file.
    prediction_file = pd.read_csv(path_prediction, sep='\t', header=0)

    #######################################################################################################################################

    #1. Check the number of columns
    if prediction_file.shape[1] != 5:
        format_error("The file should contain 5 columns, and they should be tab separated")

    #######################################################################################################################################

    #2. Check the column names
    if prediction_file.columns[0] != 'Prob_Chromosome':
        format_error("First column should be named Prob_Chromosome (case sensitive)")

    if prediction_file.columns[1] != 'Prob_Plasmid':
        format_error("Second column should be named Prob_Plasmid (case sensitive)")

    if prediction_file.columns[2] != 'Prediction':
        format_error("Third column should be named Prediction (case sensitive)")

    if prediction_file.columns[3] != 'Contig_name':
        format_error("Fourth column should be named Contig_name (case sensitive)")

    if prediction_file.columns[4] != 'Contig_length':
        format_error("Fifth column should be named Contig_length (case sensitive)")

    #######################################################################################################################################

    #3. Check the data type of every column
    if(not is_float_dtype(prediction_file.dtypes['Prob_Chromosome']) and
       not is_integer_dtype(prediction_file.dtypes['Prob_Chromosome'])):
        format_error("First column should contain numeric values")

    if max(prediction_file['Prob_Chromosome']) > 1 or min(prediction_file['Prob_Chromosome']) < 0:
        format_error("First column should contain values between 0 and 1")

    if(not is_float_dtype(prediction_file.dtypes['Prob_Plasmid']) and
       not is_integer_dtype(prediction_file.dtypes['Prob_Plasmid'])):
        format_error("Second column should contain numeric values")

    if max(prediction_file['Prob_Plasmid']) > 1 or min(prediction_file['Prob_Plasmid']) < 0:
        format_error("Second column should contain values between 0 and 1")

    if(not is_object_dtype(prediction_file.dtypes['Prediction']) and
       not is_string_dtype(prediction_file.dtypes['Prediction'])):
        format_error("Third column should contain data formatted as character")

    valid_predictions = ['Plasmid', 'Chromosome']
    if any([pred not in valid_predictions for pred in prediction_file['Prediction']]):
        format_error(f"Third column values should be either {' or '.join(valid_predictions)} (case sensitive)")

    if(not is_object_dtype(prediction_file.dtypes['Contig_name']) and
       not is_string_dtype(prediction_file.dtypes['Contig_name'])):
        format_error("Fourth column should contain data formatted as character")

    if not is_integer_dtype(prediction_file.dtypes['Contig_length']):
        format_error("Fifth column should contain integer values")

    #######################################################################################################################################

    #4. check if plasmids exist in the prediction
    if 'Plasmid' not in prediction_file['Prediction'].values:
        format_error("There are no plasmids in the prediction file, gplas can't do anything")

    #######################################################################################################################################

    #5. Check if the names of the contigs in the predictions match the ones on the FASTA file
    ##5.1 get headers from fastafile.
    with open(raw_nodes_path) as file:
        fasta_headers = [str(node[0]) for node in SimpleFastaParser(file)]

    ##5.2 Get the headers from the prediction file
    prediction_headers = prediction_file['Contig_name']

    ##5.3 See if the predictions are in the fasta headers
    #improve do we also need to check the other way around? now there can be contigs in the FASTA that are not in the prediction file
    ### use the contig file instead of the raw nodes? ask Julian about it
    comparison_output = [header in fasta_headers for header in prediction_headers]

    if not all(comparison_output):
        format_error(f"Contig names in prediction file should match exactly with those in '{raw_nodes_path}'")

    #######################################################################################################################################

    #6. All checks are successful!
    return
