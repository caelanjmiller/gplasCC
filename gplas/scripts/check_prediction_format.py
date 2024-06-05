import pandas as pd
from pandas.api.types import is_integer_dtype, is_float_dtype, is_object_dtype, is_string_dtype
from Bio.SeqIO.FastaIO import SimpleFastaParser


class PredictionFileFormatError(Exception):
    pass


def check_prediction(sample, path_prediction):
    #0. Load prediction file.
    prediction_file = pd.read_csv(path_prediction, sep='\t', header=0)

    #######################################################################################################################################

    #1. Check the number of columns
    if prediction_file.shape[1] != 5:
        raise PredictionFileFormatError("Prediction file should contain 5 tab-separated columns")

    #######################################################################################################################################

    #2. Check the column names
    if prediction_file.columns[0] != 'Prob_Chromosome':
        raise PredictionFileFormatError("First column should be named Prob_Chromosome (case sensitive)")

    if prediction_file.columns[1] != 'Prob_Plasmid':
        raise PredictionFileFormatError("Second column should be named Prob_Plasmid (case sensitive)")

    if prediction_file.columns[2] != 'Prediction':
        raise PredictionFileFormatError("Third column should be named Prediction (case sensitive)")

    if prediction_file.columns[3] != 'Contig_name':
        raise PredictionFileFormatError("Fourth column should be named Contig_name (case sensitive)")

    if prediction_file.columns[4] != 'Contig_length':
        raise PredictionFileFormatError("Fifth column should be named Contig_length (case sensitive)")

    #######################################################################################################################################

    #3. Check the data type of every column
    if(not is_float_dtype(prediction_file.dtypes['Prob_Chromosome']) and
       not is_integer_dtype(prediction_file.dtypes['Prob_Chromosome'])):
        raise PredictionFileFormatError("First column should be of type float/integer")

    if max(prediction_file['Prob_Chromosome']) > 1 or min(prediction_file['Prob_Chromosome']) < 0:
        raise PredictionFileFormatError("First column should contain values between 0 and 1")

    if(not is_float_dtype(prediction_file.dtypes['Prob_Plasmid']) and
       not is_integer_dtype(prediction_file.dtypes['Prob_Plasmid'])):
        raise PredictionFileFormatError("Second column should be of type float/integer")

    if max(prediction_file['Prob_Plasmid']) > 1 or min(prediction_file['Prob_Plasmid']) < 0:
        raise PredictionFileFormatError("Second column should contain values between 0 and 1")

    if(not is_object_dtype(prediction_file.dtypes['Prediction']) and
       not is_string_dtype(prediction_file.dtypes['Prediction'])):
        raise PredictionFileFormatError("Third column should be of type character/string")

    valid_predictions = ['Plasmid', 'Chromosome']
    if any([pred not in valid_predictions for pred in prediction_file['Prediction']]):
        raise PredictionFileFormatError(f"Third column values should be either {' or '.join(valid_predictions)} (case sensitive)")

    if(not is_object_dtype(prediction_file.dtypes['Contig_name']) and
       not is_string_dtype(prediction_file.dtypes['Contig_name'])):
        raise PredictionFileFormatError("Fourth column should be of type character/string")

    if not is_integer_dtype(prediction_file.dtypes['Contig_length']):
        raise PredictionFileFormatError("Fifth column should be of type integer")

    #######################################################################################################################################

    #4. check if plasmids exist in the prediction
    if 'Plasmid' not in prediction_file['Prediction'].values:
        raise PredictionFileFormatError("There are no plasmids in the prediction file, gplas can't do anything")

    #######################################################################################################################################

    #5. Check if the contig names in the prediction match the names in the FASTA file
    ##5.0 Get a path for fasta file.
    raw_nodes_path = f"gplas_input/{sample}_raw_nodes.fasta" 

    ##5.1 Get headers from fastafile.
    with open(raw_nodes_path) as file:
        fasta_headers = [str(node[0]) for node in SimpleFastaParser(file)]

    ##5.2 Get the headers from the prediction file
    prediction_headers = prediction_file['Contig_name']

    ##5.3 See if the predictions are in the fasta headers
    #TODO ASK Julian: do we also need to check the other way around? now there can be contigs in the FASTA that are not in the prediction file
    #TODO ASK Julian: should we use contigs.fasta instead of raw_nodes.fasta when checking this?
    comparison_output = [header in fasta_headers for header in prediction_headers]

    if not all(comparison_output):
        raise PredictionFileFormatError(f"Contig names in prediction file should match exactly with those in '{raw_nodes_path}'")

    #######################################################################################################################################

    #6. All checks are successful!
    return
