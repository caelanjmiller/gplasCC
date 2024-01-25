#!/usr/bin/Rscript
print("Running the script")
suppressMessages(library(Biostrings))

arguments = commandArgs(trailingOnly=TRUE)
isolate_name <- arguments[1]
prediction_dir <- arguments[2]

#get a path for the prediction file
prediction_file_path<-prediction_dir
#get a path for fasta file.
raw_nodes_path<-paste('gplas_input/',isolate_name,'_raw_nodes.fasta',sep='')

#Check if prediction file exists
if (file.exists(prediction_file_path)==FALSE) {
    stop ('Prediction file does not exist or name is incorrect. Please, check your input for the -P argument.')
}

#Check if fasta file exists
if (file.exists(raw_nodes_path)==FALSE) {
    stop ('Nodes file does not exist or name is incorrect. The nodes file is produced by running the "extract" command in gplas. Results will be located on the directory gplas_input/ and will be named as follows: $name_raw_nodes.fasta.' )
}

#load prediction file.
prediction_file<-read.csv(prediction_file_path, sep='\t', header=TRUE)

#1. Check the number of columns
if (ncol(prediction_file)!=5) {
  stop('Error in prediction file format. The file should contain 5 columns, and they should be tab separated.')
}

#######################################################################################################################################


#2. Check the column names
if (colnames(prediction_file)[1]!='Prob_Chromosome') {
  stop("Error in prediction file format. First column should be named Prob_Chromosome (case sensitive).")
}

if (colnames(prediction_file)[2]!='Prob_Plasmid') {
  stop("Error in prediction file format. Second column should be named Prob_Plasmid (case sensitive).")
}

if (colnames(prediction_file)[3]!='Prediction') {
  stop("Error in prediction file format. Third column should be named Prediction (case sensitive).")
}

if (colnames(prediction_file)[4]!='Contig_name') {
  stop("Error in prediction file format. Fourth column should be named Contig_name (case sensitive).")
}

if (colnames(prediction_file)[5]!='Contig_length') {
  stop("Error in prediction file format. Fifth column should be named Contig_length (case sensitive).")
}

###################################################################################################################

#3. Check the data type of every column
if (class(prediction_file[,c(1)])!='numeric' && class(prediction_file[,c(1)])!='integer'){
  stop("Error in prediction file format. First column should contain numeric values between 0 and 1.")
}

if (class(prediction_file[,c(2)])!='numeric' && class(prediction_file[,c(2)])!='integer'){
  stop("Error in prediction file format. Second column should contain numeric values between 0 and 1.")
}

if (class(prediction_file[,c(3)])!='character' && class(prediction_file[,c(3)])!='factor'){
  stop("Error in prediction file format. Third column should contain data formatted as character indicating the type of prediciton (Plasmid or Chromosome). This is a case sensitive input.")
}

if (class(prediction_file[,c(4)])!='character' && class(prediction_file[,c(4)])!='factor'){
  stop("Error in prediction file format. Fourth column should contain data formatted as character indicating the contig names. Contig names should match exactly those provided in the FASTA file provided by the extract command.")
}

if (class(prediction_file[,c(5)])!='integer'){
  stop("Error in prediction file format. Fifth column should contain integer values with lenght of contigs")
}

#4. check if plasmids exist in the prediction
if (! "Plasmid" %in% unique(prediction_file[,3])){
  stop("No plasmids are predicted, so gplas can't do anything")
}
#####################################################################################################################

#4. Check if the names of the contigs in the predictions match the ones on the FASTA file
##4.1 get headers from fastafile.
fasta_file<-readDNAStringSet(filepath = raw_nodes_path, format="fasta")
fasta_headers<-as.vector(names(fasta_file))

##4.2 Get the headers from the prediction file
prediction_headers<-as.vector(prediction_file[,c(4)])
##4.3 See if the predictions are in the fasta headers
comparision_output<-prediction_headers %in% fasta_headers

if (FALSE %in% comparision_output) {
  stop("Error in contig names. Contig names in plasmid prediction file should match exactly with those in FASTA file obtained after runnning the 'extract' command")
}

print("Congrats! Your prediction file is correctly formatted. Now we are moving on to predicting your plasmids.")
