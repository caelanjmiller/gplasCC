#!/bin/sh
awk -v min=500 'BEGIN {{RS = ">" ; ORS = ""}} length($2) >= min {{print ">"$0}}' gplas_input/ecoli_raw_nodes_unfiltered.fasta > gplas_input/ecoli_contigs.fasta