#!/bin/sh
awk '{{if($1 == "S") print ">"$1$2"_"$4"_"$5"\n"$3}}' test_ecoli.gfa 1>> gplas_input/ecoli_raw_nodes_unfiltered.fasta 2>> logs/ecoli_log_nodes.txt