#!/bin/sh
awk -F "\\t" '{{if($1 == "L") print $N}}' test_ecoli.gfa 1>> gplas_input/ecoli_raw_links.txt 2>> logs/ecoli_log_links.txt