#!/bin/sh
touch ./walks/unbinned_nodes/ecoli_solutions_unbinned.tab; for node in $(grep Unbinned ./results/normal_mode/ecoli_results_no_repeats.tab | cut -f 5 -d ' '); do \
grep -w "^${node}" ./walks/bold_mode/ecoli_solutions_bold.tab >> ./walks/unbinned_nodes/ecoli_solutions_unbinned.tab || continue; \
done