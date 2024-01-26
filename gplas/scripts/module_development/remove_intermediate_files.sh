#!/bin/bash

while getopts ":n:" opt; do
    case $opt in
        n) name=$OPTARG;;
    esac
done

#remove files
if [[ -f walks/"$name"_solutions.csv ]]; then
rm walks/"$name"_solutions*
fi

rm walks/normal_mode/"$name"_solutions*
rm walks/normal_mode/"$name"_connections*

if [[ -f walks/bold_mode/"$name"_solutions_bold.csv ]]; then
rm walks/bold_mode/"$name"_connections*
rm walks/bold_mode/"$name"_solutions*
rm walks/unbinned_nodes/"$name"_solutions*
fi

if [[ -f walks/repeats/"$name"_solutions.csv ]]; then
rm walks/repeats/"$name"_connections*
rm walks/repeats/"$name"_solutions*
fi

rm coverage/"$name"_clean*
rm coverage/"$name"_graph*
rm coverage/"$name"_estimation*
rm coverage/"$name"_repeat*
rm coverage/"$name"_initialize*
rm coverage/"$name"_isolated*

if [[ -f results/normal_mode/"$name"_results.tab ]]; then
rm results/normal_mode/"$name"_bin*
rm results/normal_mode/"$name"_results*
rm results/normal_mode/"$name"*png
fi

rm templates/"$name"_assembly*

if [[ -f results/"$name"_results_no_repeats.tab ]]; then
rm results/"$name"_*_no_repeats.tab
fi
