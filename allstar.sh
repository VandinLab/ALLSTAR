#!/bin/bash

# set default values
num_cores=1
top_k_rules=5
graph_path='PPI/FIsInGene_122921_with_annotations.txt'
max_rule_length=4
alpha=0.05
threshold_manhattan=0.01
y_value=0

# parse command-line arguments
while getopts ":p:f:s:k:y:g:l:a:t:" opt; do
  case $opt in
    p)
      num_cores=$OPTARG
      ;;
    f)
      file_path=$OPTARG
      ;;
    s)
      split=$OPTARG
      ;;
    k)
      top_k_rules=$OPTARG
      ;;
    y)
      y_value=$OPTARG
      ;;
    g)
      graph_path=$OPTARG
      ;;
    l)
      max_rule_length=$OPTARG
      ;;
    a)
      alpha=$OPTARG
      ;;
    t)
      threshold_manhattan=$OPTARG
      ;;
  esac
done

# get current timestamp
timestamp=$(date +"%Y-%m-%d_%H-%M-%S")

# Get the filename using basename
filename=$(basename "$file_path")

# set output file name
output_file_name="results/ALLSTAR_runs/ALLSTAR_${filename}_${timestamp}.txt"

echo "Calculating. Raw results in $output_file_name"
# run ALLSTAR in calculus mode with the specified number of cores
mpiexec -n $num_cores python py/allstar.py \
    -f $file_path \
    -s $split \
    -k $top_k_rules \
    -y $y_value \
    -g $graph_path \
    -l $max_rule_length \
    -a $alpha \
    -t $threshold_manhattan \
    -o $output_file_name \
    --calculus_mode
echo "Parsing results."
python py/allstar.py \
    -o $output_file_name \
    --parse_mode
echo "Analysis completed."