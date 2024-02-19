#!/bin/bash

source bio_venv/bin/activate

# set default values
num_cores=1
top_k_rules=5
graph_path='PPI/FIsInGene_122921_with_annotations.txt'
blacklist='blacklist.csv'
max_rule_length=4
alpha=0.05
threshold_manhattan=0.01
y_value=0

# parse command-line arguments
while getopts ":p:f:b:s:k:y:g:l:a:t:" opt; do
  case $opt in
    p)
      num_cores=$OPTARG
      ;;
    f)
      file_path=$OPTARG
      ;;
    b)
      blacklist=$OPTARG # file containing blacklisted elements: one word per line
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

# clean up the input file according to blacklisted elements (i.e. genes, if any)
if [ ! -s "$blacklist" ]
then
  cp "$file_path" final_input_data.csv # no blacklisted words, so no need to modify the input file
else
  # read blacklisted words into an array
  mapfile -t blacklist_words < "$blacklist"

  # modify awk so that it removes the columns containing the blacklisted words from the input file
  awk_blacklist='{ for(i=1;i<=NF;i++) { skip=0;'
  for word in "${blacklist_words[@]}"; do
    awk_blacklist+=' if ($i ~ /\<'"$word"'\>/) { skip=1; break; }'
  done
  awk_blacklist+=' if (!skip) { printf "%s%s", $i, (i==NF ? RS : FS) } } }'

  # use the modified awk
  awk -F, "$awk_blacklist" "$file_path" > final_input_data.csv
fi

input_file_name='final_input_data.csv'

echo "Calculating. Raw results in $output_file_name"
# run ALLSTAR in calculus mode with the specified number of cores
mpiexec -n $num_cores python py/allstar.py \
    -f $input_file_name \
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
