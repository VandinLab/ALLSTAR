input_dir="inputs/synthetic/g1"

for file in $input_dir/*.csv; do  
    mpiexec -n 30 python py/main.py -f "$file" -s 0 -v 11 --max_rule_len 3
done