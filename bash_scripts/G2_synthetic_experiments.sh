input_dir="inputs/synthetic/g2"
STRING="_e15"

for r in 3 4 5 6 7; do
    for file in $input_dir/*"$STRING"*.csv; do
        mpiexec -n 30 python py/main.py -f "$file" -s 0 -v 12 --max_rule_len $r
        mpiexec -n 30 python py/main.py -f "$file" -s 0 -v 12 --max_rule_len $r --use_graph_per_children
    done
done