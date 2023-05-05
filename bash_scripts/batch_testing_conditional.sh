
if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_subtype.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_subtype.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_subtype.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_subtype.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_subtype.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_subtype.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_subtype.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_subtype.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi


if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_subtype.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_subtype.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_subtype.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_subtype.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_subtype.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_subtype.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_subtype.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_subtype.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi


if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_3N.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_3N.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_3N.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_3N.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_3N.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_3N.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_3N.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_3N.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_histo.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_histo.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_histo.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_histo.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_histo.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_histo.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_histo.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_histo.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_3N.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_3N.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_3N.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_3N.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_3N.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_3N.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_3N.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_3N.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_histo.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_histo.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_histo.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_histo.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_histo.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_histo.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_histo.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_histo.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi




if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_subtype.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_subtype.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_subtype.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_subtype.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_subtype.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_subtype.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_subtype.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_subtype.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_3N.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_3N.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_3N.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_3N.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_3N.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_3N.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_3N.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_3N.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_histo.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_histo.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_histo.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_histo.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_histo.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_histo.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_histo.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_histo.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_LVI.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_LVI.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_LVI.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/aggregated_mut_LVI.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_LVI.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_LVI.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/aggregated_mut_LVI.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/aggregated_mut_LVI.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi



if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_LVI.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_LVI.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_LVI.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LOH_meth_LVI.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_LVI.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_LVI.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_LOH_meth_LVI.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LOH_meth_LVI.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi


if [ "$1" = "v3" ]
then
# Version 3
mpiexec -n 30 python py/main.py -f inputs/som_LVI.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LVI.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LVI.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
python py/objective_parser.py -f inputs/som_LVI.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children
mpiexec -n 30 python py/main.py -f inputs/som_LVI.csv -s 7 -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LVI.csv -v 3 --max_rule_len 2 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
mpiexec -n 30 python py/main.py -f inputs/som_LVI.csv -s 7 -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
python py/objective_parser.py -f inputs/som_LVI.csv -v 3 --max_rule_len 3 --use_MHT_correction --use_graph_per_children --include_unmutated_treats
fi


