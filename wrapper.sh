#!/bin/bash

## Creation of the environment
python3.10 -m venv bio_venv
source bio_venv/bin/activate

## Requirement installation
pip install -r requirements.txt

## Data preparation
mkdir -p inputs
mkdir -p inputs/synthetic
mkdir -p inputs/synthetic/g1
mkdir -p inputs/synthetic/g2
mkdir -p inputs/synthetic/g3
Rscript r/dataprep_som_LOH_meth_various_targets.R

mkdir -p results
mkdir -p reports
python py/fisher_exact.py > results/brute_force_ranks_fisher_exact.log
python py/analysis_fisher_comparison.py  # saves the comparison between effect and standard techniquest in images/fisher_ratios_effects_comparison.png

## Experimental results w/ synthetic data
python py/synthetic_ds_creator.py -v 1 -r 10 -e 100
python py/synthetic_ds_creator.py -v 2 -r 10 -e 15
python py/synthetic_ds_creator.py -v 3 -r 10 -e 100
sh bash_scripts/G1_synthetic_experiments.sh
sh bash_scripts/G2_synthetic_experiments.sh
sh bash_scripts/G3_synthetic_experiments.sh
python py/parser_synthetics.py -v 1 -i results/synth_v11/v11 -e 100 --naive
python py/parser_synthetics.py -v 1 -i results/synth_v11/v11 -e 100 --reliable
python py/parser_synthetics.py -v 1 -i results/synth_v11/v11 -e 100 --corrected
python py/parser_synthetics.py -v 2 -i results/synth_v12/v12 -e 15 --graph
python py/parser_synthetics.py -v 2 -i results/synth_v12/v12 -e 15 --no_graph
python py/parser_synthetics.py -v 3 -i results/synth_v13/v13 -e 100 --light
python py/parser_synthetics.py -v 3 -i results/synth_v13/v13 -e 100 --strong
python py/synth_exp_res_verification.py > results/synthetic_exp_result_verification.log # saves the print for counting the False Positives in the log, and the images in images/folder

# ## Experimental results w/ real data
sh bash_scripts/batch_testing.sh
