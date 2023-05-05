# ALLSTAR
Inference of ReliAble CausaL RuLes between Somatic MuTAtions and CanceR Phenotypes

This Github folder contains the code associated to the paper "ALLSTAR: Inference of ReliAble CausaL RuLes between Somatic MuTAtions and CanceR Phenotypes". 
The code runs on python 3 and R. Wrappers are written as Shell scripts.

- Linux dependencies
If you are running a Linux environment, you may need the following dependencies: libjpeg-turbo8-dev

- How to replicate the synthetic experiments
It is possible to replicate all the experiments by running wrapper.sh on the root folder.
You may need to make the wrapper executable by using chmod +x wrapper.sh. 
This script will automatically create a python environment with the required packages, install those required for R, parse the datasets to save them in the appropriate format, run the experiments code, parse the results and finally create the plots.

- Folder structure
The project is diveded in multiple folders, each one created by the wrapper method if not already present in this git folder. Here will follows a list of some of them with the respective explanation:
  * bash_scripts/ contains some of the bash scripts we used to run our experiments. Those are in term called by wrapped.sh (on the main folder) that execute the whole code
  * data/ contains the raw data that we used to generate our input datasets
  * images/ contains the images in output by our algorithm, that are the ones used in the paper as well
  * inputs/ contains the input datasets used to run the code
  * pathways/ contains the pathways description that we used to parse the datasets and merge mutated genes together
  * py/ contains the python scripts we used in our pipeline and in our experiments
  * r/ contains the R scripts we used in our pipeline
  * reports/ contains the reports for aggregated results of multiple runs, in order to inspect the agglomerative results
  * results/ contains the results of each run of the branch and bound alogrithm, each one divided in a subfolder renamed as the dataset name

- Running ALLSTAR
ALLSTAR can be launched via the omonymous Shell script allstar.sh; you may need to make it executable by chmod +x allstar.sh.
ALLSTAR takes the following parameters, some of which with default values:
  * num_cores: The number of cores to use. Default: 1. To set this parameter to a different value, use the -p option followed by the desired value, e.g., allstar.sh -p 4.
  * top_k_rules: The number of top rules to be output. Default: 5. To set this parameter to a different value, use the -k option followed by the desired value, e.g., allstar.sh -k 5.
  * graph_path: The path of the PPI graph file to be used. It must have Gene1 and Gene2 as first two columns. Default: PPI/FIsInGene_122921_with_annotations.txt. To set this parameter to a different value, use the -g option followed by the path to the desired file, e.g., allstar.sh -g path/to/graph/file.txt.
  * max_rule_length: The maximum length of a rule. Default: 4. To set this parameter to a different value, use the -l option followed by the desired value, e.g., allstar.sh -l 4.
  * alpha: The confidence level alpha to provide FWER guarantees. Default: 0.05. To set this parameter to a different value, use the -a option followed by the desired value, e.g., allstar.sh -a 0.01.
  * threshold_manhattan: The Manhattan distance threshold to be used in the rule generation. Default: 0.01. To set this parameter to a different value, use the -t option followed by the desired value, e.g., allstar.sh -t 0.01.
  * file_path: The input file's path.
  * split: The last column number referring to confounders: split+1 is the first treatment's column.
  * y_value: The target value to be evaluated.
Note that the options -p, -k, -g, -l, -a, and -t must be followed by a numerical value, while the options -f, -s, and -y must be followed by a string.

The code comes under a MIT license. Feel free to contact the authors at:
fabio.vandin@unipd.it
