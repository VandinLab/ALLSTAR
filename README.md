
# ALLSTAR: Inference of Reliable Causal Rules between Somatic Mutations and Cancer Phenotypes

Welcome to the ALLSTAR repository! 

This GitHub folder contains the code associated with the paper titled "ALLSTAR: Inference of Reliable Causal Rules between Somatic Mutations and Cancer Phenotypes". The code is designed to run on Python 3.10 and R, with wrappers provided as Shell scripts for convenience. We cannot guarantee full compatibility with Python 3.11 and above, as it may require additional dependencies and/or upgrades.

## Dependencies and Setup

**Virtual Environment** : Be sure to create and enter the virtual environment before launching ALLSTAR. You can do it by running these lines of code

```
python3.10 -m venv allstar_venv
source allstar_venv/bin/activate
```

**Linux Dependencies** : If you're using a Linux environment, ensure you have `python3.10`, `python3.10-venv`, `libopenmpi-dev`, `libpq-dev`, `graphviz-dev`, `libjpeg-turbo8-dev` installed. If your Linux distro is Ubuntu, you can install them with

```
sudo apt install -y python3.10 python3.10-venv libopenmpi-dev libpq-dev graphviz-dev libjpeg-turbo8-dev
```

**Python Dependencies** : Necessary Python dependencies are listed in `requirements.txt`. Be sure to install them with

```
pip install -r requirements.txt
```

**Replicating Experiments** : You can replicate all experiments by executing `wrapper.sh` in the root folder. Make the wrapper executable and launch it with

```
chmod +x wrapper.sh
./wrapper.sh
```

**Folder Structure** : The project is organized into multiple folders, each serving a specific purpose:

* `bash_scripts/`: Contains bash scripts used in our experiments, called by `wrapper.sh`.
* `data/`: Holds raw data used to generate input datasets.
* `images/`: Contains output images generated by our algorithm, as used in the paper.
* `inputs/`: Stores real-world input datasets used for code execution.
* `pathways/`: Holds pathway descriptions used for dataset parsing and gene merging.
* `py/`: Contains Python scripts used in our pipeline and experiments.
* `r/`: Holds R scripts used in our pipeline.
* `reports/`: Contains aggregated result reports.
* `results/`: Stores results of each run of the branch and bound algorithm, organized into subfolders named after the dataset.

# Running ALLSTAR

## Tips to Format the Input csv File

This will save you time (and pain), so read it carefully! 
ALLSTAR is built to natively accept csv (comma-separated) files, containing columns related to confounders, treatments and target, only in this rigorous order. The `--split/-s` parameter determines the index of the column representing the last confounder; the algorithm then assumes the last column as the target. 
Moreover, make sure to translate categorical variables into integers; for example:

| Confounder 1 | Treatment 1 | Treatment 2 | Target |
|----------|----------|----------|----------|
| 1 | 0 | 1 | 1 |
| 2 | 1 | 0 | 2 |
| 2 | 1 | 1 | 3 |
| 3 | 0 | 1 | 4 |
| 1 | 0 | 1 | 5 |

Remember, that ALLSTAR evaluates causal effects taking into consideration **all provided confounders**. If you have more than just one confounder, you could create more than one csv file, each one related to a specific confounder to evaliuate alone.

Therefore, make sure to build your database carefully!

## Parameters

You can then launch ALLSTAR using the `allstar.sh` Shell script. Ensure it's executable with `chmod +x allstar.sh`. ALLSTAR accepts several parameters, some with default values:

* `--file_path/-f`: Path to the input file.
* * `--y_value/-y`: Target value to evaluate.
* * `--top_k_rules/-k`: Number of top rules to output. Default: 5.
* `--graph_path/-g`: Path of the PPI graph file to use. Default: `PPI/FIsInGene_122921_with_annotations.txt`.
* `--max_rule_length/-l`: Maximum length of a rule. Default: 4.
* `--alpha/-a`: Confidence level alpha for FWER guarantees. Default: 0.05.
* `--threshold_manhattan/-t`: Manhattan distance threshold for rule generation. Default: 0.01.
* `--split/-s`: Column number referring to confounders; split+1 is the first treatment's column.
* * `--num_cores/-p`: Number of cores to use. Default: 1.

Options `-p`, `-k`, `-s`, `-l`, `-a`, and `-t` must be followed by numerical values, while `-f` and `-g` must be followed by the paths of the resources. We suggest editing your target variable so that it contains integers related to each target value

**Blacklist:** Inside the ALLSTAR folder you can find a `blacklist.csv` file. This can filled in with treatments/alterations that you may want to exclude from the analysis, without the need to edit manually the input dataset. Simply write one treatment per line, otherwise leave it blank. An example:

```
TP53_somatic
CDH1_LOH
BRCA1_meth
```

## Practical Example

The `input_file.csv` contains 5 confounders in the first 5 columns and ALLSTAR is set to investigate rules causally related to target equal to 1. The analysis will take up to 20 cores.

```
./allstar.sh -f input_file.csv -s 5 -y 1 -p 20
```

## License and Contact

This code is provided under the MIT license. For any inquiries, feel free to reach out to the authors at: [fabio.vandin@unipd.it](). We welcome your feedback and contributions!
