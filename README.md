# ENA_Database_Build
Python workflow used to process the European Nucleotide Archive (ENA) dataset, gathering chromosome neighborhood context for genes that map to UniProtKB accession IDs.

# Installation
This workflow code expects an installation of a conda distribution, such as Miniconda. 
To install Miniconda, follow https://www.anaconda.com/docs/getting-started/miniconda/install instructions.
Once conda has been installed, follow the below instructions to create the `ena_db_build` environment that can be used to run the workflow.

```
# create the env to be used to run the workflow
conda create -n ena_db_build python=3.10
conda activate ena_db_build
conda config --add channels conda-forge

# install necessary modules
conda install -y configparser gzip dask biopython
python3 -m pip install mysql-connector-python 

# add this repo's code base as modules
python3 -m pip install .
```

# Running the Dask Workflow

Having followed the above installation instructions and with the `ena_db_build` environment active, the `ena_dask_tskmgr` command is executable on the command line from anywhere. 
Similarly, the submodules `dask_tasks`, `mysql_database`, and `parse_embl` can be imported within any interactive or scipted python code. 
This setup enables the ENA database build code to be implemented on a local small compute resource (for testing) as well as on an HPC machine with more extensive compute resources. 
As of 2025-02-11, the downloaded ENA dataset is ~20 TB, consisting of millions of relatively small gzip'd files; a large storage space and access to tens to hundreds of CPU processors are required to efficiently process all of the ENA dataset. 

The `ena_db_tskmgr` workflow has numerous input arguments.
These can be seen by running: 
```
(ena_db_build) $ ena_db_build -h
usage: ena_dask_tskmgr [-h] --db-config DB_CONFIG --db-name DB_NAME --ena-paths ENA_PATHS [ENA_PATHS ...] --output-dir OUTPUT_DIR [--scheduler-file SCHEDULER_FILE] [--n-workers N_WORKERS] [--tskmgr-log-file TSKMGR_LOG_FILE]
                       [--local-scratch LOCAL_SCRATCH]

Process the ENA Database

options:
  -h, --help            show this help message and exit
  --db-config DB_CONFIG, -conf DB_CONFIG
                        file path to the config file containing database connection information; assumed format is Windows INI.
  --db-name DB_NAME, -dbn DB_NAME
                        name of the EFI database to query for retrieving IDs.
  --ena-paths ENA_PATHS [ENA_PATHS ...]
                        arbitrary number of file paths that house subdirectories to be searched for ENA related dat.gz files.
  --output-dir OUTPUT_DIR, -out OUTPUT_DIR
                        path to the common output directory within which subdirectories and associated tab-separated data files will be saved.
  --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        path string to the dask scheduler file, default = '', indicating that a local dask cluster will be spun up rather than using a pre-defined scheduler and worker population.
  --n-workers N_WORKERS, -nWorkers N_WORKERS
                        number of workers available to perform tasks, default = 2.
  --tskmgr-log-file TSKMGR_LOG_FILE, -log TSKMGR_LOG_FILE
                        path string for a logging file, default = 'dask_tskmgr.log'.
  --local-scratch LOCAL_SCRATCH, -scratch LOCAL_SCRATCH
                        path string where temp files will be written, default = '', indicating do not write temp files to storage.
```

As the name suggests, this command runs a dask workflow using the `distributed` library to efficiently utilize available compute resources to perform a large number of tasks. 
See the "Dask Workflow Overview" section for more details about the task graph. 

## Full Scale Workflow
A sample SLURM batch script is provided in `~/batch_scripts/` that can be used to run the full scale workflow on the EFI HPC machine.
Also, the `env_setup.sh` script is provided within which all important environment variables called within the batch script are defined. 
Edit the environment script for your specific compute resource. 

A sample MySQL configuration file is provided with the necessary input parameters; update this file to match the access configuration for your MySQL server.

# Dask Workflow Overview

![Batch script flowchart of Dask scheduler, workers, and client working in tandem to perform the tasks defined in the dask_tskmgr.py script](https://github.com/EnzymeFunctionInitiative/ENA_Database_Build/blob/1622e7b1e6644f77ae7c7fc6a7a84e8811de4d9a/images/ena_db_build_flowchart.png)

