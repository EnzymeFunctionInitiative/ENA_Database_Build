# ENA_Database_Build
Code used to process the European Nucleotide Archive (ENA) database.

# Installation

```
# create the env to be used to run the workflow
conda create -n ena_db_build python=3.12
conda activate ena_db_build

# install necessary modules
conda install -y configparser gzip dask biopython
python3 -m pip install mysql-connector-python 
python3 -m pip install mysql-connector-python --upgrade
# or
conda install -y --file package-list.txt

# add this repo's code base as modules
python3 -m pip install .
```

# Running the Dask Workflow

Having followed the above installation instructions, the `dask_tskmgr.py` script is callable from anywhere and the submodule packages are importable. 
With this setup, the dask workflow can be tested on a local small compute resource like a laptop or desktop as well as on a HPC machine with more extensive compute resources.
The `dask_tskmgr.py` script has multiple input arguments that can be listed via: 

```

```

## Local Tests


## Full Scale Workflow



