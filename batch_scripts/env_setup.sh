
EBROOTENA=/path/to/ENA/root/directory
EBVERSIONENA=ena_2025_02	# just an example string

CONDA_HOME=/path/to/conda/installation

# NEED TO  INITIALIZE THE CONDA ENVIRONMENT SO THE ena_db_build ENV IS IN PATH
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('$CONDA_HOME/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "$CONDA_HOME/etc/profile.d/conda.sh" ]; then
        . "$CONDA_HOME/etc/profile.d/conda.sh"
    else
        export PATH="$CONDA_HOME/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

WORKING_DIR=$PWD
OUTPUT_DIR=$PWD/TEST
SCRATCH_DIR=/scratch/$EBVERSIONENA
DB_CONFIG=/path/to/mysql/webserver/config/file
DB_NAME=efi_202412	# or the equivalent version string for the relevant MySQL EFI Database to be queried
N_WORKERS=$((SLURMNTASKS - 2)) # scheduler and client script both need at least one cpu each.

# dask architecture settings
# see https://docs.dask.org/en/latest/deploying-cli.html
scheduler_file=$WORKING_DIR/scheduler_file.json
dask_interface="eth0"
nThreads_per_worker=1
dask_scheduler_arguments="--no-dashboard --no-jupyter --no-show --scheduler-file ${SCHEDULER_FILE} --interface $dask_interface"
dask_worker_arguments="--no-dashboard --no-nanny --reconnect --nthreads $nThreads_per_worker --nworkers 1 --interface $dask_interface --scheduler-file $SCHEDULER_FILE"


