#!/bin/bash
#SBATCH --partition=efi
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --job-name="ena_db_build"
#SBATCH --nodelist=compute-9-14

date

echo $EBROOTENA		# global path to the root dir of the ENA dataset
echo $EBVERSIONENA	# version string of the ENA dataset

# NEED TO  INITIALIZE THE CONDA ENVIRONMENT SO THE ena_db_build ENV IS IN PATH
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/path/to/conda/installation/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/path/to/conda/installation/etc/profile.d/conda.sh" ]; then
        . "/path/to/conda/installation/etc/profile.d/conda.sh"
    else
        export PATH="/path/to/conda/installation/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate ena_db_build

working_dir=$PWD
scratch_dir=/scratch/$EBVERSIONENA
output_dir=$PWD/TEST
DB_CONFIG=/path/to/mysql/webserver/config/file
DB_NAME=efi_202412	# or the equivalent version string for the relevant MySQL EFI Database to be queried

################################################################################

# empty string to gather process IDs associated with the dask scheduler and worker(s) calls
dask_pids=""

echo "Spinning up the Scheduler"
SCHEDULER_FILE=$working_dir/scheduler_file.json
dask scheduler --no-dashboard --no-jupyter --no-show --scheduler-file ${SCHEDULER_FILE} --interface 'eth0' > dask_scheduler.out 2>&1 &
dask_pids="$dask_pids $!"

# spin up nWorkers workers, each with 1 core
echo "Spinning up the Workers"
srun -n 62 \
	dask worker --no-dashboard --no-nanny --reconnect --nthreads 1 --nworkers 1 --interface 'eth0' --scheduler-file ${SCHEDULER_FILE} > dask_worker.out 2>&1 &
dask_pids="$dask_pids $!"

# run the workflow on all ENA files in $EBROOTENA
echo "Starting Client Script"
python3 ena_dask_tskmgr --db-config $DB_CONFIG --db-name $DB_NAME --ena-paths $EBROOTENA/sequence $EBROOTENA/wgs/public  $EBROOTENA/wgs/suppressed  --output-dir $output_dir --local-scratch $scratch_dir --scheduler-file ${SCHEDULER_FILE} --n-workers 62 --tskmgr-log-file $working_dir/tskmgr.log

for pid in $dask_pids
do
        kill $pid
done

date

