#!/bin/bash
#SBATCH --partition=efi
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --job-name="ena_db_build_step1"
#SBATCH --nodelist=compute-9-14

date

module load ena/20241205 

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/n-z/rbdavid/Apps/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/n-z/rbdavid/Apps/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/n-z/rbdavid/Apps/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/n-z/rbdavid/Apps/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate ena_db_build

################################################################################

SCHEDULER_FILE=/private_stores/gerlt2/users/rbdavid/testing/scheduler_file.json

# print some info about nodes
echo $SLURM_NNODES
echo $SLURM_NODEID
echo $SLURM_JOB_NODELIST
echo $SLURM_TASKS_PER_NODE
echo $(scontrol show hostnames $SLURM_JOB_NODELIST)

# print some info about resources
echo $SLURM_CPUS_ON_NODE
echo $SLURM_JOB_CPUS_PER_NODE
echo $SLURM_MEM_PER_CPU

# ID the first node in the nodelist
PrimaryNode=$(echo $(scontrol show hostnames $SLURM_JOB_NODELIST) | awk '{print $1;}')
echo $PrimaryNode

# empty string to gather process IDs associated with the dask scheduler and worker(s) calls
dask_pids=""

echo "Spinning up the Scheduler"
dask scheduler --no-dashboard --no-jupyter --no-show --scheduler-file ${SCHEDULER_FILE} --interface 'eth0' > dask_scheduler.out 2>&1 &
dask_pids="$dask_pids $!"

# seems to be needed to ensure the scheduler has actually spun up before workers begin to be spun up. 
# like if writing the SCHEDULER_FILE to its path is slow, then the `dask worker` call may not start up the communications between workers and the scheduler
sleep 20

# spin up nWorkers workers, each with 1 core
echo "Spinning up the Workers"
srun -n 62 \
	dask worker --no-dashboard --no-nanny --reconnect --nthreads 1 --nworkers 1 --interface 'eth0' --scheduler-file ${SCHEDULER_FILE} > dask_worker.out 2>&1 &
dask_pids="$dask_pids $!"

# just do it again to make sure everyone is communicating with each other
sleep 20

echo "Starting Client Script"

# test the workflow on a small dataset
#python3 dask_tskmgr.py --db-config /home/n-z/rbdavid/test_efi.config --db-name efi_202412 --ena-paths /home/n-z/rbdavid/Projects/ENA_building/wgs/public /home/n-z/rbdavid/Projects/ENA_building/sequence /home/n-z/rbdavid/Projects/ENA_building/wgs/suppressed --output-dir /home/n-z/rbdavid/Projects/ENA_building/TEST/ --local-scratch /scratch/ --scheduler-file ${SCHEDULER_FILE} --n-workers 62 --tskmgr-log-file /home/n-z/rbdavid/Projects/ENA_building/testing_tskmgr.log

# run the workflow on all ENA files in /private_stores/mirror/ena/20241205
python3 dask_tskmgr.py --db-config /home/n-z/rbdavid/test_efi.config --db-name efi_202412 --ena-paths /private_stores/mirror/ena/20241205/sequence /private_stores/mirror/ena/20241205/wgs/public  /private_stores/mirror/ena/20241205/wgs/suppressed  --output-dir /private_stores/gerlt2/users/rbdavid/ENA_build/2025_02_05/ --local-scratch /scratch/ --scheduler-file ${SCHEDULER_FILE} --n-workers 62 --tskmgr-log-file /private_stores/gerlt2/users/rbdavid/ENA_build/2025_02_05/tskmgr.log

for pid in $dask_pids
do
        kill $pid
done

date

