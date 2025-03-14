#!/bin/bash
#SBATCH --partition=efi
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --job-name="ena_db_build"
#SBATCH --nodelist=compute-9-14

date

source env_setup.sh

conda activate ena_db_build

################################################################################

# empty string to gather process IDs associated with the dask scheduler and worker(s) calls
dask_pids=""

echo "Spinning up the Scheduler"
SCHEDULER_FILE=$WORKING_DIR/scheduler_file.json
dask scheduler --no-dashboard --no-jupyter --no-show --scheduler-file ${SCHEDULER_FILE} --interface 'eth0' > dask_scheduler.out 2>&1 &
dask_pids="$dask_pids $!"

# spin up N_WORKERS workers, each with 1 core
echo "Spinning up the Workers"
srun -n $N_WORKERS \
	dask worker --no-dashboard --no-nanny --reconnect --nthreads 1 --nworkers 1 --interface 'eth0' --scheduler-file ${SCHEDULER_FILE} > dask_worker.out 2>&1 &
dask_pids="$dask_pids $!"

# run the workflow on all ENA files in $EBROOTENA
echo "Starting Client Script"
python3 ena_dask_tskmgr --db-config $DB_CONFIG --db-name $DB_NAME --ena-paths $EBROOTENA/sequence $EBROOTENA/wgs/public  $EBROOTENA/wgs/suppressed  --output-dir $OUTPUT_DIR --local-scratch $SCRATCH_DIR --scheduler-file ${SCHEDULER_FILE} --n-workers $N_WORKERS --tskmgr-log-file $WORKING_DIR/tskmgr.log

for pid in $dask_pids
do
        kill $pid
done

date

