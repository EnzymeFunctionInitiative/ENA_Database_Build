#!/bin/bash
#SBATCH --partition=efi
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --job-name="pickling"
#SBATCH --nodelist=compute-9-14

date

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

python3 pickle_ena.py /private_stores/gerlt2/databases/ena/20241217/output/ena.tab /private_stores/gerlt2/users/rbdavid/ENA_build/efi_202412_2025_02_20/perl_ena.pkl

python3 pickle_ena.py /private_stores/gerlt2/users/rbdavid/ENA_build/efi_202412_2025_02_20/ena.tab /private_stores/gerlt2/users/rbdavid/ENA_build/efi_202412_2025_02_20/ena.pkl

python3 compare_pickles.py /private_stores/gerlt2/users/rbdavid/ENA_build/efi_202412_2025_02_20/perl_ena.pkl /private_stores/gerlt2/users/rbdavid/ENA_build/efi_202412_2025_02_20/ena.pkl /home/n-z/rbdavid/Projects/ENA_building/comparison_code

date

