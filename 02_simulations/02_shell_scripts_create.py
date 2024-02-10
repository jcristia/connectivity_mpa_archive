# create shell scripts for each opendrift and biology run

import os

folder = 'jobs'
base_od = 'od'
base_bi = 'bi'

od = folder + '/' + base_od
bi = folder + '/' + base_bi

# There are 9 groups of MPAs for each date
# I need an opendrift script for each of them
# I only need 1 biology script per date
dates = [
    '1101','1105','1108',
    '1401','1405','1408',
    '1701','1705','1708'
]
groups = [str(i) for i in list(range(1,10))]

for date in dates:
    for group in groups:
        name = od + date + '_' + group + '.sh'
        with open(name, 'w', newline='\n') as rsh: # use unix line endings
            rsh.write(f'''\
#!/bin/bash

#SBATCH --time=1-12:00        # time (DD-HH:MM)
#SBATCH --mem-per-cpu=128G  # memory; default unit is megabytes

#SBATCH --mail-user=jcristia10@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --job-name={base_od + date + '_' + group}
#SBATCH --output=./outputlogs/%x_%j.out

module load singularity
singularity exec --home /home/jcristia/scratch/mpaconn opendrift_mpaconn.sif  python scripts/sim{date}/{base_od + group}.py
''')


for date in dates:
    for group in groups:
        name = bi + date + '_' + group + '.sh'
        with open(name, 'w', newline='\n') as rsh:
            rsh.write(f'''\
#!/bin/bash

#SBATCH --time=1-00:00        # time (DD-HH:MM)
#SBATCH --mem-per-cpu=4000M   # memory; default unit is megabytes
#SBATCH --array=0-7

#SBATCH --mail-user=jcristia10@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --job-name={base_bi + date + '_' + group}
#SBATCH --output=./outputlogs/%x_%j.out

module load singularity
singularity exec --home /home/jcristia/scratch/mpaconn biology_mpaconn.sif  python scripts/sim{date}/{base_bi + group}.py $SLURM_ARRAY_TASK_ID
''')