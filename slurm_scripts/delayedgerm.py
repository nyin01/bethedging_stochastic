import itertools
import datetime

# parameters
# Args --> ##Args --> $filename $SLURM_ARRAY_TASK_ID $env $germ $hard $pA $pGermBH $pGermWT $bankLength $wB

all_envs = [1]
all_germs = [1]
all_hards = [1, 0]
all_pA = [0, 0.25, 0.5, 0.75, 1]
all_pGermBH = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
all_pGermWT = [0.7, 0.8, 0.9, 1]
all_bankLength = [5]
all_wB = [0.75, 0.8, 0.85, 0.9, 0.95, 0.99]

size = len(all_envs) * len(all_germs) * len(all_hards) * len(all_pA) * \
    len(all_pGermBH) * len(all_pGermWT) * len(all_bankLength) * len(all_wB)
job_name = 'delayedgerm_sim'

job_script = """#!/bin/bash

#SBATCH -J {job_name}           # Job name.

#SBATCH -N 1                               # Number of cores.
#SBATCH -t 25:00:00                         # Runtime in D-HH:MM (or use minutes).
#SBATCH --mem 4g                           # Memory in GB.
#SBATCH --array=1-{}
#SBATCH --mail-type=ALL                    # Type of email notification: BEGIN,END,FAIL,ALL.
#SBATCH --mail-user=zheng_yin@brown.edu   # Email where notifications will be sent.
#SBATCH -o {job_name}-%A.out
#SBATCH -e {job_name}-%A.err

declare -a all_envs=({})
declare -a all_germs=({})
declare -a all_hards=({})
declare -a all_pA=({})
declare -a all_pGermBH=({})
declare -a all_pGermWT=({})
declare -a all_bankLength=({})
declare -a all_wB=({})

# Get the parameters for this job
current_index=$((SLURM_ARRAY_TASK_ID-1))
envs=${{all_envs[$current_index]}}
germs=${{all_germs[$current_index]}}
hards=${{all_hards[$current_index]}}
pA=${{all_pA[$current_index]}}
pGermBH=${{all_pGermBH[$current_index]}}
pGermWT=${{all_pGermWT[$current_index]}}
bankLength=${{all_bankLength[$current_index]}}
wB=${{all_wB[$current_index]}}

now="$(date +'%m.%d.%Y')"
filename="delayedgerm_sim.$now"

module load julia/1.9.3s
julia delayedgerm_sim.jl $envs $germs $hards $pA $pGermBH $pGermWT $bankLength $wB

"""

# Generate all combinations of the parameters
combinations = list(itertools.product(all_envs, all_germs, all_hards,
                    all_pA, all_pGermBH, all_pGermWT, all_bankLength, all_wB))

# Generate the parameter arrays for each combination
all_envs_array = [str(x[0]) for x in combinations]
all_germs_array = [str(x[1]) for x in combinations]
all_hards_array = [str(x[2]) for x in combinations]
all_pA_array = [str(x[3]) for x in combinations]
all_pGermBH_array = [str(x[4]) for x in combinations]
all_pGermWT_array = [str(x[5]) for x in combinations]
all_bankLength_array = [str(x[6]) for x in combinations]
all_wB_array = [str(x[7]) for x in combinations]


# Generate the job script content
job_array_size = len(combinations)
job_script_content = job_script.format(
    job_array_size,
    ' '.join(all_envs_array),
    ' '.join(all_germs_array),
    ' '.join(all_hards_array),
    ' '.join(all_pA_array),
    ' '.join(all_pGermBH_array),
    ' '.join(all_pGermWT_array),
    ' '.join(all_bankLength_array),
    ' '.join(all_wB_array),
    script=job_script,
    job_name=job_name
)

# Write the job script to a file
with open(job_name + '.sh', 'w') as file:
    file.write(job_script_content)
