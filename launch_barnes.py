import os
import numpy as np

sim_values = np.array([200, 250, 300, 350, 400, 450])

main_str = 'gk_Barnes.py'
job_str = 'ABC_Barnes'

Ncores = 1 # Number of requested cores
Nhours = 4 # Requsted walltime (hours)

for i in range(0,len(sim_values)):
    nSim = sim_values[i]
    jobId = job_str + '_' + str(nSim)
    slurmFileName = 'run_' + jobId + '.slurm'
    
    with open(slurmFileName, 'w') as fout:
        fout.write('#!/bin/bash -l\n')
        fout.write('\n')
        fout.write('#SBATCH --job-name ' + jobId + '\n')
        fout.write('#SBATCH --nodes 1\n')
        fout.write('#SBATCH --cpus-per-task ' + str(Ncores) + '\n')
        fout.write('#SBATCH --time ' + str(Nhours) + ':00:00\n')
        fout.write('#SBATCH --output ' + 'out_' + jobId + '.txt\n')
        fout.write('#SBATCH --mem=4000\n')
        fout.write('#SBATCH --mail-type=FAIL\n')
        fout.write('#SBATCH --mail-user=louis.filstroff@aalto.fi\n')
        fout.write('\n')
        fout.write('module load anaconda\n')
        fout.write('activate TOAL\n')
        fout.write('python ' + main_str + ' ' + str(nSim))
    
    os.popen('sbatch ' + slurmFileName)
    os.remove(slurmFileName)