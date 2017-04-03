<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/mechanics_next_wordmark.png" width=64px>

# Brute-force Near-field HEDM Detector Calibration

Does a parameter study on the four detector parameters (distance, xcenter, ytilt, ztilt)
in a nf-HEDM setup. The result is a text file with the number of grid points
with the confidence above a threshold. For a good calibration, this number should
be large.

Best if run on a supercomputer. Sample job file for XSEDE Stampede is below.

```
!/bin/bash

#SBATCH -J NFCAL          # Job name
#SBATCH -o NFCAL.%j.out   # stdout; %j expands to jobid
#SBATCH -e NFCAL.%j.err   # stderr; skip to combine stdout and stderr
#SBATCH -p largemem       # queue
#SBATCH -N 1              # Number of nodes, not cores (16 cores/node)
#SBATCH -n 32             # Total number of MPI tasks (if omitted, n=N)
#SBATCH -t 20:00:00       # max time

python nf_detector_calibration_bruteforce.py
```
