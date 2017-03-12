<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/mechanics_next_wordmark.png" width=64px>
#Far-field High-energy Diffraction Microscopy Spot Segmentation Utility

A utility to read raw GE2 or GE3 files obtained during a [far-field high-energy diffraction microscopy](https://pdfs.semanticscholar.org/0c56/6a8040f5d60674063d41a3628b1da8d5270a.pdf) (ff-HEDM) experiment and segment the spots
in the data using connected component segmentation, local intensity maxima finding, and watershed segmentation. This utlity is particularly useful when
the data consists of multiple phases such that the Debye-Scherrer rings from the two phases have large intensity contrast but small
radial separation.

##Installation
###Prerequisites
1. Linux system.
2. Python. [Miniconda](https://conda.io/miniconda.html) distribution is recommended.
3. heXRD package. See [installation instruction in the Wiki](https://github.com/MechanicsNext/MechanicsNext/wiki/heXRD-on-Stampede).
4. NumPy, SciPy, and Scikit-learn packages.
5. The utility can consume large memory (~ 1 TB), depending on the number of spots in the data.
6. Access to multiple processors will significantly speed-up the processing. Hence,
running this program on a workstation or a supercomputer (e.g., Stampede) is recommended.

###Installation steps
1. Install Miniconda if necessary.
2. Clone the MechanicsNext Git repo. `git clone https://github.com/MechanicsNext/MechanicsNext`
3. Note the directory where MechanicsNext was cloned. We will call it `MECHANICSNEXT_DIR`
4. Create a new configuration file (we will call it `config.yml`) based on the examples in the `examples_spot_maxima_finder` folder.
5. Make sure necessary GE2 data files, dark files, detector configuration file, and the material file are present at appropriate locations. Generally
detector and material files are kept in the same folder as the configuration file.

##Running the utility
Run the utility using `python MECHANICSNEXT_DIR/Experimentation/HEDM_preprocessing/spot_maxima_finder.py config.yml`

###Running on Stampede supercomputer
After following the Installation steps above, the utility can be run on Stampede using the following job file. Note that
the `largemem` queue must be used due to large memory requirements.

```
#!/bin/bash

#SBATCH -J JOBNAME          # Job name
#SBATCH -o JOBNAME.%j.out   # stdout; %j expands to jobid
#SBATCH -e JOBNAME.%j.err   # stderr; skip to combine stdout and stderr
#SBATCH -p largemem         # queue
#SBATCH -N 1                # Number of nodes, not cores (16 cores/node)
#SBATCH -n 32               # Total number of MPI tasks (if omitted, n=N)
#SBATCH -t 01:00:00         # max time

python MECHANICSNEXT_DIR/Experimentation/HEDM_preprocessing/spot_maxima_finder.py config.yml
```

##Basic algorithm

<p align="center">
  <img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/Experimentation_HEDM_preprocessing_algorithm.png" width=400px>
</p>

