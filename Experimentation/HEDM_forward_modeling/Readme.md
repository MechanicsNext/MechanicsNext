<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/mechanics_next_wordmark.png" width=64px>

# Far-field High-energy Diffraction Microscopy Spot Simulation (Forward Modeling) Utility

A utility to simulate ff-HEDM area diffraction patterns from microstructures. The microstructure is specified in terms of
the position, crystal structure, and elastic strain of each material point in a discrete mesh generated over a virtual 
specimen. Diffraction patterns are simulated using forward modeling routines in [heXRD](https://github.com/praxes/hexrd).
Optionally structure factors calculated using [Pymatgen](http://pymatgen.org/) are used. The output is written to text
and GE2 files.

## Installation
### Prerequisites
1. Linux system.
2. Python. [Miniconda](https://conda.io/miniconda.html) distribution is recommended.
3. heXRD package. See [installation instruction in the Wiki](https://github.com/MechanicsNext/MechanicsNext/wiki/heXRD-on-Stampede).
4. NumPy, SciPy, and Scikit-learn packages.
5. [Pymatgen](http://pymatgen.org/#conda-install-recommended) library.
6. Access to multiple processors will significantly speed-up the processing. Hence,
running this program on a workstation or a supercomputer (e.g., Stampede) is recommended.

### Installation steps
1. Install Miniconda if necessary.
2. Clone the MechanicsNext Git repo. `git clone https://github.com/MechanicsNext/MechanicsNext`
3. Note the directory where MechanicsNext was cloned. We will call it `MECHANICSNEXT_DIR`
4. Create a new configuration file (we will call it `config.yml`) based on the examples in the `examples_synth_fwd_modeling_patterns` folder.
5. Make sure necessary detector configuration file, CIF file (for optional structure factor calculation), and the material file are present at appropriate locations. Generally
detector and material files are kept in the same folder as the configuration file.

## Running the utility
Run the utility using `python MECHANICSNEXT_DIR/Experimentation/HEDM_forward_modeling/synth_fwd_modeling_patterns.py config.yml`

### Running on Stampede supercomputer
After following the Installation steps above, the utility can be run on Stampede using the following job file.

```
#!/bin/bash

#SBATCH -J JOBNAME          # Job name
#SBATCH -o JOBNAME.%j.out   # stdout; %j expands to jobid
#SBATCH -e JOBNAME.%j.err   # stderr; skip to combine stdout and stderr
#SBATCH -p largemem         # queue
#SBATCH -N 1                # Number of nodes, not cores (16 cores/node)
#SBATCH -n 32               # Total number of MPI tasks (if omitted, n=N)
#SBATCH -t 01:00:00         # max time

python MECHANICSNEXT_DIR/Experimentation/HEDM_forward_modeling/synth_fwd_modeling_patterns.py config.yml
```

## Basic algorithm

<p align="center">
  <img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/Experimentation_HEDM_forward_modeling_algorithm.png" width=400px>
</p>


## Microstructure file format
Imagine a material volume whose far-field area diffraction patterns are to be simulated. Generate a grid of points covering the whole volume. Assign position, material, orientation, and elastic strain information to each point. The mesh along with the information associated with each point is a *microstructure*. This microstructure is described in a text file containing 14 columns.

* **Columns 1 to 3**: Location of a point in the grid. In microns. It is a good practice to create a grid that goes from -X to X, -Y to Y, -Z to Z.
* **Column 4**: Material name. Thre should be a corresponding material in the heXRD material file provided. The name should not contain spaces. Crystal structure and lattice parameters from this material will be in determining spot positions.
* **Columns 5 to 8**: Crystal orientation of the grid point in quaternions.
* **Columns 9 to 14**: Elastic right stretch tensor at the point. i.e., *U* in *F = RU*, where *F* is the deformation gradient, *U* is the stretch, and *R* is the rotation.

