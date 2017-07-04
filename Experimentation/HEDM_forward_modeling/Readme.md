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

### Annotated example

| X         | Y         | Z         | Material   | Quat 1 | Quat 2  | Quat 3  | Quat 4 | Stretch 11 | Stretch 22 | Stretch 33 | Stretch 12 | Stretch 13 | Stretch 23 |
|-----------|-----------|-----------|------------|--------|---------|---------|--------|--------|--------|--------|--------|--------|--------|
| -490.0000 | -490.0000 | -480.0000 | NiTi_cubic | 0.3737 | -0.0757 | -0.0639 | 0.9234 | 1.0097 | 1.0080 | 1.0101 | 0.0087 | 0.0024 | 0.0077 |
| -490.0000 | -490.0000 | -470.0000 | NiTi_cubic | 0.3731 | -0.0744 | -0.0636 | 0.9279 | 1.0084 | 1.0035 | 1.0101 | 0.0070 | 0.0035 | 0.0064 |
| -490.0000 | -490.0000 | -460.0000 | NiTi_cubic | 0.3766 | -0.0758 | -0.0620 | 0.9273 | 1.0095 | 1.0049 | 1.0089 | 0.0081 | 0.0016 | 0.0055 |
| -490.0000 | -490.0000 | -450.0000 | NiTi_cubic | 0.3757 | -0.0728 | -0.0628 | 0.9258 | 1.0078 | 1.0054 | 1.0120 | 0.0098 | 0.0049 | 0.0083 |
| -490.0000 | -490.0000 | -440.0000 | NiTi_cubic | 0.3731 | -0.0686 | -0.0665 | 0.9256 | 1.0070 | 1.0039 | 1.0090 | 0.0090 | 0.0049 | 0.0083 |
| -490.0000 | -490.0000 | -430.0000 | NiTi_cubic | 0.3761 | -0.0763 | -0.0682 | 0.9244 | 1.0076 | 1.0051 | 1.0098 | 0.0081 | 0.0038 | 0.0088 |
| -490.0000 | -490.0000 | -420.0000 | NiTi_cubic | 0.3776 | -0.0773 | -0.0661 | 0.9288 | 1.0083 | 1.0060 | 1.0084 | 0.0069 | 0.0026 | 0.0046 |
| -490.0000 | -490.0000 | -410.0000 | NiTi_cubic | 0.3826 | -0.0724 | -0.0701 | 0.9238 | 1.0063 | 1.0068 | 1.0094 | 0.0082 | 0.0017 | 0.0065 |
| -490.0000 | -490.0000 | -400.0000 | NiTi_cubic | 0.3802 | -0.0681 | -0.0641 | 0.9252 | 1.0085 | 1.0074 | 1.0118 | 0.0105 | 0.0055 | 0.0085 |
| -490.0000 | -490.0000 | -390.0000 | NiTi_cubic | 0.3769 | -0.0755 | -0.0679 | 0.9250 | 1.0085 | 1.0062 | 1.0108 | 0.0069 | 0.0025 | 0.0090 |
| -490.0000 | -490.0000 | -380.0000 | NiTi_cubic | 0.3789 | -0.0772 | -0.0687 | 0.9233 | 1.0069 | 1.0036 | 1.0113 | 0.0095 | 0.0016 | 0.0045 |
| -490.0000 | -490.0000 | -370.0000 | NiTi_cubic | 0.3780 | -0.0689 | -0.0616 | 0.9312 | 1.0095 | 1.0064 | 1.0088 | 0.0066 | 0.0032 | 0.0047 |
| -490.0000 | -490.0000 | -360.0000 | NiTi_cubic | 0.3812 | -0.0740 | -0.0647 | 0.9292 | 1.0071 | 1.0043 | 1.0072 | 0.0084 | 0.0019 | 0.0083 |
| -490.0000 | -490.0000 | -350.0000 | NiTi_cubic | 0.3751 | -0.0715 | -0.0667 | 0.9233 | 1.0055 | 1.0051 | 1.0082 | 0.0103 | 0.0044 | 0.0046 |
| -490.0000 | -490.0000 | -340.0000 | NiTi_cubic | 0.3815 | -0.0758 | -0.0678 | 0.9312 | 1.0065 | 1.0069 | 1.0091 | 0.0084 | 0.0056 | 0.0055 |


## Acknowledgement

This application uses Open Source components. You can find the source code of their open source projects along with license information below. We acknowledge and are grateful to these developers for their contributions to open source.
* [MTEX project](https://mtex-toolbox.github.io/)
* [MATLAB diffraction tools](https://github.com/junspark/matlab_tools) by Jun-Sang Park
* The [NCORR](http://ncorr.com/) digital image correlation package
* The [heXRD](https://github.com/praxes/hexrd) high-energy X-ray diffraction package

## License

    This file is part of MechanicsNext.

    MechanicsNext is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MechanicsNext is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MechanicsNext.  If not, see <http://www.gnu.org/licenses/>.
