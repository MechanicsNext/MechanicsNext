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

## Acknowledgement

This source code is courtesy of Darren Pagan (Cornell High Energy Synchrotron Source) and heXRD developers.

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
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
