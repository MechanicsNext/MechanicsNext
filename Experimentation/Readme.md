<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MechanicsNext_Assets/mechanics_next_wordmark.png" width=64px>

# Experimentation

These resources relate to experimental characterization of microstructure and deformation in metals.
At this time, following experimental tecniques are covered:

1. High-energy X-ray diffraction microscopy (HEDM):
2. Digital image correlation (DIC):

Following tools are available:

1. [DIC_image_processing](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/DIC_image_processing): Visualization and post-processing of strain data measured using digital image correlation (DIC) and processed with [NCORR](http://www.ncorr.com/) software.
2. [HEDM-file_utils](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_file_utils): Linux shell functions to manipulate data in GE2, GE3, or TIFF file format that is typically acquired during far-field HEDM experiments. These shell functions can be executed like any other Linux command.
3. [HEDM_forward_modeling](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_forward_modeling): Tools to simulate far-field or near-field HEDM diffraction patterns.
4. [HEDM_logging](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_logging): A Microsoft Word template for HEDM *lab notebook*.
5. [HEDM_macros](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_macros): A collection of SPEC macros to automate in-situ tension/compression HEDM experiment.
6. [HEDM_preprocessing](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_preprocessing): A utility to segment Bragg spots in HEDM data. This code accompanies the publication: Paranjape, H. M. et al., *An Algorithm to Segment Overlapping Bragg Peaks in Far-field High-energy Diffraction Microscopy Area Diffraction Patterns from Deformed Single and Multi-phase Materials*. (under peer review) Journal of Applied Crystallography.
7. [HEDM_specimen_geometry](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_specimen_geometry): Schematics of specimens used to perform in-situ tension/compression experiments at HEDM facilities. These specimens are best machined using electron discharge machining (EDM).
8. [HEDM_visualization](https://github.com/MechanicsNext/MechanicsNext/tree/master/Experimentation/HEDM_visualization): Matlab scripts to post-process and visualize output from [MIDAS](https://github.com/marinerhemant/MIDAS) and [heXRD](https://github.com/praxes/hexrd) HEDM data analysis software.

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
