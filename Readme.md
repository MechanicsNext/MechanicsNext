<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/mechanics_next_wordmark.png" width=64px>

# MechanicsNext

Mechanics of Materials Using Numerical and Experimental Techniques.

## Introduction

MechanicsNext is a framework to assist the discovery and elucidation of multi-scale mechanics phenomena in inelastically deforming materials. This framework is inspired by the belief that combined modeling and experimental efforts can furnish higher impact results compared to either modeling-only or simulation-only efforts. The framework aims to provide tools for experiment design, analysis of experimental data, design of simulations and bridges to link experiments with simulations. The tools are partitioned in four classes.

1. Experimentation: These tools include utilities to automate high-energy X-ray diffraction experiments at Argonne National Laboratory and Cornell High Energy Synchrotron Source, scripts for the post processing of high-energy X-ray diffraction (HEDM) results as well as tools for processing data from other experimental methods.
2. Modeling: These tools include utilities for tessellation, batch-creation of Abaqus simulations and processing Abaqus finite element simulation results. Abaqus user material subroutines (UMAT) for phase transformation and plasticity (microstructural as well as macroscopic/phenomenological) will be made publically available in the future.
3. Infrastructure: Tools in this class include utilities for pre-processing (HEDM) data, simulation of far-field HEDM diffraction patterns and tutorials for setting-up and using certain diffraction-related softwares on workstations and clusters.
4. Teaching and Outreach: Resources in this section consist of software demonstrations of basic materials phenomena or methods (e.g. phase change, plasticity) and software-apps that are designed for external outreach.

Requests for collaborations can be addressed to [Harshad Paranjape](mailto:contact@harshadparanjape.com). General inquiries can be made by opening an issue on this repo.

## Cotributors
1. Harshad Paranjape (Colorado School of Mines)
2. Branden Kappes (Colorado School of Mines)

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
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
