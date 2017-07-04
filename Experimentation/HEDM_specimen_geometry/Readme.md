<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MechanicsNext_Assets/mechanics_next_wordmark.png" width=64px>

# Specimen Geometry Specifications for Various Diffraction-based Experiments

1. `ALS_microlaue_custom_load_frame_sample.DXF`: Planar specimen geometry for in-situ loading experiments at the Laueu micro-diffraction beamline (12.3.2) at the Advanced Light Source (ALS). The gage length can be changed by a few millimeters.
2. `CHESS_RAMS2_design_1_long_gage.PDF`: Specimen geometry with longer gage for in-situ loading experiemnts using the RAMS2 load frame at F2 station of the Cornell High-energy Synchrotron Source (CHESS). High-energy diffraction microscopy (HEDM) experiments can be performed at F2. This geometry is not recommended. The maximum beam height (1 mm) is substantially smaller than the gage length (8 mm) and so deformation can localize outside the beam.
3. `CHESS_RAMS2_design_1_reduced_gage.pdf`: Preferred geometry to perform in-situ HEDM experiments at F2 station of CHESS. The gage is tapered and the deformation is likely to localize in the middle of the gage.
4. `CHESS_RAMS2_design_1_reduced_gage_alternate.PDF`: Alternate reduced-gage design for F2 @ CHESS.
5. `CHESS_RAMS2_design_3_long_round_gage.PDF`: Longer gage design for F2 @ CHESS. Not recommended.


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
