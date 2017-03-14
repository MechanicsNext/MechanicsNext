<img src="https://github.com/MechanicsNext/MechanicsNext/blob/master/MeshnicsNext_Assets/mechanics_next_wordmark.png" width=64px>
# Far-field High-energy Diffraction Microscopy Spot Segmentation Utility: Examples

These examples can be used as templates to run spot segmentation on other datasets.

## local_maxima_cubic_monoclinic_mixed

An example showing spot segmentation in a two-phase material. The spots from the cubic phase
are significantly brighter than the spots from the monoclinic phase. Thus an upper
intensity cutoff can be used to eliminate cubic spots while preserving the monoclinic spots.

Following files are included:

* config_local_maxima_cubic_monoclinic_mixed.yml: Configuration file.
* local_maxima_cubic_monoclinic_mixed.data: Sample text output from the program.
* material_local_maxima_cubic_monoclinic_mixed: heXRD material specification file. Can be read using the [heXRD](https://github.com/praxes/hexrd) GUI.
* detector_local_maxima_cubic_monoclinic_mixed: heXRD detector specification in the old (binary) format. Can be read using the [heXRD](https://github.com/praxes/hexrd) GUI.
* local_maxima_cubic_monoclinic_mixed.job: A sample job file for running the spot segmentation utility on the Stampede supercomputer. Can be adapted for other scheduling environments.

