# Core libraries
import copy
import logging
import os
import sys
import time
import warnings
# YAML config processing
import yaml
# Math/array magic
import numpy as np
from scipy import ndimage
# Parallelization for speed
from joblib import Parallel, delayed
import multiprocessing
# heXRD libraries for diffraction functions
from hexrd.xrd import experiment as expt
from hexrd.coreutil import initialize_experiment
from hexrd import matrixutil as mutil
from hexrd.xrd import rotations as rot
from hexrd.xrd import transforms_CAPI as xfcapi
#--

def write_ge2(filename, arr, nbytes_header=8192, pixel_type=np.uint16):
    '''
    Write a 3D array to a GE2 file.
    '''
    fid = open(filename, 'wb')
    fid.seek(nbytes_header)
    fid.write(arr.astype(pixel_type))
    fid.close()
#--

def get_diffraction_angles_MP(fwd_model_input_i):
    '''
    Parallel processing worker for eta, two-th, omega calculation
    '''

    hkls = fwd_model_input_i['hkls']
    chi = fwd_model_input_i['chi']
    rmat_c = fwd_model_input_i['rmat_c']
    bmat = fwd_model_input_i['bmat']
    wavelength = fwd_model_input_i['wavelength']
    defgrad_i = fwd_model_input_i['defgrad_i']

    omega0, omega1 = xfcapi.oscillAnglesOfHKLs(hkls, 
                                               chi, 
                                               rmat_c, 
                                               bmat, 
                                               wavelength,
                                               vInv=defgrad_i)
#    print hkls
#    print omega0
#    print omega1
    return np.concatenate((omega0, omega1), axis=0)

class Microstructure:
    '''
    Microstructural information for forward modeling. This object holds
    spatial information for material/phase, crystal orientation, lattice strain
    and material properties and detector parameters for the virtual experiment.
    '''
    def __init__(self, config, logger, datafile):
        self.cfg = config              # heXRD config object
        self.logger = logger           # logger
        self.ms_datafile = datafile    # microstructural data file
        self.ms_grid = []              # (N, 3) array of X, Y, Z positions in microns
        self.ms_material_ids = []      # (N) Array of material present at X, Y, Z
        self.ms_quaternions = []       # (N, 4) Orientation at X, Y, Z in quaternions
        self.ms_lat_strains = []       # (N, 6) Lattice strain a X, Y, Z
        self.intensity_factors = []    # (N, 1) Spot intensity will be scaled by this factor
        self.synth_angles = []         # Two-theta, eta, omega from virtual diffraction
        self.calc_xyo = []             # X, Y, projected on the detector and omega

	# Initialize detector and reader from the experiment. Really only detector is needed.
        pd, reader, detector = initialize_experiment(config)
        # need instrument cfg later on down...
        instr_cfg = get_instrument_parameters(cfg)
        detector_params = np.hstack([
            instr_cfg['detector']['transform']['tilt_angles'],
            instr_cfg['detector']['transform']['t_vec_d'],
            instr_cfg['oscillation_stage']['chi'],
            instr_cfg['oscillation_stage']['t_vec_s'],
        ])

        if instr_cfg['detector']['distortion']['function_name'] == 'GE_41RT':
            distortion = (dFuncs.GE_41RT,
                          instr_cfg['detector']['distortion']['parameters'],
                          )
        else:
            distortion = None

	self.detector = detector
	self.reader = reader
        self.detector_params = detector_params
        self.distortion = distortion

    def read_csv(self):
        ''' 
        Load microstructural data from a csv file. Required columns are
        Position_X, position_Y, position_Z, material_name, 
        orientation_quat_1, orientation_quat_2, orientation_quat_3, orientation_quat_4,
        strain_11, strain_22, strain_33, 
        strain_12, strain_23, strain_31, 
        '''

        filename = self.ms_datafile
        logger = self.logger

        try:
            ms_data = np.loadtxt(filename, dtype=None, comments='#', delimiter=',',
                                 usecols=(0,1,2,4,5,6,7,8,9,10,11,12,13),
                                 ndmin=2)
            ms_mat_data = np.loadtxt(filename, dtype='str', comments='#', delimiter=',',
                                  usecols=(3,), ndmin=2)

            ms_data = np.array(ms_data)
            ms_mat_data = np.array(ms_mat_data)
        except:
            logger.error('Could not read microstructural data from %s', filename)
            ms_data = None
            ms_material_data = None

        self.ms_grid = ms_data[:, range(0, 3)]
        self.ms_material_ids = ms_mat_data[:, 0]
        self.ms_quaternions = ms_data[:, range(3, 7)]
        self.ms_lat_defgrads = ms_data[:, range(7, 13)]

        # Try reading intensity factors (15th column)
        try:
            int_factor_data = np.loadtxt(filename, dtype=None, comments='#', delimiter=',',
                              usecols=(14,),
                              ndmin=2)
        except:
            logger.error('Could not read intensity scaling data from %s', filename)
            int_factor_data = np.ones_like(ms_data)
            int_factor_data = int_factor_data[:, 0]

        self.intensity_factors = int_factor_data


    def get_diffraction_angles(self):
        cfg = self.cfg
        logger = self.logger
	detector = self.detector
        intensity_factors = self.intensity_factors

        # Number of cores for multiprocessing
        num_cores = cfg.multiprocessing

        # Initialize a new heXRD experiment
        ws = expt.Experiment()

        cwd = cfg.working_dir

        materials_fname = cfg.material.definitions
        material_name = cfg.material.active
        detector_fname = cfg.instrument.detector.parameters_old

        # load materials
        ws.loadMaterialList(os.path.join(cwd, materials_fname))
        mat_name_list = ws.matNames

        # Create an array of all the essential parameters to be
        # sent to parallel diffraction angle calculation routine
        fwd_model_input = []

        for xyz_i, mat_name_i, quat_i, defgrad_i in zip(self.ms_grid,
                                                       self.ms_material_ids,
                                                       self.ms_quaternions,
                                                       self.ms_lat_defgrads):

	    # Set chi tilt
            if detector.chiTilt is not None:
		chi = detector.chiTilt
	    else:
	        chi = 0.0

	    # Obtain all symmetric hkls for the material
            ws.activeMaterial = mat_name_i.strip() #material_name

	    hkls = ws.activeMaterial.planeData.getSymHKLs()
	    hkls = np.transpose(np.hstack(hkls))
	    # Rotational matrix from the orientation/quaternion
            rmat_c = rot.rotMatOfQuat(quat_i)
	    # bmat
            bmat = ws.activeMaterial.planeData.latVecOps['B']
            wavelength = ws.activeMaterial.planeData.wavelength
	    # Create a dictionary of inputs to be sent to the MP worker
            fwd_model_input.append(
                {
                    'chi': chi,
                    'hkls': hkls,
                    'rmat_c': rmat_c,
                    'bmat': bmat,
                    'wavelength': wavelength,
                    'defgrad_i': defgrad_i
                    }
                )

        # Now we have the data, run eta, twoth, omega calculations in parallel
        logger.info("Starting virtual diffraction calculations using %i processors", num_cores)
        synth_angles_MP_output = Parallel(n_jobs=num_cores, verbose=5)(delayed(get_diffraction_angles_MP)(fwd_model_input_i) for fwd_model_input_i in fwd_model_input)

        intensity_factors_spot_tmp = []
        for intensity_factor_ii, synth_angles_ii in zip(intensity_factors, synth_angles_MP_output):
            int_factor_tmp = np.ones_like(synth_angles_ii[:, 0]) * intensity_factor_ii
            int_factor_tmp_t = np.reshape(int_factor_tmp, (int_factor_tmp.shape[0], 1))
            intensity_factors_spot_tmp.append(int_factor_tmp_t)

        synth_angles = np.vstack(synth_angles_MP_output)
        intensity_factors_spot = np.vstack(intensity_factors_spot_tmp)

        self.synth_angles = synth_angles
        self.intensity_factors_spot = intensity_factors_spot
        return synth_angles

    def project_angs_to_detector(self, output_txt=None, output_file=None):
        '''
        Project two-theta, eta, omega angle to detector X, Y, omega
        '''
        cfg = self.cfg
        angs = self.synth_angles
	detector = self.detector
	# Obtain X, Y from two-theta, eta
        calc_xyo = detector.angToXYO(angs[:, 0], angs[:, 1], angs[:, 2])
        calc_xyo = np.transpose(calc_xyo)
        # Write X, Y, omega data to a text file
        template = "{0:12.2f}{1:12.2f}{2:12.2f}{3:12.2f}{4:12.2f}"
        template_file = "{0:12.2f}{1:12.2f}{2:12.2f}{3:12.2f}{4:12.2f}\n"

        if output_txt:
            if output_file is not None:
	        f = open(output_file, 'w+')
	    else:
	        f = open("synth_data.out", 'w+')

            for x, y, o, tt, e, o2 in np.hstack((calc_xyo, angs)):
                #print template.format(x, y, o)
	        if o < 0.:
		    o = o + 2.0*np.pi

	        f.write(template_file.format(x, y, o, e, tt))

	    f.close()

        self.calc_xyo = calc_xyo
        return calc_xyo

    def write_xyo_to_ge2(self, output_ge2=None, omega_start=None, omega_step=None, omega_stop=None, ge2_blur_sigma=3):
        '''
	Prepare a 3D array of size OMEGAS x XSIZE x YSIZE from X, Y, omega
	and write the data to a GE2 file.
        '''
	calc_xyo = self.calc_xyo
	detector = self.detector
        intensity_factors_spot = self.intensity_factors_spot
	# Get X, Y, omega dimensions to create a 3D array
	if output_ge2 is None:
	    output_ge2 = 'ff_synth_00000.ge2'

	if omega_start is None:
	    omega_start = 0.0

        if omega_step is None:
            omega_step = 0.1

        if omega_stop is None:
            omega_stop = 360.0 

	o_dim = int(np.round((omega_stop - omega_start)/omega_step))
	x_dim = detector.get_ncols()
        y_dim = detector.get_nrows()
	# Create an empty array of appropriate size
	synth_array = np.zeros((o_dim, x_dim, y_dim))
	# Fill in intensity details at appropriate X, Y, ome positions
	for x, y, o, int_scale_factor in np.hstack((calc_xyo, intensity_factors_spot)):
            if o < 0.:
                o = o + 2.0*np.pi

            o = (o * 180.0 / np.pi) 

	    if o > omega_start and o < omega_stop:
	    	x_op = np.round(x)
	    	y_op = np.round(y)
	    	o_op = np.round((o - omega_start)/360.0*(omega_stop - omega_start)/omega_step)
		# CHeck for sane indices
	    	if x_op < 0:
		    x_op = 0
	    	if x_op >= x_dim:
		    x_op = x_dim - 1
            	if y_op < 0:
                    y_op = 0
            	if y_op >= y_dim:
                    y_op = y_dim - 1
	    	if o_op < 0:
		    o_op = 0
	    	if o_op >= ((omega_stop - omega_start)/omega_step):
		    o_op = (omega_stop - omega_start)/omega_step - 1
		# For now we are not calculating structure factor and realistic intensities for spots.
		# Can we use heXRD crystallography functions? May be. Need to find a way of getting atom positions.
		# Read in CIFs? haha fu*k no..
		synth_array[o_op][x_op][y_op] = 12000 * int_scale_factor
        # Scale intensities so that they are really between 0, 12000
        max_intensity = np.amax(synth_array)
        final_int_scale_factor = 12000. / max_intensity
        synth_array = synth_array * final_int_scale_factor
        # Add a Gauss blur
        if ge2_blur_sigma > 0:
            synth_array_blurred = ndimage.filters.gaussian_filter(synth_array, sigma=ge2_blur_sigma)
        else:
            synth_array_blurred = synth_array

	# Write the array to a GE2
	write_ge2(output_ge2, synth_array_blurred)
	# Also write a max-over frame for now.
	synth_array_max = np.amax(synth_array_blurred, axis=0)
	output_ge2_max = output_ge2.rsplit('.', 1)[0]
	write_ge2(output_ge2_max + '-max.ge2', synth_array_max)

	self.synth_array = synth_array
	return synth_array
#--
