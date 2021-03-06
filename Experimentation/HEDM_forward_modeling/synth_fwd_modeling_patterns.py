###############################################################################
# MechanicsNext: HEDM Diffraction Pattern Simulation Program
# (https://github.com/MechanicsNext/MechanicsNext/blob/master/Experimentation/HEDM_preprocessing/)
#
# Simulate far-field high-energy diffraction microscopy patterns from
# the microstructural data provided.
#
# USAGE:
#  python synth_fwd_modeling_patterns.py config_file.yml
#
# See examples_synth_fwd_modeling_patterns/ folder for example configurations.
#
# Written by Harshad Paranjape (hparanja@mines.edu) and contributors.
# All rights reserved.
# Released under GNU Lesser General Public License v2.1 or above
# https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html
#
###############################################################################
import copy
import logging
import os
import sys
import time
import warnings

import yaml

from hexrd import config
from hexrd.coreutil import make_old_detector_parfile
# Forward modeling library
from forward_modeling.fwd_modeling_from_micro import *
from forward_modeling.microstructure_generation import *

if __name__ == '__main__':
    # Read args
    if len(sys.argv) < 2:
        print 'USAGE: python synth_fwd_modeling_patterns.py config.yml'
        sys.exit(1)
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print 'USAGE: python synth_fwd_modeling_patterns.py config.yml'
        sys.exit(1)
    else:
        cfg_file = sys.argv[1]

    # Setup logger
    log_level = logging.DEBUG
    logger = logging.getLogger('hexrd')
    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    cf = logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
    ch.setFormatter(cf)
    logger.addHandler(ch)
    # load the configuration settings
    cfgs = config.open(cfg_file)

    # cfg is a list. We will loop over each cfg set.
    for cfg in cfgs:
    	logger.info('=== begin forward-modeling ===')

        # Are we doing FF or NF pattern simulation?
        try:
            hedm_mode = cfg.get('forward_modeling')['fwdmodel']['hedm_mode'].strip().lower()
        except:
            hedm_mode = 'ff'
            logger.error('HEDM mode not specified. Defaulting to far-field.')

        # If HEDM mode is nf, convert new detector parameter file to old format
        if hedm_mode == 'nf':
            try:
                with open(cfg.get('instrument')['parameters'].strip(), 'r') as stream:
                    new_det_data = yaml.load(stream)
		#

                old_det_file_name = cfg.get('instrument')['detector']['parameters_old'].strip()
                logger.info('Writing old detector configuration to %s', old_det_file_name)

	    except:
                logger.error('For near-field mode, new detector file and the name for old file is required.')

            new_det_data_copy = {
            	'tiltAngles': np.array(new_det_data['detector']['transform']['tilt_angles']),
            	'tVec_d': np.array(new_det_data['detector']['transform']['t_vec_d']).reshape(3,1),
            	'tVec_s': np.array(new_det_data['oscillation_stage']['t_vec_s']).reshape(3,1),
            	'dParams': 0.}
            make_old_detector_parfile(new_det_data_copy, \
	  		          det_origin=(float(new_det_data['detector']['pixels']['rows']) * float(new_det_data['detector']['pixels']['size'][0]) / 2., \
				              float(new_det_data['detector']['pixels']['columns']) * float(new_det_data['detector']['pixels']['size'][1]) / 2.), \
			          filename=old_det_file_name)

        # Are we doing forward modeling or generating a microstructure for fwd modeling?
        try:
            fwd_model_mode = cfg.get('forward_modeling')['modeling_mode'].strip().lower()
        except:
            logger.error('Invalid forward modeling mode. Choices are datagen and fwdmodel')

	# If 'datagen' then generate synthetic microstructure data.
	# Certain functions are predefined in the microstructure_generation library.
	# Implement other functions there.
        if fwd_model_mode == "datagen":
            try:
	        fwd_model_nipt = cfg.get('forward_modeling')['datagen']['number_of_points']
            except:
                logger.error('Invalid value for number of data points to generate for forward modeling input.')

            try:
	        fwd_model_op_file = cfg.get('forward_modeling')['datagen']['output_file_name'].strip()
            except:
                logger.info('Invalid name for text output file. Defaulting to ms-synth.csv')
                fwd_model_op_file = 'ms-synth.csv'

            logger.info('=== generating synthetic microstructural data ===')
	    logger.info('=== writing output to %s ===', fwd_model_op_file)

	    # Generate a monoclinic single crystal with mosaicity and strain spread.
	    generate_mono_grain_mosaicity(nipt=fwd_model_nipt, output_file=fwd_model_op_file, \
					  material_name='NiTi_mono', mosaicity=0.0001, defgrad_spread=0.0001)
	    # Generate a cubic single crystal with mosaicity and strain spread.
	    #generate_cubic_grain_mosaicity(nipt=fwd_model_nipt, output_file=fwd_model_op_file)
	    # Generate an ideal cubic single crystal - without mosaicity or strain.
            #generate_cubic_grain_ideal(nipt=fwd_model_nipt, output_file=fwd_model_op_file)
	    # Generate an ideal polycrystal. Each point is a grain. No mosaicity or strain.
            #generate_cubic_grains_random_ideal(nipt=fwd_model_nipt, output_file=fwd_model_op_file)

	# If 'fwdmodel' then simulate diffraction patterns
        elif fwd_model_mode == "fwdmodel":
            # Get the fwdmodel mode: centroids or strainonly
            try:
                fwdmodel_mode = cfg.get('forward_modeling')['fwdmodel']['fwdmodel_mode'].strip().lower()
                logger.info("Setting FWDMODEL mode to '%s'.", fwdmodel_mode)
            except:
                logger.info("FWDMODEL mode is not specified. Defaulting to 'centroids'.")
                fwdmodel_mode = 'centroids'

	    # Get microstructural input file name. See Github documentation for the format of this file.
            try:
                fwd_model_ip_filename = \
                    cfg.get('forward_modeling')['fwdmodel']['input_file_name'].strip()
                logger.info('Microstructural data for FWDMODEL will be read from %s.', fwd_model_ip_filename)
            except:
                fwd_model_ip_filename = 'ms-data.csv'
                logger.info('Could not get the name for microstructural input file. Defaulting to %s.', fwd_model_ip_filename)
	    # Get the file name for output data -- two-theta, eta, omega for spots.
            try:
                fwd_model_op_filename = cfg.get('forward_modeling')['fwdmodel']['output_file_name'].strip()
            except:
                logger.info('Invalid output text file name. Defaulting to synth-data.out.')
                fwd_model_op_filename = 'synth-data.out'
	    # Get CIF file name for structure factor calculation. If CIF is not given, a structure
            # factor value of 1 is used for all rings.
            try:
                cif_filename = cfg.get('material')['cif'].strip()
                logger.info('Reading CIF data from %s..', cif_filename)
            except:
                logger.info('CIF file for structure not provided. Using a uniform structure factor for rings.')
                cif_filename = None

	    # Create a Microstructure object. This stores all data related to input and simulated output.
	    # TODO: Decide if we want to read everything from cfg into Microstructure or from cfg to here.
            ms = Microstructure(cfg, logger, fwd_model_ip_filename, cif_filename)
	    # Read microstructural input froma  CSV file.
            try:
                ms_filetype = cfg.get('forward_modeling')['fwdmodel']['input_type'].strip()
            except:
                ms_filetype = 'ms'

            logger.info('Reading microstructure.')
            ms.read_csv(filetype=ms_filetype)
            logger.info('Finished reading microstructure.')

            # Forward modeling mode

            if fwdmodel_mode is 'strainonly':
	        # Obtain diffraction angles using routines implemented in heXRD. This takes into 
                # account everything -- strain, position, orientation.
                logger.info('Starting simulation of diffraction patterns. Using strain, position, orientation data.')
                ms.get_diffraction_angles()
                logger.info('Finished simulation of diffraction patterns.')
            else:
                logger.info('Starting simulation of diffraction patterns.')
                # This mode does not consider the effects of having a finite sized sample.
                ms.simulate_pattern_to_detector()
                logger.info('Finished simulation of diffraction patterns.')

            # Write output (spot angles, position etc.) to aa text file?
            try:
                output_txt_flag = cfg.get('forward_modeling')['fwdmodel']['output_txt']
            except:
                output_txt_flag = False

            if output_txt_flag is not False:
                logger.info('Writing text output to %s', fwd_model_op_filename)
	        # Project the two-theta, eta, omega angles to X, Y using heXRD detector routines.
                ms.project_angs_to_detector(output_txt=output_txt_flag, output_file=fwd_model_op_filename)
            # Write output to a GE2 file?
	    try:
	        output_ge2_flag = cfg.get('forward_modeling')['fwdmodel']['output_ge']
	    except:
		output_ge2_flag = False
            # Write output to a TIF files?
            try:
                output_tif_flag = cfg.get('forward_modeling')['fwdmodel']['output_tif']
                output_tif_stem = cfg.get('forward_modeling')['fwdmodel']['output_tif']
            except:
                output_tif_flag = False

	    if output_ge2_flag is not False:
                # Write output to GE2. First get the file name.
                try:
                    output_ge2 = cfg.get('forward_modeling')['fwdmodel']['output_ge_name']
                    logger.info('GE2 output will be written to %s.', output_ge2)
                except:
                    logger.info('Invalid name for output GE2 file. Defaulting to ff_00000.ge2.')
                    output_ge2 = 'ff_00000.ge2'
                # Get omega angular bounds (degrees)
                try:
                    omega_start = cfg.get('forward_modeling')['fwdmodel']['output_omega']['start']
                    omega_step = cfg.get('forward_modeling')['fwdmodel']['output_omega']['step']
                    omega_stop = cfg.get('forward_modeling')['fwdmodel']['output_omega']['stop']
                except:
                    logger.info('Invalid omega start, stop or end parameters. Defaulting to the range of 0 to 360 degrees with 0.1 degree step size.')
                    omega_start = 0.0
                    omega_step = 0.1
                    omega_stop = 360.0

                # Apply Gaussian blur to the spots?
                try:
                    ge2_blur_sigma = cfg.get('forward_modeling')['fwdmodel']['ge2_blur_sigma']
                except:
                    ge2_blur_sigma = 3

                logger.info('Writing GE2 output to %s', output_ge2)
                # Actually write to a GE2.
                if fwdmodel_mode is 'strainonly':
                    logger.info('Starting GE2 output for strainonly mode')
	            ms.write_xyo_to_ge2(output_ge2=output_ge2, 
                                        omega_start=omega_start, 
                                        omega_step=omega_step, 
                                        omega_stop=omega_stop, 
                                        ge2_blur_sigma=ge2_blur_sigma)
                else:
                    if hedm_mode == 'nf':
                        logger.info('Starting GE2 output for near-field data')
                        ms.write_xyo_to_ge2_v2(output_ge2=output_ge2,
                                               omega_start=omega_start,
                                               omega_step=omega_step,
                                               omega_stop=omega_stop,
                                               ge2_blur_sigma=ge2_blur_sigma,
                                               pixel_size=[0.00148, 0.00148])
                    else:
                        logger.info('Starting GE2 output for far-field data')
                        ms.write_xyo_to_ge2_v2(output_ge2=output_ge2, 
                                               omega_start=omega_start, 
                                               omega_step=omega_step, 
                                               omega_stop=omega_stop, 
                                               ge2_blur_sigma=ge2_blur_sigma) 

            if output_tif_flag is not False:
                # Write output to TIF. First get the file name.
                try:
                    output_tif_stem = cfg.get('forward_modeling')['fwdmodel']['output_tif_stem']
                except:
                    logger.info('Invalid name for output TIF file. Defaulting to nf_%05d.tif.')
                    output_tif_stem = 'nf'
                # Get omega angular bounds (degrees)
                try:
                    omega_start = cfg.get('forward_modeling')['fwdmodel']['output_omega']['start']
                    omega_step = cfg.get('forward_modeling')['fwdmodel']['output_omega']['step']
                    omega_stop = cfg.get('forward_modeling')['fwdmodel']['output_omega']['stop']
                except:
                    logger.info('Invalid omega start, stop or end parameters. Defaulting to the range of 0 to 360 degrees with 0.1 degree step size.')
                    omega_start = 0.0
                    omega_step = 0.1
                    omega_stop = 360.0

                # Apply Gaussian blur to the spots?
                try:
                    tif_blur_sigma = cfg.get('forward_modeling')['fwdmodel']['tif_blur_sigma']
                except:
                    tif_blur_sigma = 3
                logger.info('Writing GE2 output to %s', output_ge2)
                # Actually write to TIF.
                if fwdmodel_mode is 'strainonly':
                   logger.info('Not implemented yet.')
                else:
                    if hedm_mode is 'nf':
                        ms.write_xyo_to_ge2_v2(output_ge2=output_tif_stem,
                                               omega_start=omega_start,
                                               omega_step=omega_step,
                                               omega_stop=omega_stop,
                                               ge2_blur_sigma=tif_blur_sigma,
                                               pixel_size=[0.00148, 0.00148],
                                               is_tif=True)
                    else:
                        ms.write_xyo_to_ge2_v2(output_ge2=output_ge2,
                                               omega_start=omega_start,
                                               omega_step=omega_step,
                                               omega_stop=omega_stop,
                                               ge2_blur_sigma=ge2_blur_sigma)

                logger.info('Forward modeling COMPLETED')
                # DONE
        else:
            logger.error('Invalid forward modeling mode: %s. Choices are datagen and fwdmodel', fwd_model_mode)
