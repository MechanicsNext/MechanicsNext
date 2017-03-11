###############################################################################
# MechanicsNext: ff-HEDM Spot Segmentation Program
# (https://github.com/MechanicsNext/MechanicsNext/blob/master/Experimentation/HEDM_preprocessing/)
#
# Find local maxima of spots in ff-HEDM diffraction patterns.
# Return spot position, shape, and size data.
#
# USAGE:
#  python spot_maxima_finder.py config_file.yml
#
# See examples_spot_maxima_finder/ folder for example configurations.
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
# YAML and heXRD imports to read heXRD-like configuration files
import yaml
from hexrd import config
# Main program routines
from ge_processor.ge_pre_processor import *

# Run the program.
if __name__ == '__main__':
    # Read command line arguments
    if len(sys.argv) < 2:
        print 'USAGE: python spot_maxima_finder.py config.yml'
        sys.exit(1)
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print 'USAGE: python spot_maxima_finder.py config.yml'
        sys.exit(1)
    else:
        cfg_file = sys.argv[1]

    # Setup logger. Information about the progress of the program will be
    # printed to the logger.
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

    # cfg is a list. We will loop over each cfg set. i.e., we can in principal have multiple
    # spot cleanup jobs in a single configuration file.
    for cfg in cfgs:
    	# Initialize the GE pre-processor object. This object wraps all
        # spot information.
    	gepp = GEPreProcessor(cfg=cfg, logger=logger)
    	# Start analysis
    	logger.info('=== begin image-smoothing ===')
    	# Load the GE2 data. This will take several minutes, since a typical ff-HEDM
        # dataset consists of tens of GB of data.
    	gepp.load_data()
        # For spot segmentation, following terminology is used:
        # Raw data or GE2 data: 3D array consisting of X, Y, and Omega dimensions that stores intensity information.
        # Blob: Any volume in the raw data that has intensity above a threshold. This may contain several overlapping spots.
        # Local maxima: Local intensity maxima in a blob.
        # Spot: A region inside a blob that has only one intensity maximum inside it.
        #
    	# Identify blobs and local maxima, and optionally print:
        # 1. A text file with local maxima positions,
        # 2. PNG files with diagnostic images showing blobs, local maxima positions, and spots identified using the watershed algo.
        # 3. GE2 file with cleaned-up (segmented) data.
    	gepp.find_blobs()
    # END of program

