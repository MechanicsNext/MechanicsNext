###############################################################################
# Diffraction-toolkit (https://github.com/hmparanjape/diffraction-toolkit)
# GE Pre-processor Library
#
# Find local maxima of spots in HEDM diffraction patterns.
# Return spot shape, size data.
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
from ge_processor.ge_pre_processor import *

if __name__ == '__main__':
    # Read args
    if len(sys.argv) < 2:
        print 'USAGE: python spot_maxima_finder.py config.yml'
        sys.exit(1)
    elif sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print 'USAGE: python spot_maxima_finder.py config.yml'
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
    	# Initialize the GE pre-processor
    	gepp = GEPreProcessor(cfg=cfg, logger=logger)
    	# Start analysis
    	logger.info('=== begin image-smoothing ===')
    	# Load the GE2 data
    	gepp.load_data()
    	# ID blobs and local maxima, and optionally print a text file with local maxima positions,
        # PNG files and GE2 file with cleaned-up data.
    	gepp.find_blobs()

