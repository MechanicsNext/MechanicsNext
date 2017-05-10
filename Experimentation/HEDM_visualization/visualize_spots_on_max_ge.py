# Superimpose simulated ff-HEDM spots on a GE2/GE3 file
# Workflow:
# 1. Generate a max-over GE2 file using the max-over function in 
#    ../HEDM_file_utils/bash-hedm-utils.sh
# 2. Run a ff-HEDM analysis using heXRD.
# 3. Run ../HEDM_forward_modeling/synth_fwd_modeling_patterns.py to
#    simulate spot positions. You will need 'accepted_orientations.dat'
#    and instrucment, material config files.
# 4. Run this script to generate a PNG file with the intensity from GE2 file
#    and spots from simulated data.
# Usage:
#    python visualize_spots_on_max_ge.py GE_FILE.GE2 SIM_SPOT_FILE.out
#
# Harshad Paranjape (hparanja@mines.edu)
# MechanicsNext (https://github.com/MechanicsNext/MechanicsNext)
#
import copy
import logging
import os
import sys
import time
import warnings
# Pickle saves a data structure to a binary file
try:
   import cPickle as pickle
except:
   import pickle
import yaml
# Numpy and Scipy for array operations, image processing etc.
import numpy as np
import scipy.ndimage as ndimage
from scipy.ndimage.filters import gaussian_filter
# Matplotlib plots like Matlab
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
# hexrd helpers to read a config file and load GE2 data
from hexrd import config
from hexrd.coreutil import initialize_experiment
# Additional image analysis routines
from skimage.morphology import watershed
from skimage.feature import peak_local_max
# Data (spot) clustering
from sklearn.cluster import DBSCAN
# Parallelization for speed
from joblib import Parallel, delayed
import multiprocessing
from random import randint
#
# Helper to save a 2D array as an image (Contributed by Branden Kappes @ mines.edu)
def write_image(filename, arr, pts=None, minsize=None, **kwds):
    '''
    Write a 2D array (arr) to a PNG image file. Optionally superimpose
    points (pts) that are described in terms of their x, y coordinates.
    '''
    xsize, ysize = 1024, 1024

    fig = plt.figure(figsize=(20, 20), dpi=120)
    # plot the image
    ax = fig.gca()
    kwds['interpolation'] = kwds.get('interpolation', 'none')
    image = ax.imshow(arr, **kwds)
    ax.set_xlabel(r'X', fontsize='large')
    ax.set_ylabel(r'Y', fontsize='large')
    cbar = plt.colorbar(image)
    # plot any points
    if pts is not None:
        pts = np.asarray(pts)
        #ax.plot(pts[:,1], pts[:,0], 'go', markersize=9)
        ax.scatter(pts[:,1], pts[:,0], s=81, facecolors='none', edgecolors='y')
        # resize (since adding points often adds padding)
        ax.set_xlim(0, 2048)
        ax.set_ylim(0, 2048)
    fig.savefig(filename, bbox_inches='tight', pad_inches=1./3.)
    fig.clf()
    plt.close()
#--
def write_ge2(filename, arr, nbytes_header=8192, pixel_type=np.uint16):
    '''
    Write a 3D array to a single GE2 file.
    '''
    fid = open(filename, 'wb')
    fid.seek(nbytes_header)
    fid.write(arr.astype(pixel_type))
    fid.close()
#--

if __name__ == '__main__':

    if len(sys.argv) < 3:
        print 'USAGE: python visualize_spots_on_max ge_file.ge2 spot_data.out'
        sys.exit()

    ge_filename = sys.argv[1]
    spot_filename = sys.argv[2]

    try:
        image_filename = sys.argv[3]
        image_filename = image_filename.strip().lower()
    except:
        image_filename = 'spots.png'

    try:
        intensity_min = float(sys.argv[4])
    except:
        intensity_min = 0

    try:
        intensity_max = float(sys.argv[5])
    except:
        intensity_max = 1000

    f = open(ge_filename, 'rb')
    f.seek(8192)
    ge_data = np.fromfile(f, dtype=np.uint16).reshape((2048, 2048))
    f.close()

    spot_data = np.loadtxt(spot_filename)

    write_image(image_filename, ge_data, pts=spot_data[:, 0:2], vmin=intensity_min, vmax=intensity_max)

