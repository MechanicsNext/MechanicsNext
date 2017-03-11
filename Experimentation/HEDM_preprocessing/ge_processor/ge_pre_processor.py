# (https://github.com/MechanicsNext/MechanicsNext/blob/master/Experimentation/HEDM_preprocessing/)
#
# Parent: ../spot_maxima_finder.py
#
# Find local maxima of spots in ff-HEDM diffraction patterns.
# Return spot position, shape, and size data.
#
# USAGE:
#  See ../spot_maxima_finder.py
#
# See ../examples_spot_maxima_finder/ folder for example configurations.
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
        ax.plot(pts[:,1], pts[:,0], 'go', markersize=6)
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
def find_blobs_mp(ge_data, int_scale_factor, min_size, min_peak_separation, cfg):
   '''
   Multiprocessing worker subroutine for blob and local maxima finding.
   First finds blobs using a connected component type algorithm. Then
   finds the local maxima (spots). Then calculates the area corresponding to 
   each of the spots using the watershed algorithm.
   '''
   # Remove noise by masking low intensity points
   ge_mask = ge_data > (int_scale_factor * cfg.get('pre_processing')['ge_reader_threshold'])
   # Dilate the mask to merge any small spots
   ge_mask = ndimage.binary_dilation(ge_mask, ndimage.generate_binary_structure(3, 3))
   # Label each connected component (blob)
   label_ge, number_of_labels = ndimage.label(ge_mask)
   ge_labeled = np.amax(label_ge, 1)
   label_rand = str(randint(0, 1000000))
   # Get size of each blob by counting number of pixels with the same label
   blob_sizes = ndimage.sum(ge_mask, label_ge, range(1, number_of_labels + 1))
   blob_labels = np.unique(label_ge)
   # Loop over all detected regions and filter based on size/aspect ratio
   # Save filtered blob information in an array (blobs) of GEBlob object type
   blobs = []                     # Blobs
   blob_centroids = []            # Blob position (centroid)
   max_points_global = []         # Local maxima in global coordinates
   roi_global = []                # Region of interest for each blob
   roi_maxima_global = []         # Local maxima inside each ROI
   watershed_global = []          # Result of watershed segmentation
   watershed_pixel_count = []     # Size of each watershed region (spot)
   watershed_eigs = []            # Size of each spot (eigenvalues)
   watershed_eigv = []            # Principal axes of the segmented spots
   for label_num, blob_size in zip(blob_labels, blob_sizes):
        # Label 0 is for the whole image background. We do not want that
        if label_num == 0:
            continue
        # Are total pixels in a blob < min_size? Then move on
        if blob_size < min_size:
            continue
        # Get the minimal region of interest (ROI) for the blob
        slice_x, slice_y, slice_z = ndimage.find_objects(label_ge == label_num)[0]
        # Save the bounding box dimensions of the roi.
        bbox = [slice_x.stop -  slice_x.start, slice_y.stop -  slice_y.start, slice_z.stop -  slice_z.start]
        # Is any of roi bounding box dim < min? Move on
        if(min(bbox) < (min_size ** (1.0/3.0))):
            continue
        # Calculate the blob centroid coordinates in global frame
        centroid_tmp = [(slice_x.stop + slice_x.start)/2.0, (slice_y.stop + slice_y.start)/2.0 - 1024.0, (slice_z.stop + slice_z.start)/2.0 - 1024.0]
        # If the radial distance of the blob centroid from the detector center is greater than a threshold, then do not process that blob
        # This option helps restrict analysis to a few inner rings in the data. Outer rings tend to have too numerous spots
        if np.sqrt(np.power(centroid_tmp[1], 2.0) + np.power(centroid_tmp[2], 2.0)) > cfg.get('pre_processing')['radial_threshold']:
            continue
        # This blob has passed several criteria. Save its centroid.
        blob_centroids.append([(slice_x.stop + slice_x.start)/2.0, (slice_y.stop + slice_y.start)/2.0, (slice_z.stop + slice_z.start)/2.0])
        #
        # Now run local maxima finding and then watershed
        # Get data only inside an roi.
        roi = ge_data[slice_x, slice_y, slice_z]
        roi_original = ge_data[slice_x, slice_y, slice_z]
        # 
        markers = np.zeros_like(roi)
        # Find local maxima
        max_points = peak_local_max(roi, min_distance=(min_peak_separation),
                                    threshold_rel=0.01, exclude_border=False, indices=False)
        # If a maximum is weak in intensity, delete it.
        max_points[roi < int_scale_factor * cfg.get('pre_processing')['ge_reader_threshold']] = 0
        # Get local coordinates of the maxima (local to the roi boundig box)
        max_points = np.nonzero(max_points)
        # Create a marker array (3D) that has the local maxima assigned unique IDs.
        # Later we will run watershed starting from these markers
        for max_x, max_y, max_z, max_id in zip(max_points[0], max_points[1], max_points[2],
                                               range(len(max_points[0]))):
           markers[max_x][max_y][max_z] = max_id+1
           # Store global coordinates and intensity of the maxima
           max_points_global.append([max_x + slice_x.start, max_y + slice_y.start, max_z + slice_z.start, roi[max_x][max_y][max_z]])
        # Run watershed now
        labels = watershed(-roi, markers, mask=(roi>0.1*np.amax(roi)))
        roi_global.append(roi)
        # Save watershed results
        watershed_global.append(labels)
        roi_maxima = np.zeros_like(roi)
        # Get spot size and shape
        for max_x, max_y, max_z, max_id in zip(max_points[0], max_points[1], max_points[2],
                                               range(len(max_points[0]))):
           # Size = number of voxels with same watershed label
           watershed_pixel_count.append(np.sum(labels == max_id+1))
           # Save the position of local maxima
           roi_maxima[max_x][max_y][max_z] = 1
           # Shape = three eigenvalues (calculated using SVD)
           try:
              # The indices must be centered around the mean
              roi_pt_indices = np.transpose(np.nonzero(labels == max_id+1))
              roi_pt_indices = roi_pt_indices - np.mean(roi_pt_indices, axis=0)
              U, S, V = np.linalg.svd(roi_pt_indices, compute_uv=True)
              # Instead of using the eigenvalues, which are shorter than axes lengths
              # for discrete point clouds, take the lenght of each axis as the largest
              # projection from a point in the point cloud to the eigenvector for that
              # axis
              # Size = projection of roi points on the eigenvectors
              watershed_eigs.append(np.max(np.dot(roi_pt_indices, V), axis=0))
              watershed_eigv.append(np.reshape(V, 9))
           except:
              print 'Something went wrong with PCA. Setting eigenvalues to 0.'
              watershed_eigs.append([0., 0., 0.])
	      watershed_eigv.append([0., 0., 0., 0., 0., 0., 0., 0., 0.])

        # Done watershed
        roi_maxima_global.append(roi_maxima)

        # This looks like a legitimate spot. Add to blobs array
        blobs.append(GEBlob(slice_x, slice_y, slice_z, label_num, blob_size, max_points))

   ge_labeled = np.amax(label_ge, 0)

   # Return various quantities calculated/measured for each spot
   return {'blobs': blobs, 
           'label_ge': label_ge, 
           'blob_centroids': blob_centroids, 
           'local_maxima': max_points_global, 
           'roi' : roi_global, 
           'watershed' : watershed_global,
           'watershed_pixel_count': watershed_pixel_count,
           'watershed_eigs': watershed_eigs,
	   'watershed_eigv': watershed_eigv,
           'local_maxima':roi_maxima_global}
#--

# A blob is a set of pixels in an image that are connected to each other and
# have intensity above a threshold
class GEBlob:
    def __init__(self, slice_x, slice_y, slice_z,
                 blob_label, blob_size, max_points):
       # slice_* defines the bounding box of a blob in the GE2 data
       self.slice_x    = slice_x
       self.slice_y    = slice_y
       self.slice_z    = slice_z
       # Blob id
       self.blob_label = blob_label
       # Number of pixels
       self.blob_size  = blob_size
       # Local maxima in the blob ([[x], [y], [z]])
       # xyz for the local maxima are w.r.t to the blob
       # To get global xyz for the maxima, add slice_x.start
       # etc. to the xyz
       self.max_points = max_points
#--
# An object for all the GE2 pre-processing routines
class GEPreProcessor:
    '''
        Pre-processing on GE files to extract blobs,
        detect local maxima etc.
    '''

    def __init__(self, cfg, logger):
    	self.cfg                 = cfg                      # An open hexrd config file object
        self.logger              = logger                   # An open logger object
    	self.ge_data             = []                       # GE2 image data
    	self.ge_smooth_data      = []                       # Above + smoothed using Gauss
        self.ge_labeled_data     = []                       # Above + connected components labeled
        self.number_of_labels    = 0                        # number of labels = number of blobs in the data
        self.min_blob_size       = cfg.get('pre_processing')['min_blob_size']       # Blobs smaller than this are removed (user input)
        self.int_scale_factor    = 1                                                # A scale factor for intensity (auto-calculated)
        self.min_peak_separation = cfg.get('pre_processing')['min_peak_separation'] # Minimum separation in the local maxima (user input)
        self.blobs               = []                       # An array of blob objects
	self.max_points          = []                       # An array of local maxima coordinates in the blobs
        self.omega_start         = []                       # In the parallelized frame data, start omega number for each portion

        return
    #--
    def load_data(self):
        '''
            Read the config file and load appropriate GE2
            frames.
        '''
        cfg          = self.cfg
        logger       = self.logger
        omega_start  = self.omega_start
        # process the data
        pd, reader, detector = initialize_experiment(cfg)
        n_frames = reader.getNFrames()
        logger.info("Reading %d frames of data, storing values > %.1f",
                    n_frames, cfg.get('pre_processing')['ge_reader_threshold'])
        # Loop over all frames and save them in a 3D array
        # The first array dimension in frame_list is along omega
        frame_list = []
        for i in range(n_frames):
            frame = reader.read()
            frame_list.append(frame)
        # Turn the frame array into a Numpy array
        frame_list = np.array(frame_list)
        # Remove low intensity noise
        frame_list[frame_list < cfg.get('pre_processing')['ge_reader_threshold']] = 0
        # Scale the intensity to 16000
        int_scale_factor = float(2**14)/float(np.amax(frame_list))
        frame_list = frame_list*int_scale_factor
        if cfg.get('pre_processing')['print_diag_images']:
        	# Flatten along omega and write the frame array to an image
        	write_image('slice.png', np.amax(frame_list, axis=0), vmin=0)
        # Split the frame array into chunks for multiprocessing
        # Right now I get the number of processors and split data into those many chunks
        num_cores = cfg.multiprocessing
        # But we don't want to split a small dataset into way too many chunks.
        # Otherwise too many spots may get split along chunk boundaries.
        num_cores = np.round(np.min([num_cores, np.round(np.shape(frame_list)[0]/60.)]))

        frame_list_split = np.array_split(frame_list, num_cores, axis=0)
        ge_data_ang_red = ()
        omega_start.append(0)
        for array_piece in frame_list_split:
           ge_data_ang_red = ge_data_ang_red + (array_piece,)
           omega_start.append(np.shape(array_piece)[0])
        # Save the omega values at chunk boundaries
        omega_start.pop()
        omega_start = np.cumsum(omega_start)
        logger.info("Finished reading frames")

        logger.info("Split data for parallel processing at omega = %s", omega_start.tostring())

        self.ge_data          = frame_list
        self.int_scale_factor = int_scale_factor
	self.ge_data_ang_red  = ge_data_ang_red
        self.omega_start      = omega_start
        self.input_data_shape = np.shape(frame_list)

        return frame_list
    #--
    def find_blobs(self):
        '''
            Find connected componnets (blobs). Then process blobs to get local maxima.
        '''
        min_size            = self.min_blob_size
        ge_data_smooth      = self.ge_smooth_data
        cfg                 = self.cfg
        logger              = self.logger
        int_scale_factor    = self.int_scale_factor
        ge_data_ang_red     = self.ge_data_ang_red
        min_peak_separation = self.min_peak_separation
        omega_start         = self.omega_start
        # Process data chunks in parallel
        num_cores = cfg.multiprocessing
        logger.info("Starting spot finding with %d cores", num_cores)
        blobs_mp_output = Parallel(n_jobs=num_cores, verbose=5, max_nbytes=1e6)(delayed(find_blobs_mp)(ge_data, int_scale_factor, min_size, min_peak_separation, cfg) for ge_data in ge_data_ang_red)

        logger.info("Finished multiprocessing spot finding algorithm")

        blobs = []
        label_ge = []
        blob_centroids_oxy = []
        local_maxima_oxy = []
        local_maxima_oxyi = []
        spot_sizes = []
        spot_shapes = []
        spot_axes = []
        roi_counter = 0

        # Process the combined output for all chunks
        for blobs_mp_output_i, omega_start_i in zip(blobs_mp_output, omega_start):
           # Consolidate all blobs in a single array.
           for blob_i in blobs_mp_output_i['blobs']:
              blobs.append(blob_i)
           # Get roi, local maxima, and watershed information for each blob.
           for roi_i, watershed_i, local_maxima_i in zip(blobs_mp_output_i['roi'], blobs_mp_output_i['watershed'], blobs_mp_output_i['local_maxima']):
              # Print diagnostic images for roi, watershed etc. if requested in the config file.
              if cfg.get('pre_processing')['print_diag_images']:
                 write_image('watershed' + str(roi_counter) + '.png', np.amax(watershed_i, axis=0), vmin=0)
                 write_image('roi' + str(roi_counter) + '.png', np.amax(roi_i, axis=0), vmin=0)
                 write_image('local_max' + str(roi_counter) + '.png', np.amax(local_maxima_i, axis=0), vmin=0)
                 roi_counter += 1
           # Store coordinates, intensity, and spot shape/size.   
           for maxima_info, spot_size_i, spot_eigs_i, spot_eigv_i in zip(blobs_mp_output_i['local_maxima'], blobs_mp_output_i['watershed_pixel_count'], blobs_mp_output_i['watershed_eigs'], blobs_mp_output_i['watershed_eigv']):
              maxima_o, maxima_x, maxima_y, max_intensity = maxima_info
              if max_intensity > (cfg.get('pre_processing')['ge_reader_threshold']):
                 local_maxima_oxyi.append([maxima_o + omega_start_i, maxima_x, maxima_y, max_intensity])
                 local_maxima_oxy.append([maxima_o + omega_start_i, maxima_x, maxima_y])
                 spot_sizes.append(spot_size_i)
                 logger.info("Blob data: %5d %5d %5d %5d", maxima_o, omega_start_i, maxima_x, maxima_y)
                 # Consolidate array shapes
                 try:
                    spot_shapes.append([spot_eigs_i[0], spot_eigs_i[1], spot_eigs_i[2]])
		    spot_axes.append(spot_eigv_i)
                 except:
                    print 'Eig shape not correct'
                    spot_shapes.append([0., 0., 0.])
		    spot_axes.append([0., 0., 0., 0., 0., 0., 0., 0., 0.])

           #
           label_ge.append(blobs_mp_output_i['label_ge'])
           
        # Cluster spots that are very close to each other
        logger.info("Clustering spots")
        # Cluster the local minima
        local_maxima_oxy = np.array(local_maxima_oxy)
        local_maxima_oxyi = np.array(local_maxima_oxyi)
        # For now, the clustering radius is determined from the mean_peak_separation value
        eps_val = np.max([cfg.get('pre_processing')['min_peak_separation'], 2.0])
        db = DBSCAN(eps=eps_val, min_samples=1).fit(local_maxima_oxy)
        local_maxima_labels = db.labels_
        # Save spot coordinates and intensity after clustering
        o_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 0])
        x_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 1])
        y_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 2])
        i_sum = np.bincount(local_maxima_labels, weights=local_maxima_oxyi[:, 3])
        spot_sizes_sum = np.bincount(local_maxima_labels, weights=spot_sizes)

        spot_shapes = np.array(spot_shapes)
	spot_axes = np.array(spot_axes)
        # Right now I am saving "mean" spot shape and size. This is sketchy and must be revisited.
        spot_shapes_sum = zip(np.bincount(local_maxima_labels, weights=spot_shapes[:, 0]), 
                              np.bincount(local_maxima_labels, weights=spot_shapes[:, 1]),
                              np.bincount(local_maxima_labels, weights=spot_shapes[:, 2]))

        spot_axes_sum = zip(np.bincount(local_maxima_labels, weights=spot_axes[:, 0]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 1]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 2]),
			    np.bincount(local_maxima_labels, weights=spot_axes[:, 3]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 4]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 5]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 6]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 7]),
                            np.bincount(local_maxima_labels, weights=spot_axes[:, 8]))

        l_sum = np.bincount(local_maxima_labels)
        #
        local_maxima_oxyi_clustered = zip(np.divide(o_sum, l_sum), np.divide(x_sum, l_sum), np.divide(y_sum, l_sum), np.divide(i_sum, l_sum))
        local_maxima_xy = zip(np.divide(x_sum, l_sum), np.divide(y_sum, l_sum))

        spot_sizes_clustered = np.divide(spot_sizes_sum, l_sum)
        spot_shapes_clustered = np.divide(spot_shapes_sum, zip(l_sum, l_sum, l_sum))
        spot_axes_clustered = np.divide(spot_axes_sum, zip(l_sum, l_sum, l_sum, l_sum, l_sum, l_sum, l_sum, l_sum, l_sum))

        # Prune the clustered spots for spot size threshold and save final data.
        blob_centroids = []
        for blob in blobs:
           if blob.blob_size > cfg.get('pre_processing')['min_blob_size']:
              blob_centroids.append([(blob.slice_y.start + blob.slice_y.stop)/2.0, (blob.slice_z.start + blob.slice_z.stop)/2.0])
              blob_centroids_oxy.append([(blob.slice_x.start + blob.slice_x.stop)/2.0, (blob.slice_y.start + blob.slice_y.stop)/2.0, (blob.slice_z.start + blob.slice_z.stop)/2.0])

        self.blobs                       = blobs
        self.blob_centroids_oxy          = blob_centroids_oxy
        self.local_maxima_oxyi           = local_maxima_oxyi
        self.local_maxima_oxyi_clustered = local_maxima_oxyi_clustered
        self.local_maxima_clusters       = db

        logger.info("Finished clustering")
        logger.info("Found %d connected components", len(blobs))
        logger.info("Found %d spots", np.shape(local_maxima_oxyi)[0])

        # Before results are returned, an ooportunity to scale the intensities
        # The factor of 5 is because the uint16 format can store intensities upto 60k and not just 12k.
        max_vals = np.max(self.local_maxima_oxyi_clustered, axis=0)
        max_o, max_x, max_y, max_i = max_vals
        intensity_scale_factor = 12000.0 / max_i * 5.0

        # Write the data to text file and diagnostic images.
	if cfg.get('pre_processing')['print_diag_images']:
        	logger.info("Writing diagnostic images")
        	# Superimpose the centroids of blobs and the local maxima on the original data 
        	# and write to an image file.
       		write_image('slice_blob_centroids.png', np.amax(self.ge_data, axis=0), pts=blob_centroids)
        	write_image('slice_local_maxima.png', np.amax(self.ge_data, axis=0), pts=local_maxima_xy)

        # Print spot data to a text file
        if(cfg.get('pre_processing')['print_spots_info']):
	   f = open(cfg.get('analysis_name') +'_spots.data', 'w+')
           # Because I want to pretty print columns...
           template = "{0:>12}{1:>12}{2:>12}{3:>12}{4:>12}{5:>12}{6:>12}{7:>12}{8:>12}{9:>12}{10:>12}{11:>12}{12:>12}{13:>12}{14:>12}{15:>12}{16:>12}"
           template_file = "{0:>12}{1:>12}{2:>12}{3:>12}{4:>12}{5:>12}{6:>12}{7:>12}{8:>12}{9:>12}{10:>12}{11:>12}{12:>12}{13:>12}{14:>12}{15:>12}{16:>12}\n"
           print template.format("Omega", "X", "Y", "Intensity", "Size", "Eig 1", "Eig 2", "Eig 3", "Eigv 1", "Eigv 2", "Eigv 3", "Eigv 4", "Eigv 5", "Eigv 6", "Eigv 7", "Eigv 8", "Eigv 9")
	   f.write(template_file.format("Omega", "X", "Y", "Intensity", "Size", "Eig 1", "Eig 2", "Eig 3", "Eigv 1", "Eigv 2", "Eigv 3", "Eigv 4", "Eigv 5", "Eigv 6", "Eigv 7", "Eigv 8", "Eigv 9"))

           template = "{0:12.2f}{1:12.2f}{2:12.2f}{3:12.3f}{4:12.3f}{5:12.3f}{6:12.3f}{7:12.3f}{8:12.3f}{9:12.3f}{10:12.3f}{11:12.3f}{12:12.3f}{13:12.3f}{14:12.3f}{15:12.3f}{16:12.3f}"
           template_file = "{0:12.2f}{1:12.2f}{2:12.2f}{3:12.3f}{4:12.3f}{5:12.3f}{6:12.3f}{7:12.3f}{8:12.3f}{9:12.3f}{10:12.3f}{11:12.3f}{12:12.3f}{13:12.3f}{14:12.3f}{15:12.3f}{16:12.3f}\n"
           for clustered_maxima_info, spot_size, spot_shape, spot_axes in zip(local_maxima_oxyi_clustered, spot_sizes_clustered, spot_shapes_clustered, spot_axes_clustered):
              o, x, y, i = clustered_maxima_info
              print template.format(o + 1, y, x, i * intensity_scale_factor, spot_size, spot_shape[0], spot_shape[1], spot_shape[2], 
				    spot_axes[0], spot_axes[1], spot_axes[2], spot_axes[3], spot_axes[4], spot_axes[5], spot_axes[6], spot_axes[7], spot_axes[8])
	      f.write(template_file.format(o + 1, y, x, i * intensity_scale_factor, spot_size, spot_shape[0], spot_shape[1], spot_shape[2],
					   spot_axes[0], spot_axes[1], spot_axes[2], spot_axes[3], spot_axes[4], spot_axes[5], spot_axes[6], spot_axes[7], spot_axes[8]))

	   f.close()
        # Write clean data to a GE2 file
        if cfg.get('pre_processing')['print_ge']:
           logger.info("Writing GE2 files")
           # Synthesize a GE2 file based on the IDed spots. Currently this is bit of a performance bottleneck, but not too bad.
           frames_synth = np.zeros(self.input_data_shape)
           intensity_threshold = cfg.get('pre_processing')['upper_intensity_threshold']
           for o, x, y, i in local_maxima_oxyi_clustered:
              if i < intensity_threshold:
                  frames_synth[int(round(o)), int(round(x)), int(round(y))] = i * intensity_scale_factor

           # Optionally dilate spots. If spots are dilated, a Gaussian smoothing is also done.
           # Generally it's a good idea to apply 1 to 2 pixel dilation.
           grey_dilation_size = (cfg.get('pre_processing')['radius_gray_dilation_omega'],
				 cfg.get('pre_processing')['radius_gray_dilation_x'],
				 cfg.get('pre_processing')['radius_gray_dilation_y'])
           frames_synth_grey  = ndimage.morphology.grey_dilation(frames_synth, size=grey_dilation_size)
           frames_synth_gauss = ndimage.filters.gaussian_filter(frames_synth_grey, sigma=grey_dilation_size)
           write_ge2(cfg.get('analysis_name') + '_synth_spots.ge2', frames_synth_gauss)
        else:
           logger.info("Skipped writing GE2 files")

        return label_ge
    #--
