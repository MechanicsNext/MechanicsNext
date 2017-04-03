# A brute-force detector parameter study for near-field HEDM calibration
# Based on the scripts by Darren Pagan and heXRD developers.
# With additions by Ashley Bucsek.
#
# USAGE:
#   python nf_detector_calibration_bruteforce.py
#
# Refines distance, x_cen, ytilt, and ztilt by picking four random values
# for those parameters between two bounds. The bounds are user specified (see below).
# After using a parameter set, a confidance map is calculated for
# a coarse grid and then the number of points with the confidence
# greater than a cutoff is written two a text file along with
# the parameters used. The confidence cutoff and the name of the 
# text output file is user specified.
# All user options are specified at the start of the script below.
#
# Harshad Paranjape (hparanja@mines.edu)
# MechanicsNext
#
import yaml
from skimage import io as imageio
from sklearn.utils.extmath import cartesian  # For making permutations of all variable detector parameters
import copy
import cPickle as cpl
import time
import numpy as np
from scipy import sparse
import scipy.ndimage as img
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.widgets import Slider, Button, RadioButtons
import multiprocessing as mp

from hexrd import matrixutil as mutil
from hexrd.gridutil import cellIndices, cellConnectivity, cellCentroids, compute_areas
from hexrd.xrd import rotations as rot
from hexrd.xrd import material
from hexrd.xrd.crystallography import processWavelength
from hexrd.xrd.symmetry import applySym
from hexrd.xrd import transforms as xf
from hexrd.xrd.transforms_CAPI import anglesToGVec, \
                                      makeRotMatOfExpMap, makeDetectorRotMat, makeOscillRotMat, \
                                      gvecToDetectorXY, detectorXYToGvec
from hexrd.xrd.xrdutil import angularPixelSize, make_reflection_patches, \
                              simulateGVecs, _project_on_detector_plane
from hexrd import valunits

from progressbar import ProgressBar, Percentage, Bar

#
#==============================================================================
# %% ALL USER INPUTS
#==============================================================================
grain_out_file     = '/FULL_PATH_TO_FF_GRAINS_OUTPUT_FROM_HEXRD/grains.out'
det_file           = '/FULL_PATH_TO_DETECTOR_FILE/retiga_Dec15.yml'
mat_file_loc       = '/FULL_PATH_TO_THE_MATERIAL_FILE/material_filename_61keV' 
data_folder        = '/LOCATION_OF_THE_NF_DATA_FOLDER/nf/'

img_start          = 10701                      # Number of the first NF image
num_imgs           = 721                        # Number of images
#
x_ray_energy       = 61.332                     # keV
mat_name           = 'multigold61'              # Material name

ome_period_deg     = (0.,360.)                  # degrees
ome_range_deg      = [(0.,360.)]                # degrees
ome_step_deg       = 0.5                        # degrees 

max_tth            = 10.                        # degrees. if a -ve number is input, all peaks that will hit the detector are calculated
struct_fact_cutoff = 0.                         # values from 0 to 1 for cutting off weak peaks, use 0. unless you know what you're doing 
#horizontal plane (x-z) reconstruction voxel bounds in microns
h_bnds             = [-650.,650.]
h_spacing          = 10.
#vertical (y) reconstruction voxel bounds #microns
v_bnds             = [0.,0.]                    # went up to 160
v_spacing          = 1.

#image processing
num_for_median     = 250                        # num images to use for median data
threshold          = 0.
num_erosions       = 3                          # num iterations of images erosion, don't mess with unless you know what you're doing
num_dilations      = 2                          # num iterations of images erosion, don't mess with unless you know what you're doing
ome_dilation_iter  = 1                          # num iterations of 3d image stack dilations, don't mess with unless you know what you're doing

chunksize          = 200                        # chunksize for multiprocessing, don't mess with unless you know what you're doing
#
# OPTIONS FOR CALIBRATION PARAMETER STUDY
# Bounds and step for the parameter study
# syntax: parameter = [min_value, max_value, step]
# e.g. distance = [-8.0, -7.0, 0.1] will vary the distance from -8 to -7.1.
# In this bruteforce script, the step is ignored. Rather, parameter values are
# set to a random number between the bounds.
distance           = [-35  , -15  , 0.1 ]       # mm [typically [7,8]]
x_cen              = [-2   , 2    , 0.1 ]       # mm typically [-.1,.1]
ytilt              = [-0.5 , 0.5  , 0.1 ]       # input in degrees typically [-.5,.5]. NOT RADIANS
ztilt              = [-0.5 , 0.5  , 0.1 ]       # input in degrees typically [-.5,.5]. NOT RADIANS.
number_of_tests    = 2000                       # Number of four-parameter sets to study.
#
confidence_thresh  = 0.5                        # Confidence threshold. Grid points with the conf above this will be counted.
#
# Output of the parameter study will be saved to this text file.
# First four columns: distance, x_cen, ytilt, ztilt respectively. Fifth column: max confidence. Sixth column: Points with confidence
# greater than confidence_thresh
paramstudy_output_file = 'nf_param_study_results.data'
#
USE_NUMBA   = True

ref_vparams = np.r_[1., 1., 1., 0., 0., 0.]
ref_gparams = np.r_[0., 0., 0., 1., 1., 1., 0., 0., 0.]

ncpus = mp.cpu_count()

#==============================================================================
# %% ORIENTATION TESTING FUNCTIONS
#==============================================================================
def testorientation_init(params):
    global paramMP
    paramMP = params

def testorientation_reduced(params_in):
    """
    input parameters are [hkl_id, com_ome, com_eta]
    """
    rMat_d       = paramMP['rMat_d']
    chi	      = paramMP['chi']
    tVec_d	      = paramMP['tVec_d']
    tVec_s	      = paramMP['tVec_s']
    distortion   = paramMP['distortion']
    panel_dims   = paramMP['panel_dims']
    panel_buffer = paramMP['panel_buffer']
    x_col_edges  = paramMP['x_col_edges']
    y_row_edges  = paramMP['y_row_edges']
    ome_edges    = paramMP['ome_edges']
    i_dil	      = paramMP['i_dil']
    j_dil	      = paramMP['j_dil']
    nrows	      = paramMP['nrows']
    ncols	      = paramMP['ncols']
    image_stack  = paramMP['image_stack']
#    pd           = paramMP['pd']
#    detector_params = paramMP['detector_params']
#    pixel_size   = paramMP['pixel_size']
#    ome_range    = paramMP['ome_range']
#    ome_period   = paramMP['ome_period']
    all_angles   = paramMP['all_angles']

    #rMat_c = makeRotMatOfExpMap(params_in[0, :])
    #test_crds = params_in[1, :]
    #all_angles = params_in[2:, :]

    exp_maps = params_in[:3]
    rMat_c = makeRotMatOfExpMap(exp_maps)    
    test_crds = params_in[3:6]
#    gparams = np.hstack([params_in[:6], [1., 1., 1., 0., 0., 0.]])
    #sim_results = simulateGVecs(pd,
#                                detector_params,
#                                gparams,
#                                panel_dims=panel_dims,
#                                pixel_pitch=pixel_size,
#                                ome_range=ome_range,
#                                ome_period=ome_period,
#                                distortion=distortion)
    angles_cur = all_angles[int(params_in[6])]
    
    det_xy, rMat_ss = _project_on_detector_plane(angles_cur ,
                                                 rMat_d, rMat_c, chi,
                                                 tVec_d, test_crds, tVec_s, 
                                                 distortion)

    # find on spatial extent of detector
    xTest = np.logical_and(det_xy[:, 0] >= panel_dims[0][0] + panel_buffer,
                           det_xy[:, 0] <= panel_dims[1][0] - panel_buffer)
    #yTest = np.logical_and(det_xy[:, 1] >= panel_dims[0][1] + panel_buffer,
#                           det_xy[:, 1] <= panel_dims[1][1] - panel_buffer)
    
    #kludge to deal with the beam stop
    yTest1 = np.logical_and(det_xy[:, 1] >= panel_dims[0][1] + panel_buffer,
                           det_xy[:, 1] <= -0.27)
    
    yTest2 = np.logical_and(det_xy[:, 1] >= 0.27,
                           det_xy[:, 1] <= panel_dims[1][1] - panel_buffer)
    
    yTest = np.logical_or(yTest1,yTest2)

    onDetector = np.where(np.logical_and(xTest, yTest))[0]
                                                         
    # pick who's valid
    row_indices = cellIndices(y_row_edges, det_xy[onDetector, 1])
    col_indices = cellIndices(x_col_edges, det_xy[onDetector, 0])
    frame_indices = cellIndices(ome_edges, angles_cur[onDetector, 2])

    # perform check
    tmp_confidence = np.zeros(len(frame_indices), dtype=bool)
    for iref, indices in enumerate(zip(frame_indices, row_indices, col_indices)):
        i_sup = indices[1] + i_dil
        j_sup = indices[2] + j_dil      

        idx_mask = np.where(
            np.logical_and(np.logical_and(i_sup >= 0, i_sup < nrows),
                           np.logical_and(j_sup >= 0, j_sup < ncols))
                           )[0]
        #kludge to look at frames before and after predicted
        #tmp1=np.any(image_stack[indices[0]-1][i_sup[idx_mask], j_sup[idx_mask]])
        #tmp2=np.any(image_stack[indices[0]][i_sup[idx_mask], j_sup[idx_mask]])
        #tmp3=np.any(image_stack[indices[0]+1][i_sup[idx_mask], j_sup[idx_mask]])
            
        #tmp_confidence[iref] = np.any([tmp1,tmp2,tmp3])
        
        tmp_confidence[iref] = np.any(image_stack[indices[0]][i_sup[idx_mask], j_sup[idx_mask]])
        pass
    return sum(tmp_confidence)/float(len(tmp_confidence))

#
img_nums=np.arange(img_start,img_start+num_imgs,1)

#==============================================================================
# %% LOAD GRAIN DATA
#==============================================================================
ff_data=np.loadtxt(grain_out_file)

#ff_data=np.atleast_2d(ff_data[2,:])

exp_maps=ff_data[:,3:6]
t_vec_ds=ff_data[:,6:9]

n_grains=exp_maps.shape[0]

rMat_c = rot.rotMatOfExpMap(exp_maps.T)

#==============================================================================
# %% INSTRUMENT
#==============================================================================
# load config
instr_cfg = yaml.load(open(det_file, 'r'))

tiltAngles = instr_cfg['detector']['transform']['tilt_angles']
tVec_d = np.array(instr_cfg['detector']['transform']['t_vec_d']).reshape(3, 1)
chi = instr_cfg['oscillation_stage']['chi']
tVec_s = np.array(instr_cfg['oscillation_stage']['t_vec_s']).reshape(3, 1)

rMat_d = makeDetectorRotMat(tiltAngles)
rMat_s = makeOscillRotMat([chi, 0.])

pixel_size = instr_cfg['detector']['pixels']['size']

nrows = instr_cfg['detector']['pixels']['rows']
ncols = instr_cfg['detector']['pixels']['columns']

row_dim = pixel_size[0]*nrows # in mm 
col_dim = pixel_size[1]*ncols # in mm 

x_col_edges = pixel_size[1]*(np.arange(ncols+1) - 0.5*ncols)
y_row_edges = pixel_size[0]*(np.arange(nrows+1) - 0.5*nrows)[::-1]

panel_dims = [(-0.5*ncols*pixel_size[1],
               -0.5*nrows*pixel_size[0]),
              ( 0.5*ncols*pixel_size[1],
                0.5*nrows*pixel_size[0])]

# a bit overkill, but grab max two-theta from all pixel transforms
rx, ry = np.meshgrid(x_col_edges, y_row_edges)
gcrds = detectorXYToGvec(np.vstack([rx.flatten(), ry.flatten()]).T,
                         rMat_d, rMat_s,
                         tVec_d, tVec_s, np.zeros(3))
pixel_tth = gcrds[0][0]

detector_params = np.hstack([tiltAngles, tVec_d.flatten(), chi, tVec_s.flatten()])

distortion = None

panel_dims_expanded = [(-10, -10), (10, 10)]
panel_buffer = 0.02 # mm

# scan range and period
ome_period = (ome_period_deg[0]*np.pi/180.,ome_period_deg[1]*np.pi/180.)
ome_range = [(ome_range_deg[0][0]*np.pi/180.,ome_range_deg[0][1]*np.pi/180.)]
ome_step = ome_step_deg*np.pi/180.
nframes = 0
for i in range(len(ome_range)):
    del_ome = ome_range[i][1]-ome_range[i][0]
    nframes += int((ome_range[i][1]-ome_range[i][0])/ome_step)
    pass
ome_edges = np.arange(nframes+1)*ome_step

num_ori_grid_pts=1.

#==============================================================================
# %% near field map parameters
#==============================================================================
# form test grid and make main loop over spatial coordinates.  
cvec_s = 0.001*np.arange(h_bnds[0], h_bnds[1]+h_spacing,h_spacing)
#cvy = 0.001*np.arange(0, 1,1)
cvy = 0.001*np.arange(v_bnds[0], v_bnds[1]+v_spacing,v_spacing)

Xs, Ys, Zs = np.meshgrid(cvec_s, cvy, cvec_s)
test_crds = np.vstack([Xs.flatten(), Ys.flatten(), Zs.flatten()]).T
n_crds = len(test_crds)

# biggest projected diameter for 5 micron cube
max_diameter = np.sqrt(3)*0.002

row_dilation = np.round( 0.5*max_diameter/float(pixel_size[0]) )
col_dilation = np.round( 0.5*max_diameter/float(pixel_size[1]) )

i_dil, j_dil = np.meshgrid(np.arange(-row_dilation, row_dilation + 1),
                           np.arange(-col_dilation, col_dilation + 1))
i_dil = np.array([i_dil.flatten()], dtype=int)
j_dil = np.array([j_dil.flatten()], dtype=int)

#==============================================================================
# %% LOAD MATERIAL DATA
#==============================================================================
materials=cpl.load(open( mat_file_loc, "rb" ))

check=np.zeros(len(materials))
for ii in np.arange(len(materials)):
    check[ii]=materials[ii].name==mat_name

mat_used=materials[np.where(check)[0][0]]
#mat_used.atominfo=atom_info
#niti_mart.beamEnergy = valunits.valWUnit("wavelength","ENERGY",61.332,"keV")
mat_used.beamEnergy = valunits.valWUnit("wavelength","ENERGY",x_ray_energy,"keV")            
mat_used.planeData.exclusions = np.zeros(len(mat_used.planeData.exclusions), dtype=bool)
f=mat_used.planeData.calcStructFactor(mat_used.atominfo)  
weak_peaks=np.where(f<0.01*np.max(f))[0]

if max_tth>0.:
     mat_used.planeData.tThMax = np.amax(np.radians(max_tth))   
else:
    mat_used.planeData.tThMax = np.amax(pixel_tth)        

mat_used.planeData.set_exclusions(weak_peaks)
pd=mat_used.planeData

#==============================================================================
# %% MAKE MEDIAN DARK
#==============================================================================
dark_stack=np.zeros([num_for_median,2048,2048])

print('Loading data for median...')
for ii in np.arange(num_for_median):
    print('Image #: ' + str(ii))
    dark_stack[ii,:,:]=img.imread(data_folder+'nf_%0.5d.tif' %(img_nums[ii]))
    #image_stack[ii,:,:]=np.flipud(tmp_img>threshold)
print('making median...')
dark=np.median(dark_stack,axis=0)

#==============================================================================
# %% LOAD IMAGE DATA AND PROCESS
#==============================================================================
image_stack=np.zeros([img_nums.shape[0],2048,2048],dtype=bool)

print('Loading and Cleaning Images...')
for ii in np.arange(img_nums.shape[0]):
    print('Image #: ' + str(ii))
    tmp_img=img.imread(data_folder+'nf_%0.5d.tif' %(img_nums[ii]))-dark
    #image procesing
    image_stack[ii,:,:]=img.morphology.binary_erosion(tmp_img>threshold,iterations=num_erosions)
    image_stack[ii,:,:]=img.morphology.binary_dilation(image_stack[ii,:,:],iterations=num_dilations)
    
#%A final dilation that includes omega
print('Final Dilation Including Omega....')
image_stack=img.morphology.binary_dilation(image_stack,iterations=ome_dilation_iter)
    
#==============================================================================
# %% GENERATE DATA FOR INDEXING
#==============================================================================
# first evaluate diffraction angles from orientation list (fixed)
print('Generating gvecs...')
all_angles = []
for i in range(n_grains):
    gparams = np.hstack([exp_maps[i, :].flatten(), t_vec_ds[i].flatten(),ref_vparams])
    sim_results = simulateGVecs(pd,
                                detector_params,
                                gparams,
                                panel_dims=panel_dims,
                                pixel_pitch=pixel_size,
                                ome_range=ome_range,
                                ome_period=ome_period,
                                distortion=None)
    all_angles_tmp=sim_results[2]
    all_angles.append(all_angles_tmp)
    
    pass
print('done...')

input_p = np.hstack([
    np.tile(exp_maps, (n_crds, 1)),
    np.tile(test_crds, (1, n_grains)).reshape(n_crds*n_grains, 3), 
    np.tile(np.arange(n_grains),(1,n_crds)).T])

params = {
    'rMat_d':rMat_d,
    'chi':chi,
    'tVec_d':tVec_d,
    'tVec_s':tVec_s,
    'distortion':distortion,
    'panel_dims':panel_dims,
    'panel_buffer':panel_buffer,
    'x_col_edges':x_col_edges,
    'y_row_edges':y_row_edges,
    'ome_edges':ome_edges,
    'i_dil':i_dil,
    'j_dil':j_dil,
    'nrows':nrows,
    'ncols':ncols,
    'image_stack':image_stack,
    'pd':pd,
    'detector_params':detector_params,
    'pixel_size':pixel_size,
    'ome_range':ome_range,
    'ome_period':ome_period,
    'all_angles':all_angles
    }

#==============================================================================
# %% GENERATE DATA FOR INDEXING   ADDED FOR EASY CAL'S - AB
#==============================================================================
cvec_s = 0.001*np.arange(h_bnds[0], h_bnds[1]+h_spacing,h_spacing)
cvy = 0.001*np.arange(v_bnds[0], v_bnds[1]+v_spacing,v_spacing)
Xs, Ys, Zs = np.meshgrid(cvec_s, cvy, cvec_s)
test_crds = np.vstack([Xs.flatten(), Ys.flatten(), Zs.flatten()]).T
n_crds = len(test_crds)
max_diameter = np.sqrt(3)*0.002
row_dilation = np.round( 0.5*max_diameter/float(pixel_size[0]) )
col_dilation = np.round( 0.5*max_diameter/float(pixel_size[1]) )
i_dil, j_dil = np.meshgrid(np.arange(-row_dilation, row_dilation + 1),
                           np.arange(-col_dilation, col_dilation + 1))
i_dil = np.array([i_dil.flatten()], dtype=int)
j_dil = np.array([j_dil.flatten()], dtype=int)
print('Generating gvecs...')
all_angles = []
for i in range(n_grains):
    gparams = np.hstack([exp_maps[i, :].flatten(), t_vec_ds[i].flatten(),ref_vparams])
    sim_results = simulateGVecs(pd,
                                detector_params,
                                gparams,
                                panel_dims=panel_dims,
                                pixel_pitch=pixel_size,
                                ome_range=ome_range,
                                ome_period=ome_period,
                                distortion=None)
    all_angles_tmp=sim_results[2]
    all_angles.append(all_angles_tmp)
    pass
print('done...')
input_p = np.hstack([
    np.tile(exp_maps, (n_crds, 1)),
    np.tile(test_crds, (1, n_grains)).reshape(n_crds*n_grains, 3), 
    np.tile(np.arange(n_grains),(1,n_crds)).T])
params = {
    'rMat_d':rMat_d,
    'chi':chi,
    'tVec_d':tVec_d,
    'tVec_s':tVec_s,
    'distortion':distortion,
    'panel_dims':panel_dims,
    'panel_buffer':panel_buffer,
    'x_col_edges':x_col_edges,
    'y_row_edges':y_row_edges,
    'ome_edges':ome_edges,
    'i_dil':i_dil,
    'j_dil':j_dil,
    'nrows':nrows,
    'ncols':ncols,
    'image_stack':image_stack,
    'pd':pd,
    'detector_params':detector_params,
    'pixel_size':pixel_size,
    'ome_range':ome_range,
    'ome_period':ome_period,
    'all_angles':all_angles
    }
    

#==============================================================================
# %% CALIBRATION: PARAMETER STUDY
#==============================================================================
#note this will work best with a single layer

ytilt = np.radians(ytilt)
ztilt = np.radians(ztilt)

distance_array = np.random.uniform(distance[0], distance[1], number_of_tests)
x_cen_array    = np.random.uniform(x_cen[0], x_cen[1], number_of_tests)
ytilt_array    = np.random.uniform(ytilt[0], ytilt[1], number_of_tests)
ztilt_array    = np.random.uniform(ztilt[0], ztilt[1], number_of_tests)

# Create all permutations of the four parameters above that we wish to vary
# parameters = cartesian((distance_array, x_cen_array, ytilt_array, ztilt_array))
parameters = np.vstack((distance_array, x_cen_array, ytilt_array, ztilt_array))
parameters = parameters.T

Xs, Ys, Zs = np.meshgrid(cvec_s, cvy, cvec_s)
test_crds = np.vstack([Xs.flatten(), Ys.flatten(), Zs.flatten()]).T

num_parm_pts = len(parameters)

print "Shape of parameters array is:"
print np.shape(parameters)

trial_data=np.zeros([num_parm_pts,cvec_s.shape[0],cvec_s.shape[0]])

# Open a file to store parameter and confidence data
out_file = open(paramstudy_output_file, 'w+')
template_file = "{0:12.4f}{1:12.4f}{2:12.4f}{3:12.4f}{4:12.4f}{5:12.4f}\n"

for jj in np.arange(num_parm_pts):
    # Create a copy of tVec_d
    tmp_td = copy.copy(tVec_d)

    tmp_td[2]  = parameters[jj, 0]                                              # Set distance
    tmp_td[0]  = parameters[jj, 1]                                              # Set x_cen
    rMat_d_tmp = makeDetectorRotMat([0., parameters[jj, 2], parameters[jj, 3]]) # Set tilts

    print('Current parameter permutation:')
    print(parameters[jj, :])
    print('detector tilt orientation matrix:')        
    print(rMat_d_tmp)
    print('translation vector to detector centers:') 
    print(tmp_td)  

    params = {
        'rMat_d':rMat_d_tmp,
        'chi':chi,
        'tVec_d':tmp_td,
        'tVec_s':tVec_s,
        'distortion':distortion,
        'panel_dims':panel_dims,
        'panel_buffer':panel_buffer,
        'x_col_edges':x_col_edges,
        'y_row_edges':y_row_edges,
        'ome_edges':ome_edges,
        'i_dil':i_dil,
        'j_dil':j_dil,
        'nrows':nrows,
        'ncols':ncols,
        'image_stack':image_stack,
        'pd':pd,
        'detector_params':detector_params,
        'pixel_size':pixel_size,
        'ome_range':ome_range,
        'ome_period':ome_period,
        'all_angles':all_angles
        }

    confidence = None
    if ncpus > 1:
        # multiple process version

        try:
            pool = mp.Pool(ncpus, testorientation_init, (params, ))
            confidence = pool.map_async(testorientation_reduced, input_p,chunksize=chunksize)
            pool.close()
            while (True):
              if (confidence.ready()): break
              remaining = confidence._number_left
              print "Waiting for", remaining, "chunks remaining..."
              time.sleep(2.)
        except:
            print "Something went wrong for this parameter set. setting confidence to zero."
            confidence = None
    else:
        # single process version.
        global paramMP
        testorientation_init(params) # sets paramMP
        confidence = map(testorientation_reduced, input_p)
        paramMP = None # clear paramMP

    if confidence is not None:
        try:
            # Successfully calculated everything
            conf = np.r_[confidence.get()].reshape(n_crds, n_grains).T
            trial_data[jj]=np.max(conf,axis=0).reshape([cvec_s.shape[0],cvec_s.shape[0]])
        except:
            conf = np.array([0])

        # Write the result for this parameter set to the text file
        out_file.write(template_file.format(parameters[jj, 0],
                                            parameters[jj, 1], 
                                            parameters[jj, 2], 
                                            parameters[jj, 3], 
                                            np.amax(conf), (conf > confidence_thresh).sum()))
        out_file.flush()
    else:
        # Something went wrong
        out_file.write(template_file.format(parameters[jj, 0],
                                            parameters[jj, 1],
                                            parameters[jj, 2],
                                            parameters[jj, 3],
                                            0.0, 0.0))
        out_file.flush()

#
out_file.close()

