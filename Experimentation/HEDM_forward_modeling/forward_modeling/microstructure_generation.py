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


# Random cubic grains without mosaicity
def generate_cubic_grains_random_ideal(nipt=1000, output_file="ms-data-test.csv"):
    '''
    Generate a random microstructure dataset and write it to csv.
    Useful for testing the forward modeling code
    '''

    grid_size = np.floor(nipt**(1./3.))
    x = np.linspace(0, 1, grid_size)
    y = np.linspace(0, 1, grid_size)
    z = np.linspace(0, 1, grid_size)
    xmesh, ymesh, zmesh = np.meshgrid(x, y, z)

    template = "{0:12.4f},{1:12.4f},{2:12.4f},{3:>12s},{4:12.4f},{5:12.4f},{6:12.4f},{7:12.4f},{8:12.4f},{9:12.4f},{10:12.4f},{11:12.4f},{12:12.4f},{13:12.4f}\n"

    f = open(output_file, 'w+')

    for ii in range(len(x)):
        for jj in range(len(y)):
            for kk in range(len(z)):
                quat_random = np.random.rand(1, 4)
                quat_random = quat_random / np.linalg.norm(quat_random)
                def_grad_random = np.random.rand(1, 6)
                def_grad_random = def_grad_random / 1000.0

                f.write(template.format(xmesh[ii, jj, kk],                   # X coordinate
                                        ymesh[ii, jj, kk],                   # Y coordinate
                                        zmesh[ii, jj, kk],                   # Z coordinate
                                        "NiTi_cubic",                        # Phase name. Must correspond to a material in heXRD material file to be used.
                                        quat_random[0][0],                   # Orientation in quaternions (4 components)
                                        quat_random[0][1],
                                        quat_random[0][2],
                                        quat_random[0][3],
                                        (1. + def_grad_random[0][0]),        # Stretch tensor components (6. Three axial and three shear)
                                        (1. + def_grad_random[0][1]),
                                        (1. + def_grad_random[0][2]),
                                        def_grad_random[0][3],
                                        def_grad_random[0][4],
                                        def_grad_random[0][5]))
#--
# Near [1 0 0] cubic grain without mosaicity
def generate_cubic_grain_ideal(nipt=1000, output_file="ms-data-test.csv"):
    '''
    Generate a cubic single crystal microstructure dataset.
    2 mm cube sample. No mosaicity
    '''

    grid_size = np.floor(nipt**(1./3.))
    x = np.linspace(-1000, 1000, grid_size)
    y = np.linspace(-1000, 1000, grid_size)
    z = np.linspace(-1000, 1000, grid_size)
    xmesh, ymesh, zmesh = np.meshgrid(x, y, z)

    template = "{0:12.4f},{1:12.4f},{2:12.4f},{3:>12s},{4:12.4f},{5:12.4f},{6:12.4f},{7:12.4f},{8:12.4f},{9:12.4f},{10:12.4f},{11:12.4f},{12:12.4f},{13:12.4f}\n"

    f = open(output_file, 'w+')

    for ii in range(len(x)):
        for jj in range(len(y)):
            for kk in range(len(z)):
                quat_random = np.array([9.32918129e-01,        1.13125807e-01,        -3.32369294e-01,       7.99810469e-02])
                quat_random = quat_random / np.linalg.norm(quat_random)

                def_grad_random = np.random.rand(1, 6)
                def_grad_random = def_grad_random / 100.0

                f.write(template.format(xmesh[ii, jj, kk],
                                        ymesh[ii, jj, kk],
                                        zmesh[ii, jj, kk],
                                        "NiTi_cubic",
                                        quat_random[0],
                                        quat_random[1],
                                        quat_random[2],
                                        quat_random[3],
                                        (1. + def_grad_random[0][0]),
                                        (1. + def_grad_random[0][1]),
                                        (1. + def_grad_random[0][2]),
                                        def_grad_random[0][3],
                                        def_grad_random[0][4],
                                        def_grad_random[0][5]))
#--
# Near [1 0 0] cubic grain with mosaicity
def generate_cubic_grain_mosaicity(nipt=1000, output_file="ms-data-test.csv"):
    '''
    Generate a cubic single crystal microstructure dataset.
    2 mm cube sample. No mosaicity
    '''

    grid_size = np.floor(nipt**(1./3.))
    x = np.linspace(-1000, 1000, grid_size)
    y = np.linspace(-1000, 1000, grid_size)
    z = np.linspace(-1000, 1000, grid_size)
    xmesh, ymesh, zmesh = np.meshgrid(x, y, z)

    template = "{0:12.4f},{1:12.4f},{2:12.4f},{3:>12s},{4:12.4f},{5:12.4f},{6:12.4f},{7:12.4f},{8:12.4f},{9:12.4f},{10:12.4f},{11:12.4f},{12:12.4f},{13:12.4f}\n"

    f = open(output_file, 'w+')

    for ii in range(len(x)):
        for jj in range(len(y)):
            for kk in range(len(z)):
                quat_random = np.array([9.32918129e-01,        1.13125807e-01,        -3.32369294e-01,       7.99810469e-02])
                quat_dev = np.random.rand(1, 4) / 100.0
                quat_random[0] = quat_random[0] + quat_dev[0][0]
                quat_random[1] = quat_random[1] + quat_dev[0][1]
                quat_random[2] = quat_random[2] + quat_dev[0][2]
                quat_random[3] = quat_random[3] + quat_dev[0][3]
                quat_random = quat_random / np.linalg.norm(quat_random)

                def_grad_random = np.random.rand(1, 6)
                def_grad_random = def_grad_random / 100.0

                f.write(template.format(xmesh[ii, jj, kk],
                                        ymesh[ii, jj, kk],
                                        zmesh[ii, jj, kk],
                                        "NiTi_cubic",
                                        quat_random[0],
                                        quat_random[1],
                                        quat_random[2],
                                        quat_random[3],
                                        (1. + def_grad_random[0][0]),
                                        (1. + def_grad_random[0][1]),
                                        (1. + def_grad_random[0][2]),
                                        def_grad_random[0][3],
                                        def_grad_random[0][4],
                                        def_grad_random[0][5]))

# Near [1 0 0] monoclinic grain with mosaicity
def generate_mono_grain_mosaicity(nipt=1000, output_file="ms-data-test.csv", material_name='NiTi_mono', mosaicity=None, defgrad_spread=None):
    '''
    Generate a cubic single crystal microstructure dataset.
    2 mm cube sample. No mosaicity
    '''

    grid_size = np.floor(nipt**(1./3.))
    x = np.linspace(-1000, 1000, grid_size)
    y = np.linspace(-1000, 1000, grid_size)
    z = np.linspace(-1000, 1000, grid_size)
    xmesh, ymesh, zmesh = np.meshgrid(x, y, z)

    template = "{0:12.4f},{1:12.4f},{2:12.4f},{3:>12s},{4:12.4f},{5:12.4f},{6:12.4f},{7:12.4f},{8:12.4f},{9:12.4f},{10:12.4f},{11:12.4f},{12:12.4f},{13:12.4f}\n"

    f = open(output_file, 'w+')

    for ii in range(len(x)):
        for jj in range(len(y)):
            for kk in range(len(z)):
                quat_random = np.array([   -0.0041,   -0.8192,   -0.3506,    0.4538])
		if mosaicity is None:
		    mosaicity = 0.001
                quat_dev = np.random.rand(1, 4) * mosaicity
                quat_random[0] = quat_random[0] + quat_dev[0][0]
                quat_random[1] = quat_random[1] + quat_dev[0][1]
                quat_random[2] = quat_random[2] + quat_dev[0][2]
                quat_random[3] = quat_random[3] + quat_dev[0][3]
                quat_random = quat_random / np.linalg.norm(quat_random)

		if defgrad_spread is None:
		    defgrad_spread = 0.001
                def_grad_random = np.random.rand(1, 6)
                def_grad_random = def_grad_random * defgrad_spread

                f.write(template.format(xmesh[ii, jj, kk],
                                        ymesh[ii, jj, kk],
                                        zmesh[ii, jj, kk],
                                        material_name.strip(),
                                        quat_random[0],
                                        quat_random[1],
                                        quat_random[2],
                                        quat_random[3],
                                        (1. + def_grad_random[0][0]),
                                        (1. + def_grad_random[0][1]),
                                        (1. + def_grad_random[0][2]),
                                        def_grad_random[0][3],
                                        def_grad_random[0][4],
                                        def_grad_random[0][5]))
