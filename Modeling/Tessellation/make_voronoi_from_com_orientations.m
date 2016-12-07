%% Make a voronoi tessellation from grain center of mass and orientation
% This will distribute NxNxN grid points in a cube shaped specimen of size
% LxLxL microns. Then create M grains in it using Voronoi
% tessellation. The output will be written to a text file with name
% ms_file_name ready to be sent to the program for simulating virtual
% diffraction patterns.
%
% PREREQUISITES:
%  MTEX and export_fig packages are required. 
%%
grid_type = 1;                                                             % 0 | 1: 0 = auto-generate uniform grid, 1 = read mesh/connectivity files and generate grid
L         = 1000;                                                          % Sample dimensions in microns (Not used if grid_type == 1)
N         = 100;                                                           % NxNxN uniform grid will be created from the specimen (Not used if grid_type == 1)
%
mesh_data_dir          = '/Users/Harshad/Documents/Work/Research_data/HEDM_Study_of_Deformation_in_SMA/CHESS_2015_June_NiTi_SC/midas_analysis/tessellation';
node_file_name         = 'nodes.inp';
connectivity_file_name = 'connectivity.inp';
%
M         = 3;                                                             % Number of grains
com       = [   95.0597  -23.3888  -60.1245;                               % Grain center of mass. Mx3 array. Leave empty [] to generate random com.
                97.1985   47.7753   35.8656;
              -254.0877  -64.5567  154.6093];
quat      = [    0.0172    0.6982    0.1323    0.7034;                     % Grain orientation in quaternions. Mx4 array. Leave empty to generate random orientations
                 0.0476    0.0555   -0.7723   -0.6311;
                 0.6612   -0.0786    0.0742    0.7424];
sigma_o   = 1.0E-2;                                                        % Amount of orientation spread (mosaicity) to add to each grain (degree)
sigma_e   = 5.0E-3;                                                        % Amount of deformation spread to add to each grain
ms_file_name = 'ms-synth-cubic.csv';                                       % Name of the file to which the output is written. This file will be the input for forward modeling code.
%%
if(~exist('com', 'var') || isempty(com))
    com    = (rand(M, 3) - 0.5) * L;  % Random grain centroid coordinates
end
%
if(~exist('quat', 'var') || isempty(quat))
    quat   = randq(M);                % Random quaternions for grain orientations
    quat   = squeeze(double(quat));   % Convert from quat MTEX object to normal matrix
end
% Create a tessellation
if(grid_type == 0)
    tess = zeros(N*N*N, 1);               % A 3D matrix containing grain ID for each point
    tess_coords = zeros(N*N*N, 3);        % A 3D matrix containing coordinates of each point
    for ii = 1:N
        for jj = 1:N
            for kk = 1:N
                tess_coords(N*N*(ii - 1) + N*(jj - 1) + kk, :) = [ii/N*L - L/2.0 jj/N*L - L/2.0 kk/N*L - L/2.0];
            end
        end
    end
else
    % User supplied grid
    nodes_in = importdata(fullfile(mesh_data_dir, node_file_name));
    connectivity_in = importdata(fullfile(mesh_data_dir, connectivity_file_name));
    nodes(nodes_in(:, 1), :) = nodes_in(:, 2:end);
    connectivity(connectivity_in(:, 1), :) = connectivity_in(:, 2:end);
    elements = zeros(size(connectivity, 1), 3);
    for ii = 1:size(elements, 1)
        elements(ii, :) = mean(nodes(connectivity(ii, :)', :), 1);
    end
    %
    tess_coords = elements;
    tess = zeros(size(tess_coords, 1), 1);
end

% Calculate distance between each material point and grain centroids
dist = pdist2(tess_coords, com);

for ii = 1:size(tess, 1)
    % Voronoi tessellation is where each material point is assigned
    % to a grain whose COM is closest to that material point
    [~, tess(ii)] = min(dist(ii, :));
    
end
%%
% File to write microstructural data
f = fopen(ms_file_name, 'w');
strain = rand(M, 6)/100.0;        % Small random strain at each grain
for ii = 1:size(tess, 1)
    % Write microstructure data in the column format that fwdmodel
    % python code reads.
    % Strain is really right stretch tensor.
    % I am adding a small noise to orientations and strains. Hence
    % the rand() calls.
    fprintf(f, '%12.4f, %12.4f, %12.4f, %s, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f\n', ...
        tess_coords(ii, 1), tess_coords(ii, 2), tess_coords(ii, 3), 'NiTi_cubic', ...
        quat(tess(ii), 1) + rand*sigma_o, quat(tess(ii), 2) + rand*sigma_o, quat(tess(ii), 3) + rand*sigma_o, quat(tess(ii), 4) + rand*sigma_o, ...
        1.0 + strain(tess(ii), 1) + rand*sigma_e, 1.0 + strain(tess(ii), 2) + rand*sigma_e, 1.0 + strain(tess(ii), 3) + rand*sigma_e, ...
        strain(tess(ii), 4) + rand*sigma_e, strain(tess(ii), 5) + rand*sigma_e, strain(tess(ii), 6) + rand*sigma_e);
end
% Don't forget to close files.
fclose(f);
%%
% Plot COM, orientations etc? Useful for debugging.
plot_figures = 1;
%
if(plot_figures)
    % IPF with orientations
    q = quaternion(quat');
    o = orientation(q, crystalSymmetry('cubic'), specimenSymmetry('1'));
    plotIPDF(o, yvector, 'xAxisDirection','east', 'MarkerSize', 8);
    export_fig 'NiTi-A_orientations' -png -r100
    % Scatter plot with tessellation
    figure;
    scatter3(tess_coords(:, 1), tess_coords(:, 2), tess_coords(:, 3), 4, tess(:), 'filled')
    axis equal; axis vis3d;
    export_fig 'NiTi-A_tess' -png -r100
end
